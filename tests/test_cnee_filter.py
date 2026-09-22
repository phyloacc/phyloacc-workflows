#############################################################################
# Unit tests for lib/cnee_filter.py - the clean single-copy ortholog filtering
# policy applied to per-CNEE alignment FASTAs. Pure logic + the two mafutils TSV
# parsers + the FASTA masking rewrite. Run with `pytest tests/` from the repo root.
#############################################################################

import os

import pytest

import lib.cnee_filter as CF

REF = "galGal6"

# Thresholds mirroring the config defaults: split kept if gap <= 50; element kept
# if >= 4 real species remain; shared boundary if >=50% of split species (and >=3
# absolute) split at one boundary with gap spread <= 10.
THR = CF.Thresholds(
    split_max_gap_bp=50,
    min_species=4,
    shared_boundary_min_frac=0.5,
    shared_boundary_max_gap_spread_bp=10,
)


#############################################################################
# parsing helpers

def _write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def test_load_summary_by_name_not_position(tmp_path):
    p = tmp_path / "maf_fetch_summary.tsv"
    # extra columns + non-default order; parser must key by header name.
    _write(p,
           ["scaffold", "start", "end", "basename", "n.species",
            "n.split", "split.max.boundary.n.species", "split.max.boundary.gap.spread"],
           [["chr1", 100, 200, "chr1.cnee1", 45, 0, ".", "."],
            ["chr1", 300, 400, "chr1.cnee2", 44, 6, 6, 1]])
    s = CF.load_summary(str(p))
    assert set(s) == {"chr1.cnee1", "chr1.cnee2"}
    assert s["chr1.cnee1"].n_split == 0
    assert s["chr1.cnee1"].boundary_n_species is None      # '.' -> None
    assert s["chr1.cnee2"].n_split == 6
    assert s["chr1.cnee2"].boundary_n_species == 6
    assert s["chr1.cnee2"].boundary_gap_spread == 1


def test_load_loci_first_row_authoritative(tmp_path):
    p = tmp_path / "maf_fetch_loci.tsv"
    _write(p,
           ["region_id", "species", "src_scaffold", "src_size", "class", "max_gap"],
           [["chr1.cnee1", "galGal6", "chr1", 197608386, "single", 0],
            ["chr1.cnee1", "HLcolLiv2", "AKCR02000132", 7331655, "split", 111282],
            ["chr1.cnee1", "HLcolLiv2", "AKCR02000132", 7331655, "split", 111282]])  # 2nd block row
    loci = CF.load_loci(str(p))
    assert set(loci["chr1.cnee1"]) == {"galGal6", "HLcolLiv2"}
    pig = loci["chr1.cnee1"]["HLcolLiv2"]
    assert pig.cls == "split" and pig.max_gap == 111282
    assert pig.src_scaffold == "AKCR02000132" and pig.src_size == 7331655


def test_load_raises_on_missing_column(tmp_path):
    p = tmp_path / "bad_summary.tsv"
    _write(p, ["basename", "n.split"], [["chr1.cnee1", 0]])  # missing boundary cols
    with pytest.raises(ValueError, match="missing expected column"):
        CF.load_summary(str(p))


#############################################################################
# classify_species

def _info(cls, max_gap=0):
    return CF.LociInfo(species="X", cls=cls, max_gap=max_gap, src_scaffold=None, src_size=None)


@pytest.mark.parametrize("cls", ["single", "contiguous"])
def test_single_contiguous_kept(cls):
    assert CF.classify_species(_info(cls), THR, shared_boundary=False)[0] == "keep"


@pytest.mark.parametrize("cls,reason", [
    ("no_bases", CF.REASON_NO_BASES),
    ("multi_scaffold", CF.REASON_MULTI_SCAFFOLD),
    ("multi_strand", CF.REASON_MULTI_STRAND),
])
def test_always_masked_classes(cls, reason):
    decision, r = CF.classify_species(_info(cls, max_gap=None), THR, shared_boundary=False)
    assert decision == "mask" and r == reason


def test_short_split_kept_long_split_masked():
    assert CF.classify_species(_info("split", 50), THR, shared_boundary=False)[0] == "keep"   # boundary value
    assert CF.classify_species(_info("split", 9), THR, shared_boundary=False)[0] == "keep"
    d, reason = CF.classify_species(_info("split", 51), THR, shared_boundary=False)
    assert d == "mask" and reason == CF.REASON_SPLIT_LONG


def test_split_kept_when_shared_boundary_regardless_of_gap():
    # a huge gap that would normally be masked is kept under a shared-boundary element
    assert CF.classify_species(_info("split", 111282), THR, shared_boundary=True)[0] == "keep"


def test_unknown_class_masked():
    d, reason = CF.classify_species(_info("weird"), THR, shared_boundary=False)
    assert d == "mask" and reason == CF.REASON_UNKNOWN


#############################################################################
# is_shared_boundary

def _summ(n_split, b_n, spread):
    return CF.SummaryRow("r", n_split=n_split, boundary_n_species=b_n, boundary_gap_spread=spread)


def test_shared_boundary_true_when_concentrated():
    # 6 species split, all 6 at one boundary, gaps agree within 1bp (the mafutils example)
    assert CF.is_shared_boundary(_summ(6, 6, 1), THR) is True


def test_shared_boundary_false_when_scattered():
    # 6 species split but only 2 at the top boundary -> lineage-specific, not shared
    assert CF.is_shared_boundary(_summ(6, 2, 1), THR) is False


def test_shared_boundary_false_when_gaps_disagree():
    assert CF.is_shared_boundary(_summ(6, 6, 500), THR) is False


def test_shared_boundary_false_below_absolute_floor():
    # 2 species agreeing is not evidence of a shared event
    assert CF.is_shared_boundary(_summ(2, 2, 0), THR) is False


def test_shared_boundary_false_when_no_summary():
    assert CF.is_shared_boundary(None, THR) is False


#############################################################################
# decide_element

def _loci(**species_cls):
    # species_cls: name -> (class, max_gap)
    return {sp: CF.LociInfo(sp, cls, gap, None, None) for sp, (cls, gap) in species_cls.items()}


def test_element_all_clean_retained_unmasked():
    loci = {"r": _loci(galGal6=("single", 0), A=("single", 0), B=("contiguous", 0),
                       C=("single", 0), D=("single", 0))}
    d = CF.decide_element("r", {}, loci, REF, THR)
    assert d.masked_species == frozenset()
    assert d.n_kept_real == 5 and d.retained is True and d.shared_boundary is False


def test_element_masks_long_split_but_retains():
    loci = {"r": _loci(galGal6=("single", 0), A=("single", 0), B=("single", 0),
                       C=("single", 0), D=("split", 5000))}
    d = CF.decide_element("r", {}, loci, REF, THR)
    assert d.masked_species == frozenset({"D"})
    assert d.mask_reasons[CF.REASON_SPLIT_LONG] == 1
    assert d.n_kept_real == 4 and d.retained is True


def test_element_rejected_when_too_few_real_remain():
    # only 4 aligned species, one masked -> 3 remain < min_species(4)
    loci = {"r": _loci(galGal6=("single", 0), A=("single", 0), B=("single", 0),
                       C=("multi_scaffold", None))}
    d = CF.decide_element("r", {}, loci, REF, THR)
    assert d.masked_species == frozenset({"C"})
    assert d.n_kept_real == 3 and d.retained is False


def test_reference_never_masked_even_if_odd_class():
    loci = {"r": _loci(galGal6=("split", 99999), A=("single", 0), B=("single", 0),
                       C=("single", 0), D=("single", 0))}
    d = CF.decide_element("r", {}, loci, REF, THR)
    assert REF not in d.masked_species and d.retained is True


def test_element_rejected_when_reference_absent():
    loci = {"r": _loci(A=("single", 0), B=("single", 0), C=("single", 0),
                       D=("single", 0), E=("single", 0))}
    d = CF.decide_element("r", {}, loci, REF, THR)
    assert d.ref_present is False and d.retained is False


def test_shared_boundary_keeps_all_split_rows_and_flags():
    # 5 non-ref species all split with big gaps, all at one boundary -> keep all, flag
    loci = {"r": _loci(galGal6=("single", 0), A=("split", 300), B=("split", 301),
                       C=("split", 300), D=("split", 302), E=("split", 300))}
    summ = {"r": _summ(5, 5, 2)}
    d = CF.decide_element("r", summ, loci, REF, THR)
    assert d.shared_boundary is True
    assert d.masked_species == frozenset()      # nothing masked despite big gaps
    assert d.retained is True


def test_shared_boundary_still_masks_multi_scaffold():
    loci = {"r": _loci(galGal6=("single", 0), A=("split", 300), B=("split", 301),
                       C=("split", 300), D=("multi_scaffold", None), E=("split", 300))}
    summ = {"r": _summ(4, 4, 1)}
    d = CF.decide_element("r", summ, loci, REF, THR)
    assert d.shared_boundary is True
    assert d.masked_species == frozenset({"D"})   # chimera masked even in a shared element


#############################################################################
# mask_fasta_records

def _fasta(path, records):
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


def test_mask_replaces_with_Ns_same_length_in_place(tmp_path):
    p = tmp_path / "chr1.cnee1.fa"
    _fasta(p, [("galGal6.chr1:1-8(+) id:chr1.cnee1", "ACGTACGT"),
               ("HLcolLiv2.AKCR02000132:1-8(-) id:chr1.cnee1", "AC-GTACG")])
    n = CF.mask_fasta_records(str(p), str(p), {"HLcolLiv2"})
    assert n == 1
    seqs = {}
    for line in open(p):
        if line.startswith(">"):
            cur = CF.species_from_header(line.rstrip())
        else:
            seqs[cur] = line.strip()
    assert seqs["galGal6"] == "ACGTACGT"       # reference untouched
    assert seqs["HLcolLiv2"] == "NNNNNNNN"     # same length (incl. former gap column)


def test_mask_noop_copies_when_empty(tmp_path):
    src = tmp_path / "in.fa"
    dst = tmp_path / "out.fa"
    _fasta(src, [("galGal6", "ACGT"), ("A", "ACGT")])
    n = CF.mask_fasta_records(str(src), str(dst), set())
    assert n == 0
    assert open(dst).read() == open(src).read()


def test_mask_handles_multiline_sequences(tmp_path):
    p = tmp_path / "c.fa"
    with open(p, "w") as fh:
        fh.write(">A x\nACGT\nACGT\n>B\nTTTT\n")   # A spans two seq lines
    n = CF.mask_fasta_records(str(p), str(p), {"A"})
    assert n == 1
    seqs = {}
    for line in open(p):
        if line.startswith(">"):
            cur = CF.species_from_header(line.rstrip())
        else:
            seqs[cur] = line.strip()
    assert seqs["A"] == "N" * 8      # both lines collapsed, 8 residues -> 8 Ns
    assert seqs["B"] == "TTTT"


#############################################################################
# stow_sidecars

def test_stow_sidecars_moves_to_run_info(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "chr1.cnee1.fa").write_text(">A\nACGT\n")          # real output, must stay
    (outdir / "manifest.txt").write_text("chr1.cnee1.fa\n")       # must stay (.txt)
    (outdir / "maf_fetch_summary.tsv").write_text("x\n")
    (outdir / "maf_fetch_loci.tsv").write_text("y\n")
    (outdir / "maf_fetch.log").write_text("z\n")
    handled = CF.stow_sidecars(str(outdir), keep=True)
    assert set(handled) == {"maf_fetch_summary.tsv", "maf_fetch_loci.tsv", "maf_fetch.log"}
    # sidecars moved into run-info/, fasta dir left clean
    assert sorted(os.path.basename(p) for p in (outdir / "run-info").iterdir()) == \
        ["maf_fetch.log", "maf_fetch_loci.tsv", "maf_fetch_summary.tsv"]
    remaining = sorted(p.name for p in outdir.iterdir() if p.is_file())
    assert remaining == ["chr1.cnee1.fa", "manifest.txt"]     # only real outputs remain


def test_stow_sidecars_to_summary_dir_with_prefix(tmp_path):
    # Pipeline usage: move sidecars into a separate summary dir, chromosome-prefixed so
    # per-chrom sidecars sharing one dir don't collide.
    outdir = tmp_path / "fasta" / "1"
    outdir.mkdir(parents=True)
    summary = tmp_path / "summary" / "autosomes"
    summary.mkdir(parents=True)
    (summary / "1.cnee-ortho-filter.tsv").write_text("existing\n")   # pre-existing summary
    (outdir / "chr1.cnee1.fa").write_text(">A\nACGT\n")
    (outdir / "maf_fetch_summary.tsv").write_text("x\n")
    (outdir / "maf_fetch_loci.tsv").write_text("y\n")
    handled = CF.stow_sidecars(str(outdir), keep=True, dest_dir=str(summary), prefix="1.")
    assert set(handled) == {"1.maf_fetch_summary.tsv", "1.maf_fetch_loci.tsv"}
    assert (summary / "1.maf_fetch_summary.tsv").exists()
    assert (summary / "1.cnee-ortho-filter.tsv").exists()            # not clobbered
    assert sorted(p.name for p in outdir.iterdir()) == ["chr1.cnee1.fa"]   # fasta dir clean


def test_stow_sidecars_delete_when_not_kept(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "chr1.cnee1.fa").write_text(">A\nACGT\n")
    (outdir / "maf_fetch_summary.tsv").write_text("x\n")
    (outdir / "maf_fetch.log").write_text("z\n")
    handled = CF.stow_sidecars(str(outdir), keep=False)
    assert set(handled) == {"maf_fetch_summary.tsv", "maf_fetch.log"}
    assert not (outdir / "run-info").exists()
    assert sorted(p.name for p in outdir.iterdir()) == ["chr1.cnee1.fa"]


def test_stow_sidecars_noop_when_none_present(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "chr1.cnee1.fa").write_text(">A\nACGT\n")
    assert CF.stow_sidecars(str(outdir), keep=True) == []
    assert not (outdir / "run-info").exists()      # no empty run-info created


#############################################################################
# species_from_header

@pytest.mark.parametrize("header,expected", [
    (">galGal6.chr1:97715-97807(+) id:chr1.cnee1", "galGal6"),
    (">HLcolLiv2.AKCR02000033:222638-222737(-) id:chr1.cnee1", "HLcolLiv2"),
    (">HLcolLiv2:222638-222737(-)", "HLcolLiv2"),   # no scaffold (older header)
    (">galGal6", "galGal6"),                          # species-only
    ("galGal6.chr1:1-2(+)", "galGal6"),               # no leading '>'
])
def test_species_from_header(header, expected):
    assert CF.species_from_header(header) == expected


#############################################################################
# filter_cnee_directory - end-to-end over synthesized mafutils outputs
# (fasta files + loci/summary TSVs), i.e. everything cnee_alignments_chr does
# after the `mafutils fetch` call.

SPECIES = ["galGal6", "A", "B", "C", "D", "E"]


def _fasta_dir_element(outdir, name):
    with open(outdir / f"{name}.fa", "w") as fh:
        for sp in SPECIES:
            fh.write(f">{sp}.scaf:1-10(+) id:{name}\nACGTACGTAC\n")


def _loci_rows(loci_map):
    # loci_map: region -> list of (species, class, max_gap)
    rows = []
    for rid, specs in loci_map.items():
        for sp, cls, gap in specs:
            rows.append([rid, sp, "scaf", 1000, cls, gap])
    return rows


def test_filter_cnee_directory_end_to_end(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    rejects = tmp_path / "rej"
    for n in ["chr1.cnee1", "chr1.cnee2", "chr1.cnee3", "chr1.cnee4"]:
        _fasta_dir_element(outdir, n)

    # cnee1: all single           -> retained, clean, unchanged
    # cnee2: D long split (5000)  -> D blanked to N, retained (5 real >= 4)
    # cnee3: only 4 aligned, C multi_scaffold -> 3 real < 4 -> rejected
    # cnee4: A-E split big gaps at one boundary -> shared, kept intact + flagged
    loci_map = {
        "chr1.cnee1": [(sp, "single", 0) for sp in SPECIES],
        "chr1.cnee2": [("galGal6", "single", 0), ("A", "single", 0), ("B", "single", 0),
                       ("C", "single", 0), ("D", "split", 5000), ("E", "single", 0)],
        "chr1.cnee3": [("galGal6", "single", 0), ("A", "single", 0), ("B", "single", 0),
                       ("C", "multi_scaffold", ".")],
        "chr1.cnee4": [("galGal6", "single", 0), ("A", "split", 300), ("B", "split", 301),
                       ("C", "split", 300), ("D", "split", 302), ("E", "split", 300)],
    }
    lp = tmp_path / "maf_fetch_loci.tsv"
    _write(lp, ["region_id", "species", "src_scaffold", "src_size", "class", "max_gap"],
           _loci_rows(loci_map))
    sp_ = tmp_path / "maf_fetch_summary.tsv"
    _write(sp_, ["basename", "n.split", "split.max.boundary.n.species", "split.max.boundary.gap.spread"],
           [["chr1.cnee1", 0, ".", "."],
            ["chr1.cnee2", 1, 1, 0],
            ["chr1.cnee3", 0, ".", "."],
            ["chr1.cnee4", 5, 5, 2]])

    summary = CF.load_summary(str(sp_))
    loci = CF.load_loci(str(lp))
    outs = sorted(str(outdir / f"{n}.fa") for n in
                  ["chr1.cnee1", "chr1.cnee2", "chr1.cnee3", "chr1.cnee4"])
    retained, rejected, stats = CF.filter_cnee_directory(outs, str(rejects), summary, loci, REF, THR)

    assert set(retained) == {"chr1.cnee1.fa", "chr1.cnee2.fa", "chr1.cnee4.fa"}
    assert rejected == ["chr1.cnee3.fa"]
    # rejected file moved out of the main folder
    assert not (outdir / "chr1.cnee3.fa").exists()
    assert (rejects / "chr1.cnee3.fa").exists()

    def seqs(path):
        out = {}
        for line in open(path):
            if line.startswith(">"):
                cur = CF.species_from_header(line.rstrip())
            else:
                out[cur] = line.strip()
        return out

    # cnee1 untouched
    assert seqs(outdir / "chr1.cnee1.fa")["D"] == "ACGTACGTAC"
    # cnee2: D blanked, others intact
    s2 = seqs(outdir / "chr1.cnee2.fa")
    assert s2["D"] == "NNNNNNNNNN" and s2["A"] == "ACGTACGTAC" and s2["galGal6"] == "ACGTACGTAC"
    # cnee4: shared boundary -> nothing blanked despite big split gaps
    s4 = seqs(outdir / "chr1.cnee4.fa")
    assert all(s4[sp] == "ACGTACGTAC" for sp in SPECIES)

    assert stats["cnee_fetched"] == 4
    assert stats["cnee_retained"] == 3 and stats["cnee_rejected"] == 1
    assert stats["cnee_elements_all_clean"] == 1          # cnee1
    assert stats["cnee_elements_masked"] == 1             # cnee2
    assert stats["cnee_elements_shared_boundary"] == 1    # cnee4
    assert stats["cnee_rows_masked_split_long"] == 1      # cnee2 D
    assert stats["cnee_rows_masked_multi_scaffold"] == 0  # cnee3 rejected wholesale, not counted
    assert stats["cnee_rejected_min_species"] == 1        # cnee3


def test_filter_cnee_directory_rejects_when_reference_absent(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    with open(outdir / "chr1.cnee9.fa", "w") as fh:
        for sp in ["A", "B", "C", "D", "E"]:   # no reference
            fh.write(f">{sp}.scaf:1-10(+) id:chr1.cnee9\nACGTACGTAC\n")
    loci = {"chr1.cnee9": {sp: CF.LociInfo(sp, "single", 0, None, None)
                           for sp in ["A", "B", "C", "D", "E"]}}
    retained, rejected, stats = CF.filter_cnee_directory(
        [str(outdir / "chr1.cnee9.fa")], str(tmp_path / "rej"), {}, loci, REF, THR)
    assert rejected == ["chr1.cnee9.fa"]
    assert stats["cnee_rejected_ref_absent"] == 1
