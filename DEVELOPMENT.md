# Development

Contributor and maintainer notes for `phyloacc-workflows` — repo layout, the test suite, and
the release process. For installing and running the pipeline, see the [README](README.md).

## Repo layout

- `Snakefile` - entry point; includes the workflow files below based on config toggles
- `workflow/` - the active Snakemake rule files (`phylofit_models.smk`, `phylop_regions.smk`, `phastcons_cnees.smk`, `cnees.smk` - the shared, source-agnostic CNEE-building stage that turns conserved elements from phastCons and/or phyloP into CNEEs - and `align_elements.smk` - the align-only mode); `workflow/legacy/` holds earlier versions of the pipeline not used by the current `Snakefile`
- `lib/` - shared Python helpers used across the workflow files
- `utils/` - standalone scripts/tools invoked by rules
- `envs/environment.yml` - pinned dependencies for the `phyloacc_workflows` environment
- `cfgs/` - example filled-in configs from past analyses

## Running tests

`tests/` holds three tiers of tests. These are a developer/maintainer concern, not
something needed to run the pipeline itself, so `pytest` is not part of
`envs/environment.yml`. To run them, install pytest into your existing environment and
run it from the repo root:

```bash
mamba install -n phyloacc-workflows pytest    # or: pip install pytest
pytest tests/
```

- **Unit tests** (`test_intervals.py`, `test_parsing.py`) - pure-Python logic in `lib/`
  (interval merging, BED/GFF parsing, Newick tip extraction, etc. - the same code the real
  Snakemake rules import and run). No external tools needed. ~1-2s.
- **DAG/config validation** (`test_dag_validation.py`) - each test invokes a real
  `snakemake -n` subprocess against a minimal config to check that bad config values are
  rejected with the expected error. Needs `snakemake` on `PATH`. ~3-4 min.
- **Integration** (`tests/integration/`) - runs the real pipeline (real
  `mafutils`/`phyloFit`/`phastCons`, not mocked) end-to-end against a small real-data
  fixture; see `tests/integration/README.md` for the fixture's provenance. Needs
  `mafutils`/`phyloFit`/`phastCons` actually on `PATH` - skips cleanly (not a failure) if
  they're missing. ~4-5 min. Includes a soft "silver standard" check that warns (does not
  fail) if a run's summary numbers drift far from a committed reference - see that
  directory's README for how to regenerate it if the fixture or tool versions change.

Running everything together takes roughly 7-8 minutes, dominated by the last two tiers'
real subprocess/tool invocations.

## Releasing

Version bumps are manual, not automated. To cut a release:

1. On a branch, bump `version` (and the `releasedate-*`/`latest-commit-date`
   fields) in `lib/info.yaml`, and `version`/`date-released` in
   `CITATION.cff`. Open a PR and merge it to `main`.
2. From `main`, tag and push:
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
3. `.github/workflows/release.yml` creates the GitHub Release from that tag
   automatically - it does not check that the tag matches `lib/info.yaml`,
   so a mismatched tag publishes a Release with the wrong version recorded
   in-repo without any warning.

Activate the repo's pre-push hook once per clone to guard against that:

```bash
git config core.hooksPath .githooks
```

With it active, pushing a `vX.Y.Z` tag whose `lib/info.yaml` version doesn't
match is blocked locally with a clear error, before the push happens - this
is the only check in place, so activating it is strongly recommended. It's
optional, dev-only tooling (unrelated to `phyloacc_workflows setup`, which
sets up the environment for *running* the pipeline).
