<!-- Barest-minimum detectability blurb for the tutorial. Equations are written for MathJax
     ( \[ ... \] ); if the tutorial doesn't use MathJax, render them as images instead.
     Link target assumes the biologist page is reachable at the given path. -->

### Can per-site phyloP call conserved sites on your tree?

Per-site phyloP tests each base for conservation, producing a p-value:

\[ p = \operatorname{erfc}\!\left(\sqrt{\mathrm{LRT}/2}\right) \]

where \(\mathrm{LRT}\) measures how much better a slower-than-neutral (conserved) model fits the base
than the neutral model — a larger \(\mathrm{LRT}\) gives a smaller \(p\), i.e. stronger evidence of
conservation — and \(\operatorname{erfc}\) is the complementary error function. phyloP reports these
as scores, \(-\log_{10} p\).

These scores have a hard **ceiling**, which is determined mainly by the total neutral tree depth
(substitutions/site) and the number of species in the tree, with a smaller contribution from the
substitution rate itself (a model whose fastest base turns over more quickly than average raises the
ceiling). Concretely, the largest \(\mathrm{LRT}\) a site can reach — which sets that ceiling — is
capped by whichever of two limits is smaller:

\[ \mathrm{LRT}_{\max} \approx \min\!\Big(\underbrace{2\,T\,\max_b(-Q_{bb})}_{\text{depth}\,\times\,\text{rate}},\ \ \underbrace{2(n-1)\big(-\ln \min_b \pi_b\big)}_{\text{species}\,\times\,\text{composition}}\Big) \]

The first term grows with the tree depth \(T\) and the fastest base's exit rate \(\max_b(-Q_{bb})\); the
second with the number of species \(n\) and the rarest base's frequency \(\min_b\pi_b\). The ceiling
follows whichever term is smaller — so too little depth, or too few species, each cap it.

Since millions of sites are being tested, we must correct the resulting scores for multiple tests.
This results in a multiple-testing threshold that is set by the number of sites being tested and the
desired false-positive rate:

\[ \text{threshold} = \log_{10}\!\left(M/\alpha\right) \]

where \(M\) is the number of sites tested and \(\alpha\) the desired false-positive rate — the more
sites tested, the higher the bar.

If phyloP's score ceiling falls below the multiple-testing threshold, no site is detectable,
regardless of how conserved it is. In other words, conserved sites are detectable only when:

\[ -\log_{10}\operatorname{erfc}\!\left(\sqrt{\mathrm{LRT}_{\max}/2}\right)\;\ge\;\log_{10}(M/\alpha) \]

This is the condition under which the tree's best-possible score clears the correction for \(M\)
tested sites at level \(\alpha\). The ceiling rises with tree depth but eventually saturates, so two
things can leave phyloP powerless: **too little total tree depth**, or **too few species**. Deep,
taxon-rich trees are comfortably detectable; shallow or few-taxon trees need element-based
conservation instead.

[Will per-site phyloP work on your tree? →](detectability-biologist.html)
