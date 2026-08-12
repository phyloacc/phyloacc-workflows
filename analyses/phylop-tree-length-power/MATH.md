# The statistical-power ceiling of per-site phyloP — the math

Canonical derivation behind the detectability figures (`06_generic_detectability.py`,
`07_theoretical_detectability.py`) and `lib/phylop_power.py`. Everything here is
first-principles; the only empirical numbers quoted are the range of the model constant
κ for real fitted models (~1.3–1.5), and those can be dropped for a pure-theory version.

---

## 1. The detectability condition

A single genomic site is detectable by per-site phyloP under multiple-testing correction
iff the **maximum score it could possibly attain** clears the significance bar:

$$\mathrm{ceiling}(T,\text{model},n) \;\ge\; \log_{10}\!\big(M/\alpha\big)$$

- **T** — total neutral tree depth (sum of branch lengths, expected substitutions/site)
- **n** — number of taxa; **M** — number of sites tested; **α** — significance/FDR level

## 2. The bar (right-hand side): multiple testing

With M sites tested at level α, the single most-significant site must reach $p\le\alpha/M$
(the Benjamini–Hochberg rank-1 condition; equivalently Bonferroni for the top test):

$$\mathrm{bar}(M,\alpha) = -\log_{10}(\alpha/M) = \log_{10}(M/\alpha)$$

This side is exact and model-free.

## 3. The ceiling (left-hand side): exact form

phyloP scores a site as $-\log_{10}p$ from a likelihood-ratio test of a branch-scaling
parameter ρ (ρ=1 neutral, ρ<1 conserved). The strongest possible evidence is a **perfectly
conserved column** — every species the same base $b$, MLE $\hat\rho\to0$. For that column,

$$\mathrm{LRT}_b = 2\big[\log\pi_b - \log P(\text{all leaves}=b\mid\text{neutral})\big],\qquad
\mathrm{LRT}_{\max}=\max_b \mathrm{LRT}_b$$

where $P(\text{all leaves}=b\mid\text{neutral})$ is computed by **Felsenstein's pruning
algorithm** on the fitted tree (branch lengths uniformly scaled so the total length is $T$)
with the fitted rate matrix $Q$, transition probabilities $P(t)=e^{Qt}$. The p-value comes
from the $\chi^2_1$ upper tail (two-sided CONACC convention, validated vs the phyloP binary):

$$\mathrm{ceiling}(T)= -\log_{10}\operatorname{erfc}\!\Big(\sqrt{\mathrm{LRT}_{\max}/2}\Big)$$

## 4. What sets the ceiling: a race between depth and taxa

`LRT_max` is **not** a function of `T` alone. It has two limiting regimes; the ceiling is
governed by whichever constraint binds first.

### 4a. Tree-length limit (unsaturated branches) — the knob κ

Write the diagonal of the rate matrix as $-Q_{bb}=$ the total rate of evolving *away* from
base $b$ (its "exit rate"). PHAST normalizes $Q$ so the average exit rate is 1
($-\sum_b\pi_b Q_{bb}=1$; this is what puts branch lengths in subs/site). Individual bases
still differ; define

$$\boxed{\;\kappa=\max_b(-Q_{bb})\;}\qquad(\kappa\ge 1,\ \text{with }\kappa=1\text{ iff all bases evolve equally}).$$

When branches are short (many taxa), $P(t)\approx I+Qt$ and the pruning collapses to

$$\mathrm{LRT}_{\max}\;\longrightarrow\;2\,\kappa\,T \qquad(\text{many-taxon / unsaturated limit}).$$

So the model enters **only through κ**, which simply **stretches the depth axis by κ**: a
structured model detects like a $\kappa\times$ deeper tree. The most detectable site is one
conserved at the **fastest-evolving base** — conservation there is the most surprising.

**κ is defined for any reversible model, not just JC.** κ=1 is exactly Jukes–Cantor (all
rates equal). Any real asymmetry — unequal base frequencies, transition/transversion bias,
GC content — makes some bases turn over faster, so κ>1. For real nucleotide models κ≈1.3–1.5
(mammal 1.55, hamster 1.28, bird 1.50). JC is just the κ=1 special case of the same formula.

### 4b. Taxon-count limit (saturated branches) — a hard cap

`T` is the *total* expected substitutions, but a tree can only give you as many independent
observations as it has taxa. Push `T` up by making branches long, and each branch
**saturates**: its two ends become independent draws from equilibrium, $P(t)\to\pi$. Then a
perfectly conserved column is just $n$ independent "landed on base $b$" coincidences,
$P(\text{all}=b)\to\pi_b^{\,n}$, giving a **cap that depends on n, not T**:

$$\mathrm{LRT}_{\max}\;\longrightarrow\;2\,(n-1)\,(-\ln\pi_{\min})\qquad(\text{saturated limit}),$$

which for JC ($\pi_b=\tfrac14$) is $2(n-1)\ln 4$. Beyond saturation, **adding tree depth does
nothing** — you are taxon-limited. Numerically (JC), the ceiling caps at:

| taxa n | ceiling cap (−log10 p) | reached by ~T |
|---|---|---|
| 4  | 2.4  | 16 |
| 10 | 6.2  | 64 |
| 20 | 12.4 | 256 |
| 50 | 30.7 | large |

### The picture

$$\mathrm{LRT}_{\max}\ \approx\ \min\big(\,\underbrace{2\kappa T}_{\text{depth-limited}},\ \underbrace{2(n-1)(-\ln\pi_{\min})}_{\text{taxon-limited}}\,\big)$$

The crossover is at **T ≈ n** (i.e. average branch length ≈ 1 sub/site): with fewer, longer
branches you are taxon-limited; with many short branches you are depth-limited. The
"depends only on tree length" intuition is the **depth-limited regime** — true once you have
enough taxa (in practice the JC ceiling has essentially saturated the taxon effect by ~50
taxa). Sparse datasets (e.g. 15–20 species) are often *taxon-limited*: that is why the
15-taxon hamster's ceiling (7.0 at T=16) sits far below its depth-limited potential (~9.8),
despite a respectable T.

## 5. Closed forms

**Jukes–Cantor on a star tree of n tips at depth T** (purely theoretical, no fit):

$$\mathrm{LRT}_{\max}(n,T)=-2\ln\!\big(P_{\text{same}}^{\,n}+3P_{\text{diff}}^{\,n}\big),\quad
P_{\text{same}}=\tfrac14+\tfrac34e^{-4T/3n},\ P_{\text{diff}}=\tfrac14-\tfrac14e^{-4T/3n}$$

This single expression interpolates both regimes above (→ 2T as n→∞; → 2(n−1)ln4 as T→∞).

**Many-taxon limit** (used for the boundary band):

$$\mathrm{ceiling}(T)= -\log_{10}\operatorname{erfc}\!\big(\sqrt{\kappa T}\big)$$

## 6. Drawing the boundary

`ceiling(T)` is monotone increasing in T, so for each M the boundary depth is
$T^\*(M)=\mathrm{ceiling}^{-1}(\log_{10}(M/\alpha))$ (by interpolation); sites are detectable
for $T\ge T^\*(M)$. In `07`, the shaded band spans $\kappa\in[1,1.5]$ (flat → typical
nucleotide structure, many-taxon), and the dashed line is a 20-taxon JC tree, isolating the
finite-taxon (saturation) shift at κ=1.

`08_taxon_walls.py` makes the taxon effect explicit: it draws one boundary **per** taxon
count (n=8, 12, 20, ∞) at fixed κ=1, so the plane is a "minimum taxa needed" map. Each
boundary is a full detectable/undetectable line; with fewer taxa it shifts right (deeper
tree required) and flattens into a horizontal **taxon wall** at
$M_{\text{wall}}(n)=\alpha\cdot 10^{\mathrm{cap}(n)}$, $\mathrm{cap}(n)=-\log_{10}\operatorname{erfc}\!\big(\sqrt{(n-1)\ln 4}\big)$
(n=8 → M≈5e3, n=12 → 1.5e6, n=20 → 1.3e11). Above the wall no depth detects.

## 7. Caveats

- Rescaling a fitted model to an arbitrary depth stretches its branch lengths uniformly; far
  from the fitted depth this extrapolates the model shape.
- The $\chi^2_1$ p-value is phyloP's convention; the ceiling is the *best case* (perfectly
  conserved column at the most-informative base), so it is an upper bound on real per-site
  scores, which is exactly what a detectability limit needs.
- `06` draws the band from two real fitted models (illustrative, dataset-flavored); `07`
  draws it from κ and n as first-principles knobs (no fitted models).
