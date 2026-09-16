"""Task 1 と修正1・修正2 の本文反映（3群の gap、置換検定、B群への限定、置換の但し書き）。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp, insert_before

R, D, F = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/FRONTMATTER.md'

BLOCK = """
**The positive claim is limited to the pairs with a reliably measured feline state.** Of the three
groups, only group B - the 24 cat-mouse pairs built on feline cortical CKD3/4 or medullary CKD1/2 -
supports a statement that the observed value falls short of the ceiling. Group C, the 24 pairs built
on the two less reliably measured feline states, has 22 intervals containing zero and a median
difference of 0.215 against 0.119 for same-species pairs; we do not offer it as evidence in either
direction, and no claim below rests on it.

**What separates the classes is the size of the shortfall, not whether there is one.** Taken as a
yes-or-no question, falling short of the ceiling is not a cross-species phenomenon at all: the
interval excludes zero for 24 of the 39 same-species pairs (62%) and for 26 of the 48 cat-mouse pairs
(54%). The distributions differ in magnitude. Under the uncorrected ceiling the median difference is
0.119 (interquartile range 0.065-0.291) for same-species pairs, 0.403 (0.350-0.469) for cat-mouse
pairs with a reliably measured feline state and 0.215 (0.189-0.299) for the rest; under the
Spearman-Brown ceiling the three are 0.196 (0.117-0.360), 0.487 (0.437-0.561) and 0.338
(0.322-0.434). Permuting which four states are feline, exhaustively over all C(16,4) = 1,820
assignments and comparing the median difference of cross-species against same-species pairs, places
the observed gap of +0.219 against a null median of -0.012 (exact p = 0.0011); under the
Spearman-Brown ceiling the observed gap is +0.239 (p = 0.0005). These p values carry the same
qualification as the permutation test of the area under the curve above. Species is perfectly
confounded with dataset here, so the null contains assignments in which the four feline states are
split between species, which cannot occur; the stricter permutation over which dataset is designated
feline has a floor of p = 1/3. They therefore measure how unusual the observed split is among
relabellings of the states, and are not evidence of a species effect. The horizontal lines in Figure
2D mark these group medians. What the ceiling analysis establishes is a difference of degree between
cross-species and same-species shortfalls, not a qualitative property of cross-species comparison.
"""

insert_before(R, "Three limits on how far this goes.", BLOCK)

rp(D, "What this identifies is the range over which this dataset and this estimation procedure can",
      "We take the first group as the claim and treat the second as uninformative rather than as weak "
      "support for it: a median difference of 0.215 there, against 0.119 among same-species pairs, is "
      "not something we read as a shortfall. Falling short of the ceiling is not in itself a "
      "cross-species phenomenon: it happens about as often within species as across them - the "
      "interval excludes zero for 62% of same-species pairs and 54% of cat-mouse pairs - and what "
      "separates the two classes is how large the shortfall is. The median difference is 0.119 within "
      "species against 0.403 and 0.215 in the two cross-species groups, and permuting which four "
      "states are feline places the observed difference in medians, +0.219, against a null median of "
      "-0.012 (exact p = 0.0011) - a p value that, like the one attached to the area under the curve, "
      "is computed over relabellings of states drawn from datasets that are themselves confounded "
      "with species, and so describes the observed split rather than testing a species effect. The "
      "cross-species shortfall is continuous with the same-species one, and the evidence supports a "
      "difference of degree and nothing stronger. What this identifies is the range over which this "
      "dataset and this estimation procedure can")

rp(F, "Second, how far a model sits below the ceiling set by measurement reliability can be judged "
      "only where the reference state is well estimated: resampling animals rather than genes, the "
      "gap survives in the 24 cat-mouse pairs built on the two reliable feline states and cannot be "
      "resolved in 22 of the others.",
      "Second, falling below the ceiling set by measurement reliability is not peculiar to "
      "cross-species pairs, occurring about as often within species; they differ in the size of the "
      "shortfall, and the claim holds only for the 24 cat-mouse pairs with a reliably measured feline "
      "state (median 0.403 against 0.119 within species).")
print('Task 1 applied')
