"""Task 2 の本文反映（非中心化コサインへの減衰式、Figure S1）。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp, insert_before, append_section

D, M, S = 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md', 'manuscript/SUPPLEMENTARY.md'

SIM = """
This attenuation relation is derived for Pearson correlations, and cos θ here is not centred, so its
use is an approximation. Representational similarity analysis makes the same distinction between
metrics, normalising to unit norm for a cosine and to zero mean and unit variance for a correlation
[18,19]. We checked the approximation by simulating the whole procedure: two states with a known
angle between their true responses and a mean shift common to all genes, sampled at the group sizes
used here, with the ceiling estimated from the simulated samples exactly as above. The estimate is
not a strict upper bound, and how it behaves depends on the true angle. Restricting to the range of
estimated reliabilities our states actually have, r_half between 0.50 and 0.84, the observed value
never exceeded the estimate in any replicate at a true cos θ of 0.3 or 0.5 - the regime the
cross-species pairs occupy - while at 0.9 and at 1.0 it exceeded it in every replicate, and at 0.7 in
up to 46%. The size of the common mean shift, varied from 0 to 0.44 in units of the within-state
standard deviation against 0.26 in the median state of our panel, moves these figures by a few
percentage points and does not change the pattern. Two consequences follow. The estimate sits below
the value two identical states attain, so it understates the attainable similarity and the distance
we report to it is conservative. And it is not a limit for pairs that are already well aligned, which
is consistent with the same-species pairs sitting above their ceiling in Table 1. We therefore treat
it as a reference level rather than a hard limit, and rest the claims of Section 2.2 on the
animal-level interval of the paired difference rather than on a point comparison with it. Figure S1
gives the full grid.
"""

insert_before(M, "The uncorrected form applies no correction and is therefore the most conservative.",
              SIM)
rp(M, "The uncorrected form applies no correction and is therefore the most conservative.",
      "The uncorrected form applies no correction and is therefore the more conservative of the two.")

rp(D, "Ceiling estimates assume the two states' errors are independent, which fails for states "
      "sharing a cohort or control samples, and we observe exactly that failure, bounded at 0.094.",
      "Ceiling estimates assume the two states' errors are independent, which fails for states "
      "sharing a cohort or control samples, and we observe exactly that failure, bounded at 0.094. "
      "They rest on an attenuation relation derived for Pearson correlations and applied here to an "
      "uncentred cosine, which simulation shows is not a strict upper bound (Section 4.10, Figure "
      "S1): two states with the same true response exceed it in essentially every replicate, and "
      "pairs whose true alignment is high exceed it routinely. In the regime the cross-species pairs "
      "occupy, and at the reliabilities our states have, it was not exceeded in any simulated "
      "replicate, so the distance we report to it is if anything understated; but it is a reference "
      "level rather than a limit, and the claims here are made on the animal-level interval of the "
      "paired difference rather than on the point comparison.")

append_section(S, """**Figure S1.** *Does the attenuation relation bound an uncentred cosine?* Simulation of the
procedure of Section 4.10 at the group sizes used here. Each panel is one size of common mean shift,
in units of the within-state standard deviation; the median state of the panel has 0.26. Curves are
true angles between the two responses. The dashed line marks 1%. At the reliabilities our states
have (0.50-0.84) the estimate is never exceeded when the true alignment is 0.5 or below, and is
exceeded in every replicate when the two responses are identical, so it sits below the value
identical states attain (`results/revision1/cosine_attenuation_sim.tsv`).""")
print('Task 2 applied')
