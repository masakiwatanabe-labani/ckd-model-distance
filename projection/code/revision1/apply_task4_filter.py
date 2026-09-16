"""Task 4 の本文反映（猫側のフィルタ非対称性）。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

M = 'manuscript/METHODS.md'
rp(M, "Feline RNA. The deposited values are already on a log2 scale, so no further filter was "
      "applied beyond removing rows of zero variance.",
      "Feline RNA. The deposited values are already on a log2 scale, so no further filter was "
      "applied beyond removing rows of zero variance. Because this leaves the two species filtered "
      "differently, we checked what the asymmetry costs. Within Group A, on which every comparison "
      "is made, the two species occupy the same part of their own expression distributions: the "
      "median percentile rank of a Group A gene is 0.857 in the feline control samples against 0.810 "
      "and 0.828 in the two mouse control sets, and 1.4% of Group A genes fall below the 25th "
      "percentile in the feline data against 3.3-3.8% in the mouse data. Imposing a matched filter "
      "on the feline side as well - discarding genes below the 25th percentile of the feline "
      "expression distribution - raises the split-half reliability of all four feline states, by "
      "0.035, 0.008, 0.006 and 0.013, leaving their order and the split between the better- and "
      "worse-measured pair unchanged. Since the ceilings would rise with them, the filter choice is "
      "not what produces the unresolved pairs of Section 2.2.")
print('Task 4 applied')
