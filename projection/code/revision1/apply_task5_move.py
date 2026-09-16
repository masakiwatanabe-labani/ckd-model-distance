"""Task 5 の本文反映（Discussion 3.6 を Results 2.7 へ）。"""
import sys, re, pathlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp, wrap

D = pathlib.Path('manuscript/DISCUSSION.md')
R = pathlib.Path('manuscript/RESULTS.md')

t = D.read_text()
secs = re.split(r'(?m)^(?=## 3\.\d+\.)', t)
head, secs = secs[0], secs[1:]
i = next(k for k, s in enumerate(secs) if s.startswith('## 3.6.'))
body = secs[i].split('\n', 1)[1].strip()

secs[i] = """## 3.6. Extrapolation to human disease

The exploratory comparison with human tubulointerstitial profiles is reported in Section 2.7. Its
lesson is procedural rather than biological. A comparison in which two species reach the human
reference through orthologue maps of unequal quality is not a comparison of the species, and the
control that detects this - matching genes on mapping quality - reverses the sign of the effect
rather than attenuating it. We therefore treat the result as uninformative about species and retain
it only as an illustration of the kind of asymmetry that more human cohorts would not repair.

"""
D.write_text(head + ''.join(secs))

rt = R.read_text()
anchor = '\n---\n\n### Figures and tables'
assert rt.count(anchor) == 1
new = '## 2.7. Extrapolation to human disease does not survive its own control\n\n' + body + '\n'
R.write_text(rt.replace(anchor, '\n' + new + anchor))

rp('manuscript/BACKMATTER.md', "the exploratory analysis of Section 3.6", "the exploratory analysis of Section 2.7")
rp('manuscript/METHODS.md', "Because the result does not survive that control (Section 3.6)",
   "Because the result does not survive that control (Section 2.7)")
print('Task 5 applied')
