"""rev4 の OMML 数式は平文化して壊れるので、LaTeX 表記に置き換える。

置き換え先は以前の Markdown で検証済みの式と同一。数値や定義は変えていない。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

M = 'manuscript/METHODS.md'
R = 'manuscript/RESULTS.md'

FORMULAS = [
    (M, "Δg,s = 1si∈s​yg,i - 1ctrlsi∈ctrls​yg,i",
     r"$$\Delta_{g,s} \;=\; \frac{1}{|s|}\sum_{i \in s} y_{g,i} \;-\; "
     r"\frac{1}{|\mathrm{ctrl}(s)|}\sum_{i \in \mathrm{ctrl}(s)} y_{g,i}$$"),
    (M, "zg,i=yg,i-yg,ctrl-meani⋅sdi⋅,  Δg,sz=zg,s-zg,ctrl",
     r"$$z_{g,i} = \frac{(y_{g,i} - \bar{y}_{g,\mathrm{ctrl}}) - \mathrm{mean}_i(\cdot)}"
     r"{\mathrm{sd}_i(\cdot)}, \qquad \Delta^{z}_{g,s} = \bar{z}_{g,s} - \bar{z}_{g,\mathrm{ctrl}}$$"),
    (M, "α=νcosθ,  R⊥=1-cos2θ=sinθ",
     r"$$\alpha = \nu \cos\theta, \qquad R_\perp = \sqrt{1 - \cos^2\theta} = \sin\theta$$"),
    (M, "Δg,s1A=yg,s1-yg,A,  Δg,s2B=yg,s2-yg,B,  A∩B=∅",
     r"$$\Delta^{(A)}_{g,s_1} = \bar{y}_{g,s_1} - \bar{y}_{g,A}, \qquad "
     r"\Delta^{(B)}_{g,s_2} = \bar{y}_{g,s_2} - \bar{y}_{g,B}, \qquad A \cap B = \emptyset$$"),
    (M, "rXX=2 rhalf1+rhalf,  rXX:=0  if rhalf≤0",
     r"$$r_{XX} = \frac{2\,r_{\text{half}}}{1 + r_{\text{half}}}, \qquad "
     r"r_{XX} := 0 \ \text{ if } r_{\text{half}} \le 0$$"),
    (M, "ceilingXYraw=rhalf,X rhalf,Y  ceilingXYSB=rXXSB rYYSB",
     r"$$\mathrm{ceiling}^{\text{raw}}_{XY} = \sqrt{r_{\text{half},X}\, r_{\text{half},Y}} "
     r"\qquad \mathrm{ceiling}^{\text{SB}}_{XY} = \sqrt{r^{\text{SB}}_{XX}\, r^{\text{SB}}_{YY}}$$"),
    (M, "zg,s=Δg,s-meangΔ⋅,ssdgΔ⋅,s,  mP,s=1Pg∈P​zg,s",
     r"$$z_{g,s} = \frac{\Delta_{g,s} - \mathrm{mean}_g(\Delta_{\cdot,s})}"
     r"{\mathrm{sd}_g(\Delta_{\cdot,s})}, \qquad "
     r"m_{P,s} = \frac{1}{|P|}\sum_{g \in P} z_{g,s}$$"),
]


def main():
    for path, old, new in FORMULAS:
        rp(path, old, new)
    print('formulas restored:', len(FORMULAS))


if __name__ == "__main__":
    main()
