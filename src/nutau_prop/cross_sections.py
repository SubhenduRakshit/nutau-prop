import numpy as np
from numpy.polynomial.legendre import leggauss
from dataclasses import dataclass
from .constants import GF, M_N, M_W, M_Z, GEV2_TO_CM2, SIN2W_DEFAULT
from .bsm import BSMModel, delta_gv

UP_TYPES = (2, 4)          # u,c; top omitted in v0.1
DOWN_TYPES = (1, 3, 5)    # d,s,b
QCHARGE = {1:-1/3, 2:2/3, 3:-1/3, 4:2/3, 5:-1/3}

def nc_chiral_couplings(pid, sin2w):
    q = QCHARGE[abs(pid)]
    if abs(pid) in UP_TYPES:
        t3 = 0.5
    else:
        t3 = -0.5
    gL = t3 - q*sin2w
    gR = -q*sin2w
    return gL, gR, q

@dataclass
class DISConfig:
    q2_min_gev2: float = 2.0
    w2_min_gev2: float = 4.0
    nx: int = 72
    sin2w: float = SIN2W_DEFAULT

class DISCalculator:
    """
    Leading-order neutrino DIS with LHAPDF PDFs.

    Normalization follows the standard chiral-coupling form:
      d2sigma_NC/dxdy = 2 GF^2 M E / pi * P_Z *
        sum_q [ gL^2 (xq + xqbar(1-y)^2)
              +gR^2 (xq(1-y)^2 + xqbar) ].

    Charge radius shifts the vector coupling. Since gV=gL+gR and
    gA=gL-gR, delta gV is implemented as delta gL=delta gR=delta gV/2.
    """
    def __init__(self, pdf, config=None):
        self.pdf = pdf
        self.cfg = config or DISConfig()
        self._nodes, self._weights = leggauss(self.cfg.nx)

    def _x_bounds(self, E, y):
        denom = 2*M_N*E*y
        if denom <= 0:
            return None
        xmin = self.cfg.q2_min_gev2/denom
        xmax_w = 1.0 - max(0.0, self.cfg.w2_min_gev2-M_N*M_N)/denom
        xmax = min(1.0-1e-10, xmax_w)
        xmin = max(xmin, 1e-9)
        if xmin >= xmax:
            return None
        return xmin, xmax

    def _integrate_logx(self, func, xmin, xmax):
        a, b = np.log(xmin), np.log(xmax)
        t = 0.5*(self._nodes+1)*(b-a)+a
        x = np.exp(t)
        jac = x * 0.5*(b-a)
        vals = np.array([func(xx) for xx in x])
        return float(np.sum(self._weights * vals * jac))

    def d2_nc(self, E, x, y, antinu=False, model=None):
        model = model or BSMModel()
        q2 = 2*M_N*E*x*y
        if q2 < self.cfg.q2_min_gev2:
            return 0.0
        w2 = M_N*M_N + 2*M_N*E*y*(1-x)
        if w2 < self.cfg.w2_min_gev2:
            return 0.0
        prop = M_Z**4/(q2+M_Z**2)**2
        s = 0.0
        omy2 = (1-y)**2
        for pid in (1,2,3,4,5):
            xq = self.pdf.xfx_isoscalar(pid, x, q2)
            xqb = self.pdf.xfx_isoscalar(-pid, x, q2)
            gL,gR,qcharge = nc_chiral_couplings(pid, self.cfg.sin2w)
            dgv = delta_gv(qcharge, q2, model, self.cfg.sin2w, antinu=antinu)
            gL += 0.5*dgv
            gR += 0.5*dgv
            if not antinu:
                s += gL*gL*(xq+xqb*omy2) + gR*gR*(xq*omy2+xqb)
            else:
                s += gL*gL*(xq*omy2+xqb) + gR*gR*(xq+xqb*omy2)
        pref = 2*GF*GF*M_N*E/np.pi
        return pref*prop*s*GEV2_TO_CM2

    def d2_cc(self, E, x, y, antinu=False):
        q2 = 2*M_N*E*x*y
        if q2 < self.cfg.q2_min_gev2:
            return 0.0
        w2 = M_N*M_N + 2*M_N*E*y*(1-x)
        if w2 < self.cfg.w2_min_gev2:
            return 0.0
        prop = M_W**4/(q2+M_W**2)**2
        omy2 = (1-y)**2
        down = sum(self.pdf.xfx_isoscalar(pid,x,q2) for pid in DOWN_TYPES)
        up = sum(self.pdf.xfx_isoscalar(pid,x,q2) for pid in UP_TYPES)
        adown = sum(self.pdf.xfx_isoscalar(-pid,x,q2) for pid in DOWN_TYPES)
        aup = sum(self.pdf.xfx_isoscalar(-pid,x,q2) for pid in UP_TYPES)
        if not antinu:
            s = down + aup*omy2
        else:
            s = adown + up*omy2
        pref = 2*GF*GF*M_N*E/np.pi
        return pref*prop*s*GEV2_TO_CM2

    def dsdy(self, E, y, current="NC", antinu=False, model=None):
        b = self._x_bounds(E,y)
        if b is None:
            return 0.0
        xmin,xmax = b
        if current.upper()=="NC":
            f = lambda x: self.d2_nc(E,x,y,antinu=antinu,model=model)
        elif current.upper()=="CC":
            f = lambda x: self.d2_cc(E,x,y,antinu=antinu)
        else:
            raise ValueError("current must be NC or CC")
        return self._integrate_logx(f,xmin,xmax)

    def y_grid(self, y_min=1e-5, y_split=0.08, nlog=55, nlin=55):
        e1 = np.geomspace(y_min,y_split,nlog+1)
        e2 = np.linspace(y_split,1.0-1e-6,nlin+1)
        edges = np.unique(np.r_[e1,e2[1:]])
        centers = 0.5*(edges[:-1]+edges[1:])
        widths = np.diff(edges)
        return edges,centers,widths

    def spectrum(self,E,current="NC",antinu=False,model=None,**ykw):
        edges,y,dy = self.y_grid(**ykw)
        ds = np.array([self.dsdy(E,yy,current=current,antinu=antinu,model=model) for yy in y])
        sigma = float(np.sum(ds*dy))
        p = ds*dy/sigma if sigma>0 else np.zeros_like(ds)
        return dict(edges=edges,y=y,dy=dy,dsdy=ds,sigma=sigma,p=p)
