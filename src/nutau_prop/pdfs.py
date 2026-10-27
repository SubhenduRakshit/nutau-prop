class LHAPDFNucleon:
    """
    Thin LHAPDF wrapper. Isoscalar nucleon = 0.5*(p+n), using isospin
    u_n=d_p, d_n=u_p; sea counterparts are swapped similarly.
    """
    def __init__(self, setname="CT18NLO", member=0):
        try:
            import lhapdf
        except ImportError as exc:
            raise RuntimeError(
                "LHAPDF Python bindings not found. See README.md / setup_macos.sh."
            ) from exc
        self.lhapdf = lhapdf
        self.setname = setname
        self.member = int(member)
        self.pdf = lhapdf.mkPDF(setname, self.member)

    @property
    def q2_min(self):
        try:
            return float(self.pdf.q2Min)
        except Exception:
            return 1.0

    @property
    def q2_max(self):
        try:
            return float(self.pdf.q2Max)
        except Exception:
            return 1e12

    def xfx_proton(self, pid, x, q2):
        return float(self.pdf.xfxQ2(int(pid), float(x), float(q2)))

    @staticmethod
    def neutron_pid(pid):
        sign = 1 if pid >= 0 else -1
        a = abs(pid)
        if a == 2:
            return sign*1
        if a == 1:
            return sign*2
        return pid

    def xfx_isoscalar(self, pid, x, q2):
        p = self.xfx_proton(pid, x, q2)
        n = self.xfx_proton(self.neutron_pid(pid), x, q2)
        return 0.5*(p+n)
