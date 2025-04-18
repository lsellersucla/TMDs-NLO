import lhapdf
import numpy as np
import sys

try:
    from Collins_fit.TMDs.Numerical.FBT import FBT
    from Collins_fit.TMDs.tmds import TMDs
    from Collins_fit.kits.share import share
except:
    sys.path.append("../")
    sys.path.append("../..")
    from TMDs.Numerical.FBT import FBT
    from TMDs.tmds import TMDs
    from kits.share import share

import warnings
warnings.filterwarnings("ignore")

class UDIS:

    def __init__(self, loops, logs):
        self.ffpip = share['LHAPDF']['FFs' ]
        self.pdfp  = share['LHAPDF']['PDFs']
        self.tmds = TMDs(loops, logs, None)

        self.scheme = share['scheme']
        self.fbt0 = FBT(0)

    #-- Import unpolarized TMD PDF from tmds.py
    def TMDPDF(self,x,Q,b):
        return self.tmds.TMDPDF(x,Q,b)

    #-- Import unpolarized TMD FF  from tmds.py
    def TMDFF( self,z,Q,b):
        return self.tmds.TMDFF( z,Q,b)

    #-- The unpolarized cross section in b-space
    def dsigmadxdydzd2PhT_b_SIDIS(self,b,x,z,Q):
        Fgg,FDD,FUU,FSS,FSB,FUB,FDB = self.TMDPDF(x,Q,b)
        Dgg,DDD,DUU,DSS,DSB,DUB,DDB = self.TMDFF(z,Q,b)
        als = self.tmds.alphas(Q)
        CF = 4./3.
        if self.scheme == 'MSbar':
            HQ = 1+als*CF/2./np.pi*(-8.+np.pi**2./6.)
        elif self.scheme == 'JCC':
            HQ = 1+als*CF/2./np.pi*(-8.)
        elif self.scheme == 'CSS':
            HQ = 1
        eu2 = 4./9.
        ed2 = 1./9.
        return HQ*(eu2*(FUU*DUU+FUB*DUB)+ed2*(FDD*DDD+FDB*DDB+FSS*DSS+FSB*DSB))

    #-- The unpolarized cross section in momentum space (ds/dx dy dz d2PhT)
    def dsigmadxdydzd2PhT_p_SIDIS(self,PhT,x,y,z,Q):
        sigma0 = 2.*np.pi*self.tmds.alphaEM(Q)**2./Q**4.*(1.+(1.-y)**2.)/y # sigma0 for dx dy dz d2PhT
        FUUb = lambda b: b*self.dsigmadxdydzd2PhT_b_SIDIS(b,x,z,Q)
        return sigma0*self.fbt0.fbt(FUUb, PhT/z, 20,Q)

    #-- The unpolarized cross section in momentum space (ds/dx dQ^2 dz d2PhT)
    def dsigmadxdQ2dzd2PhT_p_SIDIS(self,PhT,x,y,z,Q):
        sigma0 = 2.*np.pi*self.tmds.alphaEM(Q)**2./Q**4.*(1.+(1.-y)**2.)  # sigma0 for dx dQ^2 dz d2PhT
        FUUb = lambda b: b*self.dsigmadxdydzd2PhT_b_SIDIS(b,x,z,Q)
        return sigma0*self.fbt0.fbt(FUUb, PhT/z, 20,Q)

