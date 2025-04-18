import lhapdf
import numpy as np
from scipy.special import gamma,digamma,beta
from scipy.integrate import quad,fixed_quad
import sys
sys.path.append("../")
sys.path.append("../..")
from TMDs.Evolve.anom import ANOM
from TMDs.Evolve.alpha import ALPHAS
from TMDs.Numerical.FBT import FBT
from kits.share import share
from TMDs.tmds import TMDs

import warnings
warnings.filterwarnings("ignore")

class UDIA:

    def __init__(self,loops,logs,collFF = None):
        if (collFF == None) and (collPDF == None):
            self.ffpip = share['LHAPDF']['FFs' ]
        self.loops = loops
        self.logs  = logs
        self.ALFS = ALPHAS(loops)
        self.anom = ANOM(loops,logs)
        self.tmds = TMDs(loops,logs)
        self.fbt0 = FBT(0)

    #-- Strong coupling constant
    def alphas(self,Q):
        return self.ALFS.alfs(Q)

    #-- EM coupling constant
    def alphaEM(self,Q):
        return self.ALFS.alfE(Q)

    #-- Import unpolarized TMD FF  from tmds.py
    def TMDFF( self,z,Q,b):
        return self.TMDs.TMDFF( z,Q,b)
    
    #-- The DIA unpolarized cross section in b-space
    def dsigmadz1dz2d2PhTdcosH_b_DIA(self,b,z1,z2,Q):
        Dgg1,DDD1,DUU1,DSS1,DSB1,DUB1,DDB1 = self.TMDFF(z1,Q,b)
        Dgg2,DDD2,DUU2,DSS2,DSB2,DUB2,DDB2 = self.TMDFF(z2,Q,b)
        als = self.alphas(Q)
        CF = 4./3.
        HQ = 1+als*CF/2./np.pi*(-8.+7.*np.pi**2./2.) # H is not absorbed into the TMDs
        eu2 = 4./9.
        ed2 = 1./9.
        return HQ*(eu2*(DUU1*DUB2+DUB1*DUU2)+ed2*(DDD1*DDB2+DDB1*DDD2+DSS1*DSB2+DSB1*DSS2))
    
    #-- The unpolarized DIA cross section in momentum space (ds/dz1 dz2 d2PhT dcosQ d2PhT)
    def dsigmadz1dz2d2PhTdcosH_p_DIA(self,PhT,z1,z2,theta,Q):
        Nc = 3.
        sigma0 = Nc*np.pi*self.alphaEM(Q)**2./(2.*Q**2.)*(1+np.cos(theta)**2.)
        ZUUb = lambda b: b*self.dsigmadz1dz2d2PhTdcosH_b_DIA(b,z1,z2,Q)
        return sigma0*self.fbt0.fbt(FUUb, PhT/z1, 20,Q)

