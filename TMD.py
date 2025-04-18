import lhapdf
import numpy as np
from scipy.integrate import quad

try:
    from Collins_fit.TMDs.Evolve.alpha import ALPHAS
except:
    from Evolve.alpha import ALPHAS

class TMDPDF:
    def __init__(self,order, x_cut, pdf_set, member,error,scheme):
        self.order  = order
        self.pdf    = lhapdf.mkPDF(pdf_set, member)
        self.alphaS = ALPHAS(order)
        self.x_cut  = x_cut
        self.error  = error

        ## constants
        self.CF = 4.0 / 3.0
        self.CA = 3.0
        self.TF = 1.0 / 2.0

        ## choose scheme
        self.scheme = scheme

    def alf4pi(self,Q):
        return self.alphaS.alfs(Q)/4./np.pi

    def PDF(self,x,Q):
        if Q<1.:
            Q = 1.
        if x>self.x_cut:
            x = self.x_cut
        g = self.pdf.xfxQ( 0,x,Q)/x
        d = self.pdf.xfxQ( 1,x,Q)/x
        u = self.pdf.xfxQ( 2,x,Q)/x
        s = self.pdf.xfxQ( 3,x,Q)/x
        c = self.pdf.xfxQ( 4,x,Q)/x
        b = self.pdf.xfxQ( 5,x,Q)/x
        bb= self.pdf.xfxQ(-5,x,Q)/x
        cb= self.pdf.xfxQ(-4,x,Q)/x
        sb= self.pdf.xfxQ(-3,x,Q)/x
        ub= self.pdf.xfxQ(-2,x,Q)/x
        db= self.pdf.xfxQ(-1,x,Q)/x
        return g,d,u,s,c,b,bb,cb,sb,ub,db

    def PDFg(self,x,Q):
        if Q<1.:
            Q = 1.
        if x>self.x_cut:
            x = self.x_cut
        g = self.pdf.xfxQ( 0,x,Q)/x
        return g

    def PDFqsum(self,x,Q):
        if Q<1.:
            Q =	1.
        if x>self.x_cut:
            x =	self.x_cut
        d = self.pdf.xfxQ( 1,x,Q)/x
        u = self.pdf.xfxQ( 2,x,Q)/x
        s = self.pdf.xfxQ( 3,x,Q)/x
        c = self.pdf.xfxQ( 4,x,Q)/x
        b = self.pdf.xfxQ( 5,x,Q)/x
        bb= self.pdf.xfxQ(-5,x,Q)/x
        cb= self.pdf.xfxQ(-4,x,Q)/x
        sb= self.pdf.xfxQ(-3,x,Q)/x
        ub= self.pdf.xfxQ(-2,x,Q)/x
        db= self.pdf.xfxQ(-1,x,Q)/x
        return d+u+s+c+b+bb+cb+sb+ub+db

    def pqqNLO(self,x):
        return (1.+x**2.)/(1.-x)

    def pqgNLO(self,x):
        xb = 1.-x
        return 1.-2.*x*xb

    def pgqNLO(self,x):
        xb = 1.-x
        return (1.+xb**2.)/x

    def pggNLO(self,x):
        xb = 1.-x
        return (1.-x*xb)**2./(x*(1.-x))

    def QQintNLO(self,x,Lmu,Lz):
        CF = 4./3.
        xb = 1.-x
        return CF*(-2.*Lmu*self.pqqNLO(x)+2*xb)

    def QQdivNLO(self,x,Lmu,Lz):
        ## subtraction term from the plus distribution
        CF = 4./3.
        return 4.*CF*Lmu/(-1.+x)

    def QQdelta(self,Lmu,Lz):
        ## delta(1-x) terms
        CF = 4./3.
        if self.scheme == 'MSbar':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.-np.pi**2./6.)
        elif self.scheme == 'JCC':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.)
        elif self.scheme == 'CSS':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.-8.0)
        return kernel

    def QGintNLO(self,x,Lmu,Lz):
        Tr = 1./2.
        xb = 1.-x
        return Tr*(-2.*Lmu*self.pgqNLO(x)+4.*x*xb)

    def GQintNLO(self,x,Lmu,Lz):
        CF = 4./3.
#        return CF*(-2.*Lmu*self.pqgNLO(x)+2.*x)
        return CF*(-2.*Lmu*self.pgqNLO(x)+2.*x)

    def GGintNLO(self,x,Lmu,Lz):
        CA = 3.
        return CA*(-4.*Lmu*self.pggNLO(x))

    def GGdivNLO(self,x,Lmu,Lz):
        ## subtraction term from the plus distribution
        CA = 3.
        return (4*CA*Lmu)/(-1 + x)

    def GGdelta(self,Lmu,Lz):
        ## delta(1-x) terms
        CA = 3.
        if self.scheme == 'MSbar':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu - np.pi**2./6.)
        elif self.scheme == 'JCC':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu)
        elif self.scheme == 'CSS':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu-8.0)
        return kernel

    def CxPDF(self,x,Q,Lmu,Lz):
        glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo = self.PDF(x,Q)
        if self.order == 0:
            vals = [glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo]
        elif self.order == 1:
            as4pi      = self.alf4pi(Q)
            ggintnlo   = lambda xp: as4pi / xp * (self.GGintNLO(xp,Lmu,Lz)*self.PDFg(x/xp,Q) - \
                                                  self.GGdivNLO(xp,Lmu,Lz)*glo + \
                                                  self.GQintNLO(xp,Lmu,Lz)*self.PDFqsum(x/xp,Q))
            ggnloint   = quad(ggintnlo, x, self.x_cut, epsrel=self.error)[0]
            ggnlodel   = as4pi*self.GGdelta(Lmu,Lz)
            g = glo+ggnloint+glo*ggnlodel

            valslo = np.array([glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo])
            vals = [g]
            for i in range(1,11):
                lo = valslo[i]
                intnlo = lambda xp: as4pi / xp * (self.QQintNLO(xp,Lmu,Lz)*self.PDF(x/xp,Q)[i] - \
                                                  self.QQdivNLO(xp,Lmu,Lz)*lo + \
                                                  self.QGintNLO(xp,Lmu,Lz)*self.PDFg(x/xp,Q))
                nloint = quad(intnlo, x, self.x_cut, epsrel=self.error)[0]
                nlodel = as4pi*self.QQdelta(Lmu,Lz)
                vals.append(lo+nloint+lo*nlodel)

        return np.array(vals)


class TMDFF:
    def __init__(self,order,z_cut, ff_set, member,error,scheme):
        self.order    = order
        self.pdf = lhapdf.mkPDF(ff_set, member)
        self.alphaS = ALPHAS(order)
        self.z_cut = z_cut
        self.error = error

        ## choose scheme
        self.scheme = scheme

    def alf4pi(self,Q):
        return self.alphaS.alfs(Q)/4./np.pi

    def FF(self,z,Q):
       	if Q<1.:
       	    Q =	1.
        if z>self.z_cut:
            z = self.z_cut
        g = self.pdf.xfxQ( 0,z,Q)/z
        d = self.pdf.xfxQ( 1,z,Q)/z
        u = self.pdf.xfxQ( 2,z,Q)/z
        s = self.pdf.xfxQ( 3,z,Q)/z
        c = self.pdf.xfxQ( 4,z,Q)/z
        b = self.pdf.xfxQ( 5,z,Q)/z
        bb= self.pdf.xfxQ(-5,z,Q)/z
        cb= self.pdf.xfxQ(-4,z,Q)/z
        sb= self.pdf.xfxQ(-3,z,Q)/z
        ub= self.pdf.xfxQ(-2,z,Q)/z
        db= self.pdf.xfxQ(-1,z,Q)/z
        return g,d,u,s,c,b,bb,cb,sb,ub,db

    def FFg(self,z,Q):
       	if Q<1.:
       	    Q =	1.
        if z>self.z_cut:
            z = self.z_cut
        g = self.pdf.xfxQ( 0,z,Q)/z
        return g

    def FFqsum(self,z,Q):
       	if Q<1.:
       	    Q =	1.
        if z>self.z_cut:
            z = self.z_cut
        d = self.pdf.xfxQ( 1,z,Q)/z
        u = self.pdf.xfxQ( 2,z,Q)/z
        s = self.pdf.xfxQ( 3,z,Q)/z
        c = self.pdf.xfxQ( 4,z,Q)/z
        b = self.pdf.xfxQ( 5,z,Q)/z
        bb= self.pdf.xfxQ(-5,z,Q)/z
        cb= self.pdf.xfxQ(-4,z,Q)/z
        sb= self.pdf.xfxQ(-3,z,Q)/z
        ub= self.pdf.xfxQ(-2,z,Q)/z
        db= self.pdf.xfxQ(-1,z,Q)/z
        return d+u+s+c+b+bb+cb+sb+ub+db

    def pqqNLO(self,z):
        return (1.+z**2.)/(1.-z)

    def pqgNLO(self,z):
        zb = 1.-z
        return 1.-2.*z*zb

    def pgqNLO(self,z):
        zb = 1.-z
        return (1.+zb**2.)/z

    def pggNLO(self,z):
        zb = 1.-z
        return (1.-z*zb)**2./(z*(1.-z))

    def QQintNLO(self,z,Lmu,Lz):
        CF = 4./3.
        zb = 1.-z
        return CF*(-2.*Lmu*self.pqqNLO(z)+2*zb)

    def QQdivNLO(self,z,Lmu,Lz):
        ## subtraction term from the plus distribution
        CF = 4./3.
        return 4.*CF*Lmu/(-1.+z)

    def QQdelta(self,Lmu,Lz):
        ## delta(1-z) terms
        CF = 4./3.
        if self.scheme == 'MSbar':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.-np.pi**2./6.)
        elif self.scheme == 'JCC':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.)
        elif self.scheme == 'CSS':
            kernel = CF*(2.*Lz*Lmu-Lmu**2.-8.0)
        return kernel

    def GQintNLO(self,z,Lmu,Lz):
        CF = 4./3.
#        return CF*(-2.*Lmu*self.pqgNLO(z)+2.*z)
        return CF*(-2.*Lmu*self.pgqNLO(z)+2.*z)

    def QGintNLO(self,z,Lmu,Lz):
        Tr = 1./2.
        zb = 1.-z
        return Tr*(-2.*Lmu*self.pgqNLO(z)+4.*z*zb)

    def GGintNLO(self,z,Lmu,Lz):
        CA = 3.
        return CA*(-4.*Lmu*self.pggNLO(z))

    def GGdivNLO(self,z,Lmu,Lz):
        ## subtraction term from the plus distribution
        CA = 3.
        return (4*CA*Lmu)/(-1 + z)

    def GGdelta(self,Lmu,Lz):
        ## delta(1-z) terms
        CA = 3.
        if self.scheme == 'MSbar':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu - np.pi**2./6.)
        elif self.scheme == 'JCC':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu)
        elif self.scheme == 'CSS':
            kernel = CA*(2.*Lz*Lmu -Lmu*Lmu-8.0)
        return kernel

    def CxFF(self,z,Q,Lmu,Lz):
        Lmuz  = Lmu-np.log(z**2.)
        Lmuzp = lambda zp: Lmu-np.log(zp**2.)
        glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo = self.FF(z,Q)
        if self.order == 0:
            vals = [glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo]
        elif self.order == 1:
            as4pi      = self.alf4pi(Q)
            ggintnlo   = lambda zp: as4pi / zp * (self.GGintNLO(zp,Lmu,Lz)*self.FFg(z/zp,Q) - \
                                                  self.GGdivNLO(zp,Lmu,Lz)*glo + \
                                                  self.QGintNLO(zp,Lmu,Lz)*self.FFqsum(z/zp,Q))
            ggnloint   = quad(ggintnlo, z, self.z_cut, epsrel=self.error)[0]
            ggnlodel   = as4pi*self.GGdelta(Lmu,Lz)
            g = glo+ggnloint+glo*ggnlodel

            valslo = np.array([glo,dlo,ulo,slo,clo,blo,bblo,cblo,sblo,ublo,dblo])
            vals = [g]
            for i in range(1,11):
                lo = valslo[i]
                # Note that the integration of the Log(zp) term in QQintNLO is well defined as z-> 1. Thus we only regulate the Lmu term with QQdivNLO.
                intnlo = lambda zp: as4pi / zp * (self.QQintNLO(zp,Lmuzp(zp),Lz)*self.FF(z/zp,Q)[i] - \
                                                  self.QQdivNLO(zp,Lmu,Lz)*lo + \
                                                  self.GQintNLO(zp,Lmuzp(zp),Lz)*self.FFg(z/zp,Q))
                nloint = quad(intnlo, z, self.z_cut, epsrel=self.error)[0]
                # Note that the delta function has zp = 1. Thus we just use Lmu
                nlodel = as4pi*self.QQdelta(Lmu,Lz)
                vals.append(lo+nloint+lo*nlodel)
        return np.array(vals)/z/z
