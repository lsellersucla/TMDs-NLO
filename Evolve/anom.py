import numpy as np
import sys
sys.path.append("../../")
from scipy.integrate import quad
from TMDs.Evolve.alpha import ALPHAS

#try:
#    from Collins_fit.TMDs.Evolve.alpha import ALPHAS
#except:
#    from TMDs.Evolve.alpha import ALPHAS



class ANOM:
    def __init__(self,order,logorder):
        self.order     = order
        self.logorder  = logorder
        self.alphaS    = ALPHAS(order)
        self.c0        = 1.12291896713377
        self.c02       = 1.260947006748773
        self.CA        = 3.
        self.CF        = 4./3.
        self.mc        = 1.5
        self.mb        = 4.7

        self.CuspG     = {}
        self.CuspG[3]  = [0.,12.00000000000000, 109.5647471869277, 1391.0059561317298, 17659.9]
        self.CuspG[4]  = [0.,12.00000000000000, 96.23141385359438, 966.3897931924722 , 9704.83]
        self.CuspG[5]  = [0.,12.00000000000000, 82.89808052026105, 538.2180746976597 , 3494.4 ]

        self.NCuspG    = {}
        self.NCuspG[3] = [0.,-18.000000000000000,-227.8995118764504 ,-3957.378545284555 ]
        self.NCuspG[4] = [0.,-16.666666666666668,-200.7014703660655 ,-3206.1850399832547]
        self.NCuspG[5] = [0.,-15.333333333333334,-173.50342885568062,-2485.5387880396356]

        self.DtermG    = {}
        self.DtermG[11]= [0.,0.,0., 6.000000000000000, 6.000000000000000, 6.000000000000000]
        self.DtermG[10]= [0.,0.,0., 0.               , 0.               , 0.               ]
        self.DtermG[22]= [0.,0.,0., 27.00000000000000, 25.00000000000000, 23.00000000000000]
        self.DtermG[21]= [0.,0.,0., 54.78237359346385, 48.11570692679718, 41.44904026013051]
        self.DtermG[20]= [0.,0.,0.,-35.45916979810888,-41.68139202033110,-47.90361424255332]
        self.DtermG[33]= [0.,0.,0., 162.0000000000000, 138.8888888888889, 117.5555555555556]
        self.DtermG[32]= [0.,0.,0., 685.0413623411746, 554.9642243899765, 433.7759753276673]
        self.DtermG[31]= [0.,0.,0., 57.23792169990449,-211.4949704092826,-465.4130477036549]
        self.DtermG[30]= [0.,0.,0.,-586.2020781721192,-688.1321632033369,-769.6024478518811]

        self.CuspQ     = {}
        self.CuspQ[3]  = [0.,5.333333333333333, 48.69544319419009, 618.22486939187980, 7848.82]
        self.CuspQ[4]  = [0.,5.333333333333333, 42.76951726826417, 429.50657475220990, 4313.26]
        self.CuspQ[5]  = [0.,5.333333333333333, 36.84359134233824, 239.20803319895987, 1553.06]

        self.NCuspQ    = {}
        self.NCuspQ[3] = [0.,-8.00000000000000,-29.243530284415503,-738.2562930508085]
        self.NCuspQ[4] = [0.,-8.00000000000000,-14.050795508138547,-491.96573445169145]
        self.NCuspQ[5] = [0.,-8.00000000000000, 1.1419392681384102,-249.38756710544408]

        self.DtermQ    = {}
        self.DtermQ[11]= [0.,0.,0., 2.666666666666667, 2.666666666666667, 2.666666666666667]
        self.DtermQ[10]= [0.,0.,0., 0.               , 0.               , 0.               ]
        self.DtermQ[22]= [0.,0.,0., 12.00000000000000, 11.11111111111111, 10.22222222222222]
        self.DtermQ[21]= [0.,0.,0., 24.34772159709504, 21.38475863413208, 18.42179567116912]
        self.DtermQ[20]= [0.,0.,0.,-15.75963102138172,-18.52506312014716,-21.29049521891259]
        self.DtermQ[33]= [0.,0.,0., 72.00000000000000, 61.72839506172840, 52.24691358024691]
        self.DtermQ[32]= [0.,0.,0., 304.4628277071887, 246.6507663955451, 192.7893223678521]
        self.DtermQ[31]= [0.,0.,0., 25.43907631106866,-93.99776462634783,-206.8502434238466]
        self.DtermQ[30]= [0.,0.,0.,-260.5342569653863,-305.8365169792609,-342.0455323786138]

    def alf4pi(self,Q):
        return self.alphaS.alfs(Q)/4./np.pi

    def aspi_int_aux(self,mu,poww):
        return 1./mu*self.alf4pi(mu)**poww

    def aspi_int(self,mui,muf,poww):
        intg = lambda mu: self.aspi_int_aux(mu,poww)
        return quad(intg, mui, muf, epsabs=0.00, epsrel=0.05)[0]

    def aspi_log_int_aux(self,mu,Q,poww):
        #print(    mu,self.alf4pi(mu), poww,2.*np.log(Q/mu))
        return 1./mu*self.alf4pi(mu)**poww*2.*np.log(Q/mu)

    def aspi_log_int(self,mui,muf,Q,poww):
        #print(mui,muf)
        intg = lambda mu: self.aspi_log_int_aux(mu,Q,poww)
        return quad(intg, mui, muf, epsabs=0.00, epsrel=0.05)[0]

    def Dterm_q(self,b,mu):
        DtermQ = self.DtermQ
        mc    = self.mc
        mb    = self.mb
        c02   = self.c02
        order = self.logorder

        if mu<mc:
            nf = 3
        elif mu<mb:
            nf = 4
        else:
            nf = 5

        lperp = np.log(mu*mu*b*b/c02)
        alf = self.alf4pi(mu)
        alf2= alf* alf
        alf3= alf2*alf

        res = 0.
        if order >= 2:
            res = alf*(                                                        DtermQ[11][nf]*lperp                 )
        if order >= 3:
            res += alf2*(                           DtermQ[22][nf]*lperp**2. + DtermQ[21][nf]*lperp + DtermQ[20][nf])
        if order >= 4:
            res += alf3*(DtermQ[33][nf]*lperp**3. + DtermQ[32][nf]*lperp**2. + DtermQ[31][nf]*lperp + DtermQ[30][nf])
        return res

    def ExpGamma_q(self,muii,muff):
        mc = self.mc
        mb = self.mb

        NCuspQ = self.NCuspQ
        CuspQ  = self.CuspQ
        order  = self.logorder
        if   muff > muii:
            muf = muff
            mui = muii
        elif muff < muii:
            muf = muii
            mui = muff

        if ((mui < mc) and (muf < mc)):
            if order > 0:
                value = CuspQ[3][1]*self.aspi_log_int(mui,muf,muff,1)
                #print(self.aspi_log_int(1.0,1.2,1.2,1))
                #print(self.aspi_log_int(mui,muf,muff,1))
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[3][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspQ[3][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui >= mc) and (muf >= mc) and (mui < mb)  and (muf < mb)):
            if order > 0:
                value = CuspQ[4][1]*self.aspi_log_int(mui,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[4][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspQ[4][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui >= mb) and (muf >= mb)):
            if order > 0:
                value = CuspQ[5][1]*self.aspi_log_int(mui,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[5][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspQ[5][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui < mc) and (muf >= mc) and (muf < mb)):
            if order > 0:
                value = CuspQ[3][1]*self.aspi_log_int(mui,mc,muff,1)+CuspQ[4][1]*self.aspi_log_int(mc,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[3][i]*self.aspi_log_int(mui,mc,muff,i)+CuspQ[4][i]*self.aspi_log_int(mc,muf,muff,i)+NCuspQ[3][i-1]*self.aspi_int(mui,mc,i-1)+NCuspQ[4][i-1]*self.aspi_int(mc,muf,i-1)

        elif ((mui < mc) and (muf >= mb)):
            if order > 0:
                value = CuspQ[3][1]*self.aspi_log_int(mui,mc,muff,1)+CuspQ[4][1]*self.aspi_log_int(mc,mb,muff,1)+CuspQ[5][1]*self.aspi_log_int(mb,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[3][i]*self.aspi_log_int(mui,mc,muff,i)+CuspQ[4][i]*self.aspi_log_int(mc,mb,muff,i)+CuspQ[5][i]*self.aspi_log_int(mb,muf,muff,i)+NCuspQ[3][i-1]*self.aspi_int(mui,mc,i-1)+NCuspQ[4][i-1]*self.aspi_int(mc,mb,i-1)+NCuspQ[5][i-1]*self.aspi_int(mb,muf,i-1)

        elif ((mui >= mc) and (mui < mb) and (muf >= mb)):
            if order > 0:
                value = CuspQ[4][1]*self.aspi_log_int(mui,mb,muff,1)+CuspQ[5][1]*self.aspi_log_int(mb,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspQ[4][i]*self.aspi_log_int(mui,mb,muff,i)+CuspQ[5][i]*self.aspi_log_int(mb,muf,muff,i)+NCuspQ[4][i-1]*self.aspi_int(mui,mb,i-1)+NCuspQ[5][i-1]*self.aspi_int(mb,muf,i-1)

        if   muff > muii:
            return np.exp(-value)
        elif muff < muii:
            return np.exp( value)
        else:
            return 1.

    def Dterm_g(self,b,mu):
        DtermG = self.DtermG
        mc    = self.mc
        mb    = self.mb
        c02   = self.c02
        order = self.logorder

        if mu<mc:
            nf = 3
        elif mu<mb:
            nf = 4
        else:
            nf = 5

        lperp = np.log(mu*mu*b*b/c02)
        alf = self.alf4pi(mu)
        alf2= alf* alf
        alf3= alf2*alf

        res = 0.
        if order >= 2:
            res = alf*(                                                        DtermG[11][nf]*lperp                 )
        if order >= 3:
            res += alf2*(                           DtermG[22][nf]*lperp**2. + DtermG[21][nf]*lperp + DtermG[20][nf])
        if order >= 4:
            res += alf3*(DtermG[33][nf]*lperp**3. + DtermG[32][nf]*lperp**2. + DtermG[31][nf]*lperp + DtermG[30][nf])
        return res

    def ExpGamma_g(self,muii,muff):
        mc = self.mc
        mb = self.mb

        NCuspG = self.NCuspG
        CuspG  = self.CuspG
        order  = self.logorder
        if   muff > muii:
            muf = muff
            mui = muii
        elif muff < muii:
            muf = muii
            mui = muff

        if ((mui < mc) and (muf < mc)):
            if order > 0:
                value = CuspG[3][1]*self.aspi_log_int(mui,muf,muff,1)
                #print(self.aspi_log_int(1.0,1.2,1.2,1))
                #print(self.aspi_log_int(mui,muf,muff,1))
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[3][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspG[3][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui >= mc) and (muf >= mc) and (mui < mb)  and (muf < mb)):
            if order > 0:
                value = CuspG[4][1]*self.aspi_log_int(mui,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[4][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspG[4][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui >= mb) and (muf >= mb)):
            if order > 0:
                value = CuspG[5][1]*self.aspi_log_int(mui,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[5][i]*self.aspi_log_int(mui,muf,muff,i)+NCuspG[5][i-1]*self.aspi_int(mui,muf,i-1)

        elif ((mui < mc) and (muf >= mc) and (muf < mb)):
            if order > 0:
                value = CuspG[3][1]*self.aspi_log_int(mui,mc,muff,1)+CuspG[4][1]*self.aspi_log_int(mc,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[3][i]*self.aspi_log_int(mui,mc,muff,i)+CuspG[4][i]*self.aspi_log_int(mc,muf,muff,i)+NCuspG[3][i-1]*self.aspi_int(mui,mc,i-1)+NCuspG[4][i-1]*self.aspi_int(mc,muf,i-1)

        elif ((mui < mc) and (muf >= mb)):
            if order > 0:
                value = CuspG[3][1]*self.aspi_log_int(mui,mc,muff,1)+CuspG[4][1]*self.aspi_log_int(mc,mb,muff,1)+CuspG[5][1]*self.aspi_log_int(mb,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[3][i]*self.aspi_log_int(mui,mc,muff,i)+CuspG[4][i]*self.aspi_log_int(mc,mb,muff,i)+CuspG[5][i]*self.aspi_log_int(mb,muf,muff,i)+NCuspG[3][i-1]*self.aspi_int(mui,mc,i-1)+NCuspG[4][i-1]*self.aspi_int(mc,mb,i-1)+NCuspG[5][i-1]*self.aspi_int(mb,muf,i-1)

        elif ((mui >= mc) and (mui < mb) and (muf >= mb)):
            if order > 0:
                value = CuspG[4][1]*self.aspi_log_int(mui,mb,muff,1)+CuspG[5][1]*self.aspi_log_int(mb,muf,muff,1)
            if order > 1:
                for i in range(2,order+1):
                    value += CuspG[4][i]*self.aspi_log_int(mui,mb,muff,i)+CuspG[5][i]*self.aspi_log_int(mb,muf,muff,i)+NCuspG[4][i-1]*self.aspi_int(mui,mb,i-1)+NCuspG[5][i-1]*self.aspi_int(mb,muf,i-1)

        if   muff > muii:
            return np.exp(-value)
        elif muff < muii:
            return np.exp( value)
        else:
            return 1.

    def anomq(self,b,mui,muf,Qi,Qf):
        res1 = self.ExpGamma_q(mui,muf)
        res2 = self.Dterm_q(b,mui)
        return res1*np.exp(-res2*2.*np.log(Qf/Qi))

    def anomg(self,b,mui,muf,Qi,Qf):
        res1 = self.ExpGamma_g(mui,muf)
        res2 = self.Dterm_g(b,mui)
        return res1*np.exp(-res2*2.*np.log(Qf/Qi))


