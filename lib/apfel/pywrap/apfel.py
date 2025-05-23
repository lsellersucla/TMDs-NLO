# This file was automatically generated by SWIG (https://www.swig.org).
# Version 4.3.0
#
# Do not make changes to this file unless you know what you are doing - modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _apfel
else:
    import _apfel

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "this":
            set(self, name, value)
        elif name == "thisown":
            self.this.own(value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _apfel.delete_SwigPyIterator

    def value(self):
        return _apfel.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _apfel.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _apfel.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _apfel.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _apfel.SwigPyIterator_equal(self, x)

    def copy(self):
        return _apfel.SwigPyIterator_copy(self)

    def next(self):
        return _apfel.SwigPyIterator_next(self)

    def __next__(self):
        return _apfel.SwigPyIterator___next__(self)

    def previous(self):
        return _apfel.SwigPyIterator_previous(self)

    def advance(self, n):
        return _apfel.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _apfel.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _apfel.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _apfel.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _apfel.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _apfel.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _apfel.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _apfel:
_apfel.SwigPyIterator_swigregister(SwigPyIterator)

def new_doubles(nelements):
    return _apfel.new_doubles(nelements)

def delete_doubles(ary):
    return _apfel.delete_doubles(ary)

def doubles_getitem(ary, index):
    return _apfel.doubles_getitem(ary, index)

def doubles_setitem(ary, index, value):
    return _apfel.doubles_setitem(ary, index, value)

def InitializeAPFEL():
    return _apfel.InitializeAPFEL()

def EvolveAPFEL(Q0, Q):
    return _apfel.EvolveAPFEL(Q0, Q)

def DeriveAPFEL(Q):
    return _apfel.DeriveAPFEL(Q)

def CachePDFsAPFEL(Q0):
    return _apfel.CachePDFsAPFEL(Q0)

def xPDF(i, x):
    return _apfel.xPDF(i, x)

def xPDFxQ(i, x, Q):
    return _apfel.xPDFxQ(i, x, Q)

def dxPDF(i, x):
    return _apfel.dxPDF(i, x)

def xPDFj(i, x):
    return _apfel.xPDFj(i, x)

def xgamma(x):
    return _apfel.xgamma(x)

def xgammaj(x):
    return _apfel.xgammaj(x)

def dxgamma(x):
    return _apfel.dxgamma(x)

def xPDFall(x, xf):
    return _apfel.xPDFall(x, xf)

def xPDFallPhoton(x, xf):
    return _apfel.xPDFallPhoton(x, xf)

def xPDFxQall(x, Q, xf):
    return _apfel.xPDFxQall(x, Q, xf)

def xLepton(i, x):
    return _apfel.xLepton(i, x)

def xLeptonj(i, x):
    return _apfel.xLeptonj(i, x)

def ExternalEvolutionOperator(fname, i, j, x, beta):
    return _apfel.ExternalEvolutionOperator(fname, i, j, x, beta)

def ExternalEvolutionMatrixEv2Ev(i, j, alpha, beta):
    return _apfel.ExternalEvolutionMatrixEv2Ev(i, j, alpha, beta)

def ExternalEvolutionMatrixEv2Ph(i, j, alpha, beta):
    return _apfel.ExternalEvolutionMatrixEv2Ph(i, j, alpha, beta)

def ExternalEvolutionMatrixPh2Ph(i, j, alpha, beta):
    return _apfel.ExternalEvolutionMatrixPh2Ph(i, j, alpha, beta)

def ComputeExternalSplittingFunctions(fname, pt, nf, x, beta):
    return _apfel.ComputeExternalSplittingFunctions(fname, pt, nf, x, beta)

def ExternalSplittingFunctions(i, j):
    return _apfel.ExternalSplittingFunctions(i, j)

def LHAPDFgrid(Nrep, Qin, fname):
    return _apfel.LHAPDFgrid(Nrep, Qin, fname)

def LHAPDFgridDerivative(Nrep, fname):
    return _apfel.LHAPDFgridDerivative(Nrep, fname)

def AlphaQCD(Q):
    return _apfel.AlphaQCD(Q)

def AlphaQED(Q):
    return _apfel.AlphaQED(Q)

def HeavyQuarkMass(arg1, arg2):
    return _apfel.HeavyQuarkMass(arg1, arg2)

def GetThreshold(arg1):
    return _apfel.GetThreshold(arg1)

def GetMaxFlavourAlpha():
    return _apfel.GetMaxFlavourAlpha()

def GetMaxFlavourPDFs():
    return _apfel.GetMaxFlavourPDFs()

def HeavyQuarkThreshold(arg1):
    return _apfel.HeavyQuarkThreshold(arg1)

def NPDF(i, N):
    return _apfel.NPDF(i, N)

def Ngamma(N):
    return _apfel.Ngamma(N)

def LUMI(i, j, S):
    return _apfel.LUMI(i, j, S)

def xGrid(alpha):
    return _apfel.xGrid(alpha)

def nIntervals():
    return _apfel.nIntervals()

def GetVersion():
    return _apfel.GetVersion()

def CleanUp():
    return _apfel.CleanUp()

def EnableWelcomeMessage(arg1):
    return _apfel.EnableWelcomeMessage(arg1)

def EnableEvolutionOperator(arg1):
    return _apfel.EnableEvolutionOperator(arg1)

def EnableLeptonEvolution(arg1):
    return _apfel.EnableLeptonEvolution(arg1)

def LockGrids(arg1):
    return _apfel.LockGrids(arg1)

def SetTimeLikeEvolution(arg1):
    return _apfel.SetTimeLikeEvolution(arg1)

def SetPolarizedEvolution(arg1):
    return _apfel.SetPolarizedEvolution(arg1)

def SetFastEvolution(arg1):
    return _apfel.SetFastEvolution(arg1)

def EnableMassRunning(arg1):
    return _apfel.EnableMassRunning(arg1)

def SetSmallxResummation(arg1, la):
    return _apfel.SetSmallxResummation(arg1, la)

def SetAlphaQCDRef(alpharef, Qref):
    return _apfel.SetAlphaQCDRef(alpharef, Qref)

def SetAlphaQEDRef(alpharef, Qref):
    return _apfel.SetAlphaQEDRef(alpharef, Qref)

def SetAlphaEvolution(evol):
    return _apfel.SetAlphaEvolution(evol)

def SetLambdaQCDRef(lambdaref, nref):
    return _apfel.SetLambdaQCDRef(lambdaref, nref)

def SetPDFEvolution(evolp):
    return _apfel.SetPDFEvolution(evolp)

def SetEpsilonTruncation(eps):
    return _apfel.SetEpsilonTruncation(eps)

def SetQLimits(Qmin, Qmax):
    return _apfel.SetQLimits(Qmin, Qmax)

def SetFFNS(nfl):
    return _apfel.SetFFNS(nfl)

def SetGridParameters(i, np, deg, x):
    return _apfel.SetGridParameters(i, np, deg, x)

def SetQGridParameters(npq, degq):
    return _apfel.SetQGridParameters(npq, degq)

def SetLHgridParameters(nx, nxm, xmin, xm, xmax, nq2, q2min, q2max):
    return _apfel.SetLHgridParameters(nx, nxm, xmin, xm, xmax, nq2, q2min, q2max)

def SetExternalGrid(i, np, deg, x):
    return _apfel.SetExternalGrid(i, np, deg, x)

def SetMaxFlavourAlpha(nf):
    return _apfel.SetMaxFlavourAlpha(nf)

def SetMaxFlavourPDFs(nf):
    return _apfel.SetMaxFlavourPDFs(nf)

def SetMSbarMasses(mc, mb, mt):
    return _apfel.SetMSbarMasses(mc, mb, mt)

def SetMassScaleReference(Qc, Qb, Qt):
    return _apfel.SetMassScaleReference(Qc, Qb, Qt)

def SetMassMatchingScales(kmc, kmb, kmt):
    return _apfel.SetMassMatchingScales(kmc, kmb, kmt)

def SetNumberOfGrids(n):
    return _apfel.SetNumberOfGrids(n)

def SetPDFSet(name):
    return _apfel.SetPDFSet(name)

def SetPerturbativeOrder(pto):
    return _apfel.SetPerturbativeOrder(pto)

def GetPerturbativeOrder():
    return _apfel.GetPerturbativeOrder()

def GetMuF():
    return _apfel.GetMuF()

def GetMuF0():
    return _apfel.GetMuF0()

def SetPoleMasses(mc, mb, mt):
    return _apfel.SetPoleMasses(mc, mb, mt)

def SetTauMass(masst):
    return _apfel.SetTauMass(masst)

def SetRenFacRatio(ratio):
    return _apfel.SetRenFacRatio(ratio)

def SetReplica(nr):
    return _apfel.SetReplica(nr)

def SetTheory(theory):
    return _apfel.SetTheory(theory)

def EnableNLOQEDCorrections(arg1):
    return _apfel.EnableNLOQEDCorrections(arg1)

def SetVFNS():
    return _apfel.SetVFNS()

def ListFunctions():
    return _apfel.ListFunctions()

def CheckAPFEL():
    return _apfel.CheckAPFEL()

def InitializeAPFEL_DIS():
    return _apfel.InitializeAPFEL_DIS()

def ComputeStructureFunctionsAPFEL(Q0, Q):
    return _apfel.ComputeStructureFunctionsAPFEL(Q0, Q)

def CacheStructureFunctionsAPFEL(Q0):
    return _apfel.CacheStructureFunctionsAPFEL(Q0)

def SetMassScheme(ms):
    return _apfel.SetMassScheme(ms)

def SetPolarizationDIS(pol):
    return _apfel.SetPolarizationDIS(pol)

def SetProcessDIS(pr):
    return _apfel.SetProcessDIS(pr)

def SetProjectileDIS(lept):
    return _apfel.SetProjectileDIS(lept)

def SetTargetDIS(tar):
    return _apfel.SetTargetDIS(tar)

def SelectCharge(selch):
    return _apfel.SelectCharge(selch)

def ExternalDISOperator(SF, ihq, i, x, beta):
    return _apfel.ExternalDISOperator(SF, ihq, i, x, beta)

def F2light(x):
    return _apfel.F2light(x)

def F2charm(x):
    return _apfel.F2charm(x)

def F2bottom(x):
    return _apfel.F2bottom(x)

def F2top(x):
    return _apfel.F2top(x)

def F2total(x):
    return _apfel.F2total(x)

def FLlight(x):
    return _apfel.FLlight(x)

def FLcharm(x):
    return _apfel.FLcharm(x)

def FLbottom(x):
    return _apfel.FLbottom(x)

def FLtop(x):
    return _apfel.FLtop(x)

def FLtotal(x):
    return _apfel.FLtotal(x)

def F3light(x):
    return _apfel.F3light(x)

def F3charm(x):
    return _apfel.F3charm(x)

def F3bottom(x):
    return _apfel.F3bottom(x)

def F3top(x):
    return _apfel.F3top(x)

def F3total(x):
    return _apfel.F3total(x)

def g1light(x):
    return _apfel.g1light(x)

def g1charm(x):
    return _apfel.g1charm(x)

def g1bottom(x):
    return _apfel.g1bottom(x)

def g1top(x):
    return _apfel.g1top(x)

def g1total(x):
    return _apfel.g1total(x)

def gLlight(x):
    return _apfel.gLlight(x)

def gLcharm(x):
    return _apfel.gLcharm(x)

def gLbottom(x):
    return _apfel.gLbottom(x)

def gLtop(x):
    return _apfel.gLtop(x)

def gLtotal(x):
    return _apfel.gLtotal(x)

def g4light(x):
    return _apfel.g4light(x)

def g4charm(x):
    return _apfel.g4charm(x)

def g4bottom(x):
    return _apfel.g4bottom(x)

def g4top(x):
    return _apfel.g4top(x)

def g4total(x):
    return _apfel.g4total(x)

def StructureFunctionxQ(proc, sf, comp, x, Q):
    return _apfel.StructureFunctionxQ(proc, sf, comp, x, Q)

def SetZMass(massz):
    return _apfel.SetZMass(massz)

def SetWMass(massw):
    return _apfel.SetWMass(massw)

def SetProtonMass(massp):
    return _apfel.SetProtonMass(massp)

def SetSin2ThetaW(sw):
    return _apfel.SetSin2ThetaW(sw)

def SetCKM(vud, vus, vub, vcd, vcs, vcb, vtd, vts, vtb):
    return _apfel.SetCKM(vud, vus, vub, vcd, vcs, vcb, vtd, vts, vtb)

def SetPropagatorCorrection(dr):
    return _apfel.SetPropagatorCorrection(dr)

def SetEWCouplings(vd, vu, ad, au):
    return _apfel.SetEWCouplings(vd, vu, ad, au)

def SetGFermi(gf):
    return _apfel.SetGFermi(gf)

def SetRenQRatio(ratioR):
    return _apfel.SetRenQRatio(ratioR)

def SetFacQRatio(ratioF):
    return _apfel.SetFacQRatio(ratioF)

def EnableDynamicalScaleVariations(arg1):
    return _apfel.EnableDynamicalScaleVariations(arg1)

def EnableIntrinsicCharm(arg1):
    return _apfel.EnableIntrinsicCharm(arg1)

def GetZMass():
    return _apfel.GetZMass()

def GetWMass():
    return _apfel.GetWMass()

def GetProtonMass():
    return _apfel.GetProtonMass()

def GetSin2ThetaW():
    return _apfel.GetSin2ThetaW()

def GetCKM(u, d):
    return _apfel.GetCKM(u, d)

def GetGFermi():
    return _apfel.GetGFermi()

def GetSIATotalCrossSection(pto, q, comp):
    return _apfel.GetSIATotalCrossSection(pto, q, comp)

def EnableTargetMassCorrections(arg1):
    return _apfel.EnableTargetMassCorrections(arg1)

def EnableDampingFONLL(arg1):
    return _apfel.EnableDampingFONLL(arg1)

def SetDampingPowerFONLL(arg1, arg2, arg3):
    return _apfel.SetDampingPowerFONLL(arg1, arg2, arg3)

def ComputeChargesDIS(q2, bq, dq, bqt):
    return _apfel.ComputeChargesDIS(q2, bq, dq, bqt)

def F2LO(x, q):
    return _apfel.F2LO(x, q)

def FKSimulator(x, q, y, i, beta):
    return _apfel.FKSimulator(x, q, y, i, beta)

def SetFKObservable(obs):
    return _apfel.SetFKObservable(obs)

def GetFKObservable():
    return _apfel.GetFKObservable()

def FKObservables(x, q, y):
    return _apfel.FKObservables(x, q, y)

def ComputeFKTables(inputfile, outputpath, Q0, flmap):
    return _apfel.ComputeFKTables(inputfile, outputpath, Q0, flmap)

def ComputeHardCrossSectionsDY(datafile, outputfile):
    return _apfel.ComputeHardCrossSectionsDY(datafile, outputfile)

def EnableSFNLOQEDCorrections(arg1):
    return _apfel.EnableSFNLOQEDCorrections(arg1)

def LHAPDFgridStructureFunctions(Nrep, Qin, fname):
    return _apfel.LHAPDFgridStructureFunctions(Nrep, Qin, fname)

def SetScaleVariationProcedure(svp):
    return _apfel.SetScaleVariationProcedure(svp)

