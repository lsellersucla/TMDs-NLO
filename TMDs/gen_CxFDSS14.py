from TMD import TMDPDF,TMDFF
import numpy as np
from scipy.integrate import quad

scheme = "MSbar"

order    = 0
xcut     = 0.999
tmdlo    = TMDFF(order,xcut,"DSS14PI",0,1e-3,scheme)

order    = 1
tmdnlo   = TMDFF(order,xcut,"DSS14PI",1,1e-3,scheme)

Lmu      = 0.
Lz       = 0.

xx=np.array([0.100000e-01,0.200000e-01,0.300000e-01,0.400000e-01,0.500000e-01,0.600000e-01,0.700000e-01,0.800000e-01,0.900000e-01,0.950000e-01,0.100000e+00,0.125000e+00,0.150000e+00,0.175000e+00,0.200000e+00,0.225000e+00,0.250000e+00,0.275000e+00,0.300000e+00,0.325000e+00,0.350000e+00,0.375000e+00,0.400000e+00,0.450000e+00,0.500000e+00,0.550000e+00,0.600000e+00,0.650000e+00,0.700000e+00,0.750000e+00,0.800000e+00,0.850000e+00,0.900000e+00,0.930000e+00,0.100000e+01])
QQ=np.array([0.100000e+01,0.111803e+01,0.122474e+01,0.158114e+01,0.200000e+01,0.252982e+01,0.316228e+01,0.387298e+01,0.500000e+01,0.632456e+01,0.800000e+01,0.141421e+02,0.134164e+02,0.178885e+02,0.240832e+02,0.316228e+02,0.424264e+02,0.979796e+02,0.761577e+02,0.100000e+03,0.134164e+03,0.178885e+03,0.240832e+03,0.316228e+03])

f = open('CxFDSS14pinlo/CxFDSS14pinlo_0000.dat', 'w')

#f = open("data.dat","w")

f.write("PdfType: central\n")
f.write("Format: lhagrid1\n")
f.write("---\n")

stt = ""
for i in range(len(xx)-1):
    x = str(np.format_float_scientific(xx[i], precision = 6, exp_digits=2)).replace("e","e")
    while len(x)<12:
        x = x.replace("e","0e")
    deli = " "
    stt = stt+x+deli
x = str(np.format_float_scientific(xx[-1], precision = 6, exp_digits=2)).replace("e","e")
while len(x)<12:
    x = x.replace("e","0e")
stt = stt+x+"\n"
f.write(stt)

stt = ""
for i in range(len(QQ)-1):
    Q = str(np.format_float_scientific(QQ[i], precision = 6, exp_digits=2)).replace("e","e")
    while len(Q)<12:
        Q = Q.replace("e","0e")
    deli = " "
    stt = stt+Q+deli
x = str(np.format_float_scientific(QQ[-1], precision = 6, exp_digits=2)).replace("e","e")
while len(Q)<12:
    Q = Q.replace("e","0e")
stt = stt+Q+"\n"
f.write(stt)

f.write("-5 -4 -3 -2 -1 1 2 3 4 5 21\n")

for i in range(len(xx)):
    x = xx[i]
    for j in range(len(QQ)):
        Q = QQ[j]
        vals = x*tmdnlo.CxFF(x,Q,Lmu,Lz)

        if vals[-5]>=0:
            vm5 = str(np.format_float_scientific(vals[-5], precision = 6, exp_digits=2))
        else:
            vm5 = str(np.format_float_scientific(vals[-5], precision = 5, exp_digits=2))
        if vals[-4]>=0:
            vm4 = str(np.format_float_scientific(vals[-4], precision = 6, exp_digits=2))
        else:
            vm4 = str(np.format_float_scientific(vals[-4], precision = 5, exp_digits=2))
        if vals[-3]>=0:
            vm3 = str(np.format_float_scientific(vals[-3], precision = 6, exp_digits=2))
        else:
            vm3 = str(np.format_float_scientific(vals[-3], precision = 5, exp_digits=2))
        if vals[-2]>=0:
            vm2 = str(np.format_float_scientific(vals[-2], precision = 6, exp_digits=2))
        else:
            vm2 = str(np.format_float_scientific(vals[-2], precision = 5, exp_digits=2))
        if vals[-1]>=0:
            vm1 = str(np.format_float_scientific(vals[-1], precision = 6, exp_digits=2))
        else:
            vm1 = str(np.format_float_scientific(vals[-1], precision = 5, exp_digits=2))
        if vals[1]>=0:
            vp1 = str(np.format_float_scientific(vals[ 1], precision = 6, exp_digits=2))
        else:
            vp1 = str(np.format_float_scientific(vals[ 1], precision = 5, exp_digits=2))
        if vals[2]>=0:
            vp2 = str(np.format_float_scientific(vals[ 2], precision = 6, exp_digits=2))
        else:
            vp2 = str(np.format_float_scientific(vals[ 2], precision = 5, exp_digits=2))
        if vals[3]>=0:
            vp3 = str(np.format_float_scientific(vals[ 3], precision = 6, exp_digits=2))
        else:
            vp3 = str(np.format_float_scientific(vals[ 3], precision = 5, exp_digits=2))
        if vals[4]>=0:
            vp4 = str(np.format_float_scientific(vals[ 4], precision = 6, exp_digits=2))
        else:
            vp4 = str(np.format_float_scientific(vals[ 4], precision = 5, exp_digits=2))
        if vals[5]>=0:
            vp5 = str(np.format_float_scientific(vals[ 5], precision = 6, exp_digits=2))
        else:
            vp5 = str(np.format_float_scientific(vals[ 5], precision = 5, exp_digits=2))
        if vals[0]>=0:
            v00 = str(np.format_float_scientific(vals[ 0], precision = 6, exp_digits=2))
        else:
            v00 = str(np.format_float_scientific(vals[ 0], precision = 5, exp_digits=2))

        while len(vm5)<12:
            vm5 = vm5.replace("e","0e")
        while len(vm4)<12:
            vm4 = vm4.replace("e","0e")
        while len(vm3)<12:
            vm3 = vm3.replace("e","0e")
        while len(vm2)<12:
            vm2 = vm2.replace("e","0e")
        while len(vm1)<12:
            vm1 = vm1.replace("e","0e")
        while len(vp1)<12:
            vp1 = vp1.replace("e","0e")
        while len(vp2)<12:
            vp2 = vp2.replace("e","0e")
        while len(vp3)<12:
            vp3 = vp3.replace("e","0e")
        while len(vp4)<12:
            vp4 = vp4.replace("e","0e")
        while len(vp5)<12:
            vp5 = vp5.replace("e","0e")
        while len(v00)<12:
            v00 = v00.replace("e","0e")
        deli= " "
        output = deli+vm5+deli+vm4+deli+vm3+deli+vm2+deli+vm1+deli+vp1+deli+vp2+deli+vp3+deli+vp4+deli+vp5+deli+v00+"\n"
        output = output.replace("e","e")
        f.write(output)

f.write("---\n")
f.close()


f = open('CxFDSS14pinlo/CxFDSS14pinlo_0001.dat', 'w')

#f = open("data.dat","w")

f.write("PdfType: error\n")
f.write("Format: lhagrid1\n")
f.write("---\n")

stt = ""
for i in range(len(xx)-1):
    x = str(np.format_float_scientific(xx[i], precision = 6, exp_digits=2)).replace("e","e")
    while len(x)<12:
        x = x.replace("e","0e")
    deli = " "
    stt = stt+x+deli
x = str(np.format_float_scientific(xx[-1], precision = 6, exp_digits=2)).replace("e","e")
while len(x)<12:
    x = x.replace("e","0e")
stt = stt+x+"\n"
f.write(stt)

stt = ""
for i in range(len(QQ)-1):
    Q = str(np.format_float_scientific(QQ[i], precision = 6, exp_digits=2)).replace("e","e")
    while len(Q)<12:
        Q = Q.replace("e","0e")
    deli = " "
    stt = stt+Q+deli
x = str(np.format_float_scientific(QQ[-1], precision = 6, exp_digits=2)).replace("e","e")
while len(Q)<12:
    Q = Q.replace("e","0e")
stt = stt+Q+"\n"
f.write(stt)

f.write("-5 -4 -3 -2 -1 1 2 3 4 5 21\n")

for i in range(len(xx)):
    x = xx[i]
    for j in range(len(QQ)):
        Q = QQ[j]
        vals = x*tmdnlo.CxFF(x,Q,Lmu,Lz)

        if vals[-5]>=0:
            vm5 = str(np.format_float_scientific(vals[-5], precision = 6, exp_digits=2))
        else:
            vm5 = str(np.format_float_scientific(vals[-5], precision = 5, exp_digits=2))
        if vals[-4]>=0:
            vm4 = str(np.format_float_scientific(vals[-4], precision = 6, exp_digits=2))
        else:
            vm4 = str(np.format_float_scientific(vals[-4], precision = 5, exp_digits=2))
        if vals[-3]>=0:
            vm3 = str(np.format_float_scientific(vals[-3], precision = 6, exp_digits=2))
        else:
            vm3 = str(np.format_float_scientific(vals[-3], precision = 5, exp_digits=2))
        if vals[-2]>=0:
            vm2 = str(np.format_float_scientific(vals[-2], precision = 6, exp_digits=2))
        else:
            vm2 = str(np.format_float_scientific(vals[-2], precision = 5, exp_digits=2))
        if vals[-1]>=0:
            vm1 = str(np.format_float_scientific(vals[-1], precision = 6, exp_digits=2))
        else:
            vm1 = str(np.format_float_scientific(vals[-1], precision = 5, exp_digits=2))
        if vals[1]>=0:
            vp1 = str(np.format_float_scientific(vals[ 1], precision = 6, exp_digits=2))
        else:
            vp1 = str(np.format_float_scientific(vals[ 1], precision = 5, exp_digits=2))
        if vals[2]>=0:
            vp2 = str(np.format_float_scientific(vals[ 2], precision = 6, exp_digits=2))
        else:
            vp2 = str(np.format_float_scientific(vals[ 2], precision = 5, exp_digits=2))
        if vals[3]>=0:
            vp3 = str(np.format_float_scientific(vals[ 3], precision = 6, exp_digits=2))
        else:
            vp3 = str(np.format_float_scientific(vals[ 3], precision = 5, exp_digits=2))
        if vals[4]>=0:
            vp4 = str(np.format_float_scientific(vals[ 4], precision = 6, exp_digits=2))
        else:
            vp4 = str(np.format_float_scientific(vals[ 4], precision = 5, exp_digits=2))
        if vals[5]>=0:
            vp5 = str(np.format_float_scientific(vals[ 5], precision = 6, exp_digits=2))
        else:
            vp5 = str(np.format_float_scientific(vals[ 5], precision = 5, exp_digits=2))
        if vals[0]>=0:
            v00 = str(np.format_float_scientific(vals[ 0], precision = 6, exp_digits=2))
        else:
            v00 = str(np.format_float_scientific(vals[ 0], precision = 5, exp_digits=2))

        while len(vm5)<12:
            vm5 = vm5.replace("e","0e")
        while len(vm4)<12:
            vm4 = vm4.replace("e","0e")
        while len(vm3)<12:
            vm3 = vm3.replace("e","0e")
        while len(vm2)<12:
            vm2 = vm2.replace("e","0e")
        while len(vm1)<12:
            vm1 = vm1.replace("e","0e")
        while len(vp1)<12:
            vp1 = vp1.replace("e","0e")
        while len(vp2)<12:
            vp2 = vp2.replace("e","0e")
        while len(vp3)<12:
            vp3 = vp3.replace("e","0e")
        while len(vp4)<12:
            vp4 = vp4.replace("e","0e")
        while len(vp5)<12:
            vp5 = vp5.replace("e","0e")
        while len(v00)<12:
            v00 = v00.replace("e","0e")
        deli= " "
        output = deli+vm5+deli+vm4+deli+vm3+deli+vm2+deli+vm1+deli+vp1+deli+vp2+deli+vp3+deli+vp4+deli+vp5+deli+v00+"\n"
        output = output.replace("e","e")
        f.write(output)

f.write("---\n")
f.close()



