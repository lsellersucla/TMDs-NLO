from TMD import TMDPDF
import numpy as np

order    = 0
xcut     = 0.999
scheme   = "MSbar"
tmdlo    = TMDPDF(order,xcut,"CT14nlo",0,1e-3, scheme)

order    = 1
tmdnlo   = TMDPDF(order,xcut,"CT14nlo",0,1e-3, scheme)

Lmu      = 0.
Lz       = 0.

xx=np.array([0.100000e-06, 0.106977e-06, 0.114440e-06, 0.122424e-06, 0.130965e-06, 0.140102e-06, 0.149876e-06, 0.160332e-06, 0.171518e-06, 0.183484e-06, 0.196285e-06, 0.209979e-06, 0.224628e-06, 0.240300e-06, 0.257064e-06, 0.274999e-06, 0.294184e-06, 0.314708e-06, 0.336664e-06, 0.360152e-06, 0.385278e-06, 0.412157e-06, 0.440912e-06, 0.471672e-06, 0.504579e-06, 0.539781e-06, 0.577439e-06, 0.617725e-06, 0.660821e-06, 0.706923e-06, 0.756243e-06, 0.809002e-06, 0.865443e-06, 0.925821e-06, 0.990412e-06, 0.105951e-05, 0.113343e-05, 0.121250e-05, 0.129709e-05, 0.138758e-05, 0.148439e-05, 0.158795e-05, 0.169873e-05, 0.181725e-05, 0.194403e-05, 0.207966e-05, 0.222475e-05, 0.237996e-05, 0.254600e-05, 0.272362e-05, 0.291363e-05, 0.311691e-05, 0.333436e-05, 0.356698e-05, 0.381584e-05, 0.408205e-05, 0.436684e-05, 0.467150e-05, 0.499741e-05, 0.534606e-05, 0.571903e-05, 0.611802e-05, 0.654485e-05, 0.700145e-05, 0.748992e-05, 0.801246e-05, 0.857145e-05, 0.916945e-05, 0.980916e-05, 0.104935e-04, 0.112256e-04, 0.120088e-04, 0.128466e-04, 0.137428e-04, 0.147016e-04, 0.157272e-04, 0.168245e-04, 0.179982e-04, 0.192539e-04, 0.205972e-04, 0.220341e-04, 0.235714e-04, 0.252159e-04, 0.269751e-04, 0.288570e-04, 0.308702e-04, 0.330239e-04, 0.353278e-04, 0.377925e-04, 0.404291e-04, 0.432497e-04, 0.462671e-04, 0.494949e-04, 0.529480e-04, 0.566419e-04, 0.605936e-04, 0.648210e-04, 0.693433e-04, 0.741810e-04, 0.793563e-04, 0.848927e-04, 0.908153e-04, 0.971511e-04, 0.103929e-03, 0.111180e-03, 0.118936e-03, 0.127234e-03, 0.136110e-03, 0.145606e-03, 0.155765e-03, 0.166632e-03, 0.178257e-03, 0.190693e-03, 0.203997e-03, 0.218229e-03, 0.233454e-03, 0.249741e-03, 0.267164e-03, 0.285803e-03, 0.305742e-03, 0.327073e-03, 0.349891e-03, 0.374302e-03, 0.400415e-03, 0.428350e-03, 0.458235e-03, 0.490204e-03, 0.524403e-03, 0.560988e-03, 0.600126e-03, 0.641995e-03, 0.686784e-03, 0.734698e-03, 0.785955e-03, 0.840787e-03, 0.899446e-03, 0.962196e-03, 0.102932e-02, 0.110114e-02, 0.117796e-02, 0.126014e-02, 0.134805e-02, 0.144210e-02, 0.154271e-02, 0.165034e-02, 0.176548e-02, 0.188865e-02, 0.202041e-02, 0.216136e-02, 0.231215e-02, 0.247346e-02, 0.264603e-02, 0.283063e-02, 0.302811e-02, 0.323937e-02, 0.346536e-02, 0.370713e-02, 0.396576e-02, 0.424243e-02, 0.453841e-02, 0.485504e-02, 0.519375e-02, 0.555610e-02, 0.594372e-02, 0.635839e-02, 0.680199e-02, 0.727654e-02, 0.778419e-02, 0.832726e-02, 0.890822e-02, 0.952971e-02, 0.101946e-01, 0.109058e-01, 0.116666e-01, 0.124806e-01, 0.133513e-01, 0.142827e-01, 0.152792e-01, 0.163452e-01, 0.174855e-01, 0.187054e-01, 0.200104e-01, 0.214064e-01, 0.228998e-01, 0.244975e-01, 0.262066e-01, 0.280349e-01, 0.299908e-01, 0.320831e-01, 0.343214e-01, 0.367158e-01, 0.392774e-01, 0.420176e-01, 0.449490e-01, 0.480849e-01, 0.514395e-01, 0.550283e-01, 0.588674e-01, 0.629743e-01, 0.673677e-01, 0.720677e-01, 0.770955e-01, 0.824742e-01, 0.882281e-01, 0.943834e-01, 0.100968e+00, 0.108012e+00, 0.115548e+00, 0.123609e+00, 0.132233e+00, 0.141458e+00, 0.151327e+00, 0.161884e+00, 0.173178e+00, 0.185260e+00, 0.198185e+00, 0.212012e+00, 0.226803e+00, 0.242626e+00, 0.259553e+00, 0.277661e+00, 0.297032e+00, 0.317755e+00, 0.339923e+00, 0.363638e+00, 0.389008e+00, 0.416147e+00, 0.445180e+00, 0.476238e+00, 0.509463e+00, 0.545007e+00, 0.583029e+00, 0.623705e+00, 0.667218e+00, 0.713767e+00, 0.763564e+00, 0.816834e+00, 0.873821e+00, 0.934784e+00, 0.100000e+01])
QQ=np.array([0.130000e+01, 0.166682e+01, 0.213715e+01, 0.274019e+01, 0.351340e+01, 0.450478e+01, 0.577590e+01, 0.740569e+01, 0.949536e+01, 0.121747e+02, 0.156100e+02, 0.200147e+02, 0.256623e+02, 0.329035e+02, 0.421879e+02, 0.540921e+02, 0.693553e+02, 0.889254e+02, 0.114018e+03, 0.146190e+03, 0.187441e+03, 0.240331e+03, 0.308145e+03, 0.395095e+03, 0.506580e+03, 0.649522e+03, 0.832798e+03, 0.106779e+04, 0.136909e+04, 0.175541e+04, 0.225073e+04, 0.288582e+04, 0.370012e+04, 0.474419e+04, 0.608286e+04, 0.779927e+04, 0.100000e+05])

f = open('CxFCT14nlo/CxFCT14nlo_0000.dat', 'w')

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
        vals = x*tmdnlo.CxPDF(x,Q,Lmu,Lz)

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
