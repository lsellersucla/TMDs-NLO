#!/usr/bin/env python
import numpy as np
import sys
import time

try:
    from Collins_fit.kits.kits import check_path
    from Collins_fit.TMDs.TMDs.TMD import TMDPDF, TMDFF
    from Collins_fit.TMDs.Evolve.alpha import ALPHAS
except:
    from TMD import TMDPDF, TMDFF
    from alpha import ALPHAS

def generate_info(choice, pdf_set, scheme):
    Qs = np.logspace(np.log10(1.295), 5, num = 37, base = 10.0)
    name = './base.info'
    with open(name, 'r') as _:
        lines = _.readlines()

    if choice == 'pdf':
        particle = 2212
    elif choice == 'ff':
        if 'pi' in pdf_set.lower():
            particle = 211
        elif 'ka' in pdf_set.lower():
            particle = 321
    x_min = 0.05
    Q_min = 1.0
    alpha_s = ALPHAS(1)

    ## get alpha_s values
    Qs_txt = ''
    alpha_s_txt = ''
    for i in range(len(Qs) - 1):
        Qs_txt      += '{:.6e}, '.format(Qs[i])
        alpha_s_txt += '{:.6e}, '.format(alpha_s.alfs(Qs[i]))
    Qs_txt      += '{:.6e}'.format(Qs[-1])
    alpha_s_txt += '{:.6e}'.format(alpha_s.alfs(Qs[-1]))

    ## make substitutions
    lines[0 ] = lines[0 ].replace('<description>', 'C times {}s'.format(choice.upper()))
    lines[2 ] = lines[2 ].replace('<authors>'    , 'UCLA')
    lines[3 ] = lines[3 ].replace('<reference>'  , 'None')
    lines[6 ] = lines[6 ].replace('<n_member>'   , '2')
    lines[7 ] = lines[7 ].replace('<particle>'   , str(particle))
    lines[15] = lines[15].replace('<x_min>'      , str(x_min))
    lines[17] = lines[17].replace('<Q_min>'      , str(Q_min))
    lines[29] = lines[29].replace('<alpha_s_Qs>' , Qs_txt)
    lines[30] = lines[30].replace('<alpha_s>'    , alpha_s_txt)

    ## save info file
    name = f'Cx{pdf_set}-{scheme}/Cx{pdf_set}-{scheme}.info'
    with open(name, 'w') as _:
        _.writelines(lines)
    return

def generate_output(choice, theory, xs, Qs, order, pdf_set, scheme, Lmu = 0.0, Lz = 0.0):
    name = f'Cx{pdf_set}-{scheme}/Cx{pdf_set}-{scheme}'
    check_path(name[: name.rindex('/') + 1])
    print('generating ', name[name.rindex('/') + 1 :])
    if order == 0:
        name += f'_0000.dat'
    elif order == 1:
        name += f'_0001.dat'
    f = open(name, 'w')

    f.write("PdfType: central\n")
    f.write("Format: lhagrid1\n")
    f.write("---\n")

    x_string = ''
    for i in range(len(xs) - 1):
        x_string += '{:+.6e} '.format(xs[i])
    x_string += '{:+.6e}\n'.format(xs[-1])
    f.write(x_string)

    Q_string = ''
    for i in range(len(Qs) - 1):
        Q_string += '{:+.6e} '.format(Qs[i])
    Q_string += '{:+.6e}\n'.format(Qs[-1])
    f.write(Q_string)

    f.write('-5 -4 -3 -2 -1 1 2 3 4 5 21\n')
    indices = list(range(-5, 0)) + list(range(1, 6)) + [0]

    for i in range(len(xs)):
        if ((i + 1) % 30) == 0:
            print(f'\t{i+1:03d}/{len(xs)}')
        x = xs[i]
        for j in range(len(Qs)):
            Q = Qs[j]
            if choice == 'pdf':
                values = x * theory.CxPDF(x, Q, Lmu, Lz)
            elif choice == 'ff':
                values = x * theory.CxFF (x, Q, Lmu, Lz)

            output = ''
            for index in indices[:-1]:
                output += '{:+.6e} '.format(values[index])
            output += '{:+.6e}\n'.format(values[-1])
            f.write(output)

    f.write("---\n")
    f.close()
    return

if __name__ == '__main__':
    ## x and Q values
    xs = np.logspace(-7, -1, num = 201, base = 10.0)
    xs = xs[:-1].copy()
    xs = np.append(xs, np.linspace(0.1, 1.0, 40))
    Qs = np.logspace(np.log10(1.3), 4, num = 40, base = 10.0)

    x_cut   = 0.999
    member  = 0
    error   = 1e-3

    Lmu     = 0.0
    Lz      = 0.0

    ## choose scheme and PDFs/FFs
    scheme  = 'MSbar'
    # scheme  = 'JCC'
    # scheme  = 'CSS'
    choice = 'pdf'
    # choice = 'ff'
    print(f'generating C times {choice.upper()}s under {scheme} scheme')
    print()

    t_0 = time.time()
    if choice == 'pdf':
        pdf_set = 'CT14nlo'
        order   = 0
        theory  = TMDPDF(order, x_cut, pdf_set, member, error, scheme)
        generate_output(choice, theory, xs, Qs, order, pdf_set, scheme, Lmu = Lmu, Lz = Lz)

        order   = 1
        theory  = TMDPDF(order, x_cut, pdf_set, member, error, scheme)
        generate_output(choice, theory, xs, Qs, order, pdf_set, scheme, Lmu = Lmu, Lz = Lz)

        generate_info(choice, pdf_set, scheme)
    elif choice == 'ff':
        ff_set = 'DSS14PI'
        order  = 0
        theory = TMDFF(order, x_cut, ff_set, member, error, scheme)
        generate_output(choice, theory, xs, Qs, order, ff_set, scheme, Lmu = Lmu, Lz = Lz)

        order  = 1
        theory = TMDFF(order, x_cut, ff_set, member, error, scheme)
        generate_output(choice, theory, xs, Qs, order, ff_set, scheme, Lmu = Lmu, Lz = Lz)

        generate_info(choice, ff_set, scheme)

    et  = time.time() - t_0
    et /= 60.0
    print(f'time used: {et:.2f} minutes')


