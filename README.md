# TMDs-NLO
Code that produces TMDs up to NLO matching order and NLL evolution order

Before running, we must make 2 python wrappers for 2 softwares, namely apfel (for strong coupling) and kernel_q (for perturbative evolution). Here are the steps:

1. In TMDs-NLO/lib/Evolution run 'make'
2. To clean and remake the above, simply remove the resulting .so compiled library
3. In TMDs-NLO/lib/apfel run 'cmake .' and then 'make install'
4. To clean and remake the above, use: 'rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake Makefile CTestTestfile.cmake _deps lib'

Other tips:
1. Do the above in your working conda environment so that the compiled wrappers match your working python version
2. Example LHAPDF grids are included for collinear fragmentation functions, though local grids can be included (see example.ipynb)
