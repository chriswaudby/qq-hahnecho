#!/bin/csh

bruk2pipe -verb -in ./ser                                              \
  -bad 0.0 -ext -aswap -AMX -decim 1312 -dspfvs 21 -grpdly 76          \
  -xN              4096  -yN               294  -zN               128  \
  -xT              2048  -yT               294  -zT               64   \
  -xMODE            DQD  -yMODE           Real  -zMODE        Complex  \
  -xSW        15243.902  -ySW          294.000  -zSW         1912.046  \
  -xOBS         950.450  -yOBS           1.000  -zOBS         238.995  \
  -xCAR           0.400  -yCAR           0.000  -zCAR          16.700  \
  -xLAB              1H  -yLAB             Tau  -zLAB             13C  \
  -ndim               3  -aq2D         Complex                         \
-out cube.fid

# run Julia script to apply receiver phase cycling
./proc.jl

# relaxation times
set tauList = (0.1 1.0 2.0 3.0 5.0 7.0 10.0 13.0 16.0 22.0 29.0 37.0 46.0 56.0)

nmrPipe -in cubeQQ.fid -fn TP                         \
| nmrPipe  -fn ZTP                                    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -83 -p1 0.00 -di -verb          \
| nmrPipe  -fn EXT -x1 1ppm -xn -0.7ppm -sw           \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb                                 \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -zf 2                               \
| nmrPipe  -fn FT -alt -neg                           \
| nmrPipe  -fn PS -p0 58.00 -p1 180.00 -di -verb      \
| pipe2xyz -out ft/test%03d.ft2 -y -ov

sortPlanes.com -in ./ft/test%03d.ft2 -out ./ft/test%03d.ft2 -tau $tauList -title
xyz2pipe -in ft/test%03d.ft2 >cubeQQ.ft

nmrPipe -in cubeDQ.fid -fn TP                         \
| nmrPipe  -fn ZTP                                    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -83 -p1 0.00 -di -verb          \
| nmrPipe  -fn EXT -x1 1ppm -xn -0.7ppm -sw           \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb                                 \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -zf 2                               \
| nmrPipe  -fn FT -alt -neg                           \
| nmrPipe  -fn PS -p0 123.00 -p1 180.00 -di -verb     \
| pipe2xyz -out ft/test%03d.ft2 -y -ov

sortPlanes.com -in ./ft/test%03d.ft2 -out ./ft/test%03d.ft2 -tau $tauList -title
xyz2pipe -in ft/test%03d.ft2 >cubeDQ.ft

