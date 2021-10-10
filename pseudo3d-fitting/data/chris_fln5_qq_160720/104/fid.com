#!/bin/csh

set tauList = (1.1 50.0 100.0 150.0)

bruk2pipe -verb -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1248 -dspfvs 21 -grpdly 76  \
  -xN              3072  -yN                 4  -zN               340  \
  -xT              1536  -yT                 4  -zT               170  \
  -xMODE            DQD  -yMODE           Real  -zMODE    States-TPPI  \
  -xSW        16025.641  -ySW            4.000  -zSW         3017.502  \
  -xOBS         800.200  -yOBS           1.000  -zOBS         201.214  \
  -xCAR           0.400  -yCAR           0.000  -zCAR          16.700  \
  -xLAB              1H  -yLAB             Tau  -zLAB             13C  \
  -ndim               3  -aq2D         Complex                         \
| nmrPipe -fn TP \
| nmrPipe -fn ZTP \
| nmrPipe -fn TP \
| pipe2xyz -x -out ./fid/test%03d.fid -ov

sortPlanes.com -in ./fid/test%03d.fid -out ./fid/test%03d.fid -tau $tauList -title

