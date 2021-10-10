#!/bin/csh


bruk2pipe -verb -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1312 -dspfvs 21 -grpdly 76  \
  -xN              4096  -yN               294  -zN               128  \
  -xT              2048  -yT               294  -zT               64  \
  -xMODE            DQD  -yMODE           Real  -zMODE        Complex  \
  -xSW        15243.902  -ySW          294.000  -zSW         1912.046  \
  -xOBS         950.450  -yOBS           1.000  -zOBS         238.995  \
  -xCAR           0.400  -yCAR           0.000  -zCAR          16.700  \
  -xLAB              1H  -yLAB             Tau  -zLAB             13C  \
  -ndim               3  -aq2D         Complex                         \
-out cube.fid
