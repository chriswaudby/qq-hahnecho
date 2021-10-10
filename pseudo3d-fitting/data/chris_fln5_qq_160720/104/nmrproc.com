#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

xyz2pipe -in fid/test%03d.fid -x -verb \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 93.00 -p1 0.00 -di -verb         \
| nmrPipe -fn EXT -x1 1ppm -xn -0.7ppm -sw                                    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb \
| nmrPipe  -fn EM -lb 10 -c 0.5 \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 0.00 -p1 0.00 -di -verb         \
| pipe2xyz -out ft/test%03d.ft2 -y -ov
