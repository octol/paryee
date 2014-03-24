# Simple script to collect the test results and produce a suitable overall
# value per gridsize (N). For now this entails taking the minimum value for
# each N.
#
# Usage:
#
#   python3 gather_data.py tests_scaling > tests_scaling.tsv
#   python3 gather_data.py tests_perf 8 > tests_perf_8.tsv
#
import sys
from numpy import *

samples = 4

assert len(sys.argv) >= 2, 'Please specify: basename!'
basename = sys.argv[1]
if len(sys.argv) >= 3:
    nodes = sys.argv[2]

for i in range(1,samples+1):
    filename = basename + str(i)
    if len(sys.argv) >= 3:
        filename += '_' + str(nodes)
    filename += '.tsv'

    data = genfromtxt(filename, delimiter=' ')
    if i == 1:
        data_min = data
    else:
        data_min = minimum(data_min,data)

savetxt('/dev/stdout', data_min, delimiter=' ') 
