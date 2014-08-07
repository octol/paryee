# Simple script to collect the test results and produce a suitable overall
# value per gridsize (N). For now this entails taking the minimum value for
# each N.
#
# Usage:
#
#   python3 gather_data.py tests_scaling1.tsv tests_scaling2.tsv > tests_scaling.tsv
#
import sys
import numpy

first_file = True

for f in sys.argv[1:]:
    data = numpy.genfromtxt(f, delimiter=' ')

    if first_file:
        data_min = data
        first_file = False
    else:
        data_min = numpy.minimum(data_min, data)

numpy.savetxt('/dev/stdout', data_min, delimiter=' ')

