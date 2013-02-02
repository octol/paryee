from numpy import * 

data1 = genfromtxt('tests_swiftsure/tests_perf_8.tsv', delimiter=' ')
data2 = genfromtxt('tests_asuka/tests_perf_6.tsv', delimiter=' ')
data3 = genfromtxt('tests_europa/tests_perf_8.tsv', delimiter=' ')
data4 = genfromtxt('tests_europa_gcc/tests_perf_8.tsv', delimiter=' ')

assert(size(data1,0) == size(data2,0))
assert(size(data1,0) == size(data3,0))
assert(size(data1,0) == size(data4,0))
assert(size(data1,1) == size(data2,1))
assert(size(data1,1) == size(data3,1))
assert(size(data1,1) == size(data4,1))

N = data1[:,0]

# Use yee_mpi2 (position 4)
comparison = column_stack((N, data1[:,4], data2[:,4], data3[:,4], data4[:,4]))
savetxt('/dev/stdout', comparison, delimiter=' ')

