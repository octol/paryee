% Gather data from the different tests runs and create data suitable for
% comparison plots
clear all

A{1} = dlmread('tests_swiftsure/tests_perf_4.tsv',' ');
A{2} = dlmread('tests_asuka/tests_perf_4.tsv',' ');
A{3} = dlmread('tests_europa/tests_perf_4.tsv',' ');

assert(all(size(A{1}) == size(A{2})))
assert(all(size(A{1}) == size(A{3})))

B(:,1) = A{1}(:,1); % N: number of cells

% Use yee_mpi2
B(:,2) = A{1}(:,5); 
B(:,3) = A{2}(:,5); 
B(:,4) = A{3}(:,5); 

dlmwrite('/dev/stdout',B,' ');

