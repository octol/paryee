% Read and suitably filter the computed results.
%function gather_scaling_data
clear all

samples = 4;

for ii = 1:samples
    A{ii} = dlmread(['tests_scaling' num2str(ii) '.tsv'],' ');
end

% Decompose data into more manageable form

% The min function does not work over cell arrays, hence we need to transform
% to arrays.
for ii = 1:samples
    yee_omp(:,ii) = A{ii}(:,2);
    yee_pthr(:,ii) = A{ii}(:,3);
    yee_mpi(:,ii) = A{ii}(:,4);
    yee_mpi2(:,ii) = A{ii}(:,5);
end

% Filter over different simulations
yee_min(:,1) = A{1}(:,1);
yee_min(:,2) = min(yee_omp,[],2);
yee_min(:,3) = min(yee_pthr,[],2);
yee_min(:,4) = min(yee_mpi,[],2);
yee_min(:,5) = min(yee_mpi2,[],2);

dlmwrite('/dev/stdout',yee_min,' ');

