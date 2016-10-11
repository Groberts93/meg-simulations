%George Roberts 02/09/16

%Resolution script
sc_correlation_tests
datastruct = dir('sequencedata/*deg*.mat');
ldata = length(datastruct);
n_corr_big = zeros(ldata, 100);

for n_test = 1:100
    
    bf_corr_pts
    n_corr_big(:,n_test) = n_corr;
    delete sequencedata/* ;
    sc_correlation_tests
    disp(num2str(n_test));
end

save('T_35ft_9CM_REGGED', 'n_corr_big')