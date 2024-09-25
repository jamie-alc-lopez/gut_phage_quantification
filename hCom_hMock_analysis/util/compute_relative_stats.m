function [cent,d_std] = compute_relative_stats(V,x)

%This function computes (1) the percentile of the value x with respect to the set
%of data V and (2) the number of standard deviations from the mean of V the
%value x is

%Compute centile
n_less = sum(V < x);
n_eq = sum(V == x);
cent = 100*(n_less + 0.5*n_eq)/length(V);

%Compute number of standard deviations
d_std = abs(mean(V) - mean(x))/std(V);

end