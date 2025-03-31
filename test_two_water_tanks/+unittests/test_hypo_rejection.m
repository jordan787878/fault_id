function test_hypo_rejection()

% This unit test shows that the false positive (rejecting a null hypothesis
% even though it is true) is controlled by alpha

N_trials = 10000;
df = 4;
alpha = 0.05;
N = 1;

count = 0;
for j = 1:N_trials
    chi_crit = chi2inv(1-alpha, N*df)/N;
    X = mvnrnd(zeros(df,1), eye(df), N);
    X_aug = reshape(X', [], 1);
    chi = (X_aug' * X_aug)/N;
    if(chi < chi_crit)
        count = count + 1;
    end
end
disp(100*count/N_trials)


end