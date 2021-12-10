% mle.m
load('TrainingSamplesDCT_subsets_8.mat')
load('Alpha.mat')
load('Prior_2.mat') % Change to Prior_2 to check 2nd strategy

[rg, cg] = size(D1_BG); % Change to D2, D3 or D4 for other datasets
[rc, cc] = size(D1_FG);

% MLE Prior
prior_BG = rg / (rc + rg);
prior_FG = rc / (rc + rg);

mu_BG = sum(D1_BG) / rg;
mu_FG = sum(D1_FG) / rc;

sigma_BG = cov(D1_BG);
sigma_FG = cov(D1_FG);

mle_alpha_err = [];
err = segment_cheetah(prior_FG, prior_BG, mu_FG, mu_BG, sigma_FG, sigma_BG, 1);
for i = 1:9
    mle_alpha_err = [mle_alpha_err err];
end

semilogx(alpha, mle_alpha_err)
