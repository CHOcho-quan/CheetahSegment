% map.m
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

map_alpha_err = [];
for i = 1:9
    a = alpha(i);
    Sigma0 = diag(a * W0);

    % Posterior Mean \prop P(D|mu)P(mu)
    mu1_BG = Sigma0 * inv(Sigma0 + sigma_BG / rg) * mu_BG.' + sigma_BG * inv(Sigma0 + sigma_BG / rg) * mu0_BG.' / rg;
    mu1_FG = Sigma0 * inv(Sigma0 + sigma_FG / rc) * mu_FG.' + sigma_FG * inv(Sigma0 + sigma_FG / rc) * mu0_FG.' / rc;

    sigma1_BG = Sigma0 * inv(Sigma0 + sigma_BG / rg) * sigma_BG / rg;
    sigma1_FG = Sigma0 * inv(Sigma0 + sigma_FG / rc) * sigma_FG / rc;
    
    err = segment_cheetah(prior_FG, prior_BG, mu1_FG.', mu1_BG.', sigma_FG, sigma_BG, i);
    map_alpha_err = [map_alpha_err err];
end

semilogx(alpha, map_alpha_err)
