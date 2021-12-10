load('TrainingSamplesDCT_8_new.mat')

[rc, cc] = size(TrainsampleDCT_BG);
[rg, cg] = size(TrainsampleDCT_FG);

prior_hist = [];
for i = 1:rc
    prior_hist = [prior_hist 0];
end
for i = 1:rg
    prior_hist = [prior_hist 1];
end

figure(1)
histogram(prior_hist, 2, 'Normalization', 'probability');
