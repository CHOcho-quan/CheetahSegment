load('TrainingSamplesDCT_8_new.mat')

[rc, cc] = size(TrainsampleDCT_BG);
[rg, cg] = size(TrainsampleDCT_FG);

cheetah_mus_ori = [];
cheetah_sigmas_ori = [];
grass_mus_ori = [];
grass_sigmas_ori = [];
b_distances = [];
for selected = 1:64
    % Dealing with cheetah
    cheetah_mu = 0.0;
    cheetah_sigma = 0.0;

    for i = 1:rg
        feature = TrainsampleDCT_FG(i, :);
        cheetah_mu = cheetah_mu + feature(selected);
    end
    cheetah_mu = cheetah_mu / rg;
    cheetah_mus_ori = [cheetah_mus_ori cheetah_mu];

    for i = 1:rg
        feature = TrainsampleDCT_FG(i, :);
        cheetah_sigma = cheetah_sigma + (feature(selected) - cheetah_mu)^2;
    end
    cheetah_sigma = cheetah_sigma / rg;
    cheetah_sigmas_ori = [cheetah_sigmas_ori cheetah_sigma];

    % Dealing with grass
    grass_mu = 0.0;
    grass_sigma = 0.0;
    for i = 1:rc
        feature = TrainsampleDCT_BG(i, :);
        grass_mu = grass_mu + feature(selected);
    end
    grass_mu = grass_mu / rc;
    grass_mus_ori = [grass_mus_ori grass_mu];

    for i = 1:rc
        feature = TrainsampleDCT_BG(i, :);
        grass_sigma = grass_sigma + (feature(selected) - grass_mu)^2;
    end
    grass_sigma = grass_sigma / rc;
    grass_sigmas_ori = [grass_sigmas_ori grass_sigma];
    
    % Calculate Bhattahcharyya Distance to Identify Good Feature
    b_dis = (grass_mu - cheetah_mu)^2 / (grass_sigma + cheetah_sigma);
    b_distances = [b_distances b_dis];
end

[BD, ind] = sort(b_distances, 'descend');
cheetah_mus = cheetah_mus_ori(ind);
cheetah_sigmas = cheetah_sigmas_ori(ind);
grass_mus = grass_mus_ori(ind);
grass_sigmas = grass_sigmas_ori(ind);

% plot best 8
for i = 1:8
    cur_cheetah_mu = cheetah_mus(i);
    cur_cheetah_sigma = cheetah_sigmas(i);
    cur_grass_mu = grass_mus(i);
    cur_grass_sigma = grass_sigmas(i);
    
    % Plot Graph
    x = [min(-0.1, min(cur_cheetah_mu - 3 * sqrt(cur_cheetah_sigma), cur_grass_mu - 3 * sqrt(cur_grass_sigma))):0.001:max(0.1, max(cur_cheetah_mu + 3 * sqrt(cur_cheetah_sigma), cur_grass_mu + 3 * sqrt(cur_grass_sigma)))];
    y = mynormpdf(x, cur_cheetah_mu, sqrt(cur_cheetah_sigma));
    y2 = mynormpdf(x, cur_grass_mu, sqrt(cur_grass_sigma));
    subplot(2,4,i);
    figure(1)
    plot(x, y, x, y2);
    legend({'cheetah','grass'},'Location','southwest');
    title("Best-8");
end
hold on

% plot worst 8
for i = 57:64
    cur_cheetah_mu = cheetah_mus(i);
    cur_cheetah_sigma = cheetah_sigmas(i);
    cur_grass_mu = grass_mus(i);
    cur_grass_sigma = grass_sigmas(i);
    
    % Plot Graph
    x = [min(-0.1, min(cur_cheetah_mu - 3 * sqrt(cur_cheetah_sigma), cur_grass_mu - 3 * sqrt(cur_grass_sigma))):0.001:max(0.1, max(cur_cheetah_mu + 3 * sqrt(cur_cheetah_sigma), cur_grass_mu + 3 * sqrt(cur_grass_sigma)))];
    y = mynormpdf(x, cur_cheetah_mu, sqrt(cur_cheetah_sigma));
    y2 = mynormpdf(x, cur_grass_mu, sqrt(cur_grass_sigma));
    figure(2)
    subplot(2,4,i-56);
    plot(x, y, x, y2);
    legend({'cheetah','grass'},'Location','southwest');
    title("Worst-8");
end
hold off
