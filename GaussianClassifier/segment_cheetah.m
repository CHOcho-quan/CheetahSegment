% Preprocess
warning('off','all')
cheetah_original = imread('cheetah.bmp');
[r, c] = size(cheetah_original);
cheetah_original = im2double(cheetah_original);
cheetah = padarray(cheetah_original, [4, 4], 'replicate', 'both');

Pcheetah = 250 / (1053 + 250);
Pgrass = 1053 / (1053 + 250);

res = [];
best8 = ind(1:8);
cheetah_mu8 = cheetah_mus(1:8);
Sigma_cheetah = diag(zeros(8));

for i = 1:rg
    feature = TrainsampleDCT_FG(i, :);
    feature = feature(best8);
    Sigma_cheetah = Sigma_cheetah + (feature - cheetah_mu8).'*(feature - cheetah_mu8);
end
Sigma_cheetah = Sigma_cheetah / rg;

grass_mu8 = grass_mus(1:8);
Sigma_grass = diag(zeros(8));
for i = 1:rc
    feature = TrainsampleDCT_BG(i, :);
    feature = feature(best8);
    Sigma_grass = Sigma_grass + (feature - grass_mu8).'*(feature - grass_mu8);
end
Sigma_grass = Sigma_grass / rc;

for i = 5:r+4
    tmp = [];
    for j = 5:c+4
        area = cheetah([i - 3:i + 4], [j - 3:j + 4]);
        dct_res = dct2(area);
        feat = zigzag(dct_res);
        
        % Calculate by best feature
        best_feat = feat(best8);
        d_cheetah = (best_feat - cheetah_mu8) * inv(Sigma_cheetah) * (best_feat - cheetah_mu8).' + log(det(Sigma_cheetah)) - 2 * log(Pcheetah);
        d_grass = (best_feat - grass_mu8) * inv(Sigma_grass) * (best_feat - grass_mu8).' + log(det(Sigma_grass)) - 2 * log(Pgrass);
        
        if d_cheetah >= d_grass
            tmp = [tmp 0];
        else
            tmp = [tmp 255];
        end
    end
    res = [res; tmp];
end
img = imagesc(res);

mask = imread("cheetah_mask.bmp");
error = 0; total = r * c;
for i = 1:r
    for j = 1:c
        mask(i, j);
        if mask(i, j) ~= res(i, j)
            error = error + 1;
        end
    end
end

error_rate = error / total
