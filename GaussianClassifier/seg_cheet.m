% Preprocess
warning('off','all')
cheetah_original = imread('cheetah.bmp');
[r, c] = size(cheetah_original);
cheetah_original = im2double(cheetah_original);
cheetah = padarray(cheetah_original, [4, 4], 'replicate', 'both');

Pcheetah = 250 / (1053 + 250);
Pgrass = 1053 / (1053 + 250);

res = [];
cmu = cheetah_mus_ori;
Sigma_cheetah = diag(zeros(64));

for i = 1:rg
    feature = TrainsampleDCT_FG(i, :);
    Sigma_cheetah = Sigma_cheetah + (feature - cmu).' * (feature - cmu);
end
Sigma_cheetah = Sigma_cheetah / rg;

gmu = grass_mus_ori;
Sigma_grass = diag(zeros(64));
for i = 1:rc
    feature = TrainsampleDCT_BG(i, :);
    Sigma_grass = Sigma_grass + (feature - gmu).' * (feature - gmu);
end
Sigma_grass = Sigma_grass / rc;

for i = 5:r+4
    tmp = [];
    for j = 5:c+4
        area = cheetah([i - 3:i + 4], [j - 3:j + 4]);
        dct_res = dct2(area);
        feat = zigzag(dct_res);
        
        % Calculate D(f, c) & D(g, c)
        d_cheetah = (feat - cmu) * inv(Sigma_cheetah) * (feat - cmu).' + log(det(Sigma_cheetah)) - 2 * log(Pcheetah);
        d_grass = (feat - gmu) * inv(Sigma_grass) * (feat - gmu).' + log(det(Sigma_grass)) - 2 * log(Pgrass);
        
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
