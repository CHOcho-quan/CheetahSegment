% segment_cheetah.m
function output = segment_cheetah(Pc, Pg, mu_c, mu_g, sigma_c, sigma_g, k)
    cheetah_original = imread('cheetah.bmp');
    [r, c] = size(cheetah_original);
    cheetah_original = im2double(cheetah_original);
    cheetah = padarray(cheetah_original, [4, 4], 'replicate', 'both');
    
    res = [];
    for i = 5:r+4
        tmp = [];
        for j = 5:c+4
            area = cheetah([i - 3:i + 4], [j - 3:j + 4]);
            dct_res = dct2(area);
            feat = zigzag(dct_res);

            % Calculate D(f, c) & D(g, c)
            d_cheetah = (feat - mu_c) * inv(sigma_c) * (feat - mu_c).' + log(det(sigma_c)) - 2 * log(Pc);
            d_grass = (feat - mu_g) * inv(sigma_g) * (feat - mu_g).' + log(det(sigma_g)) - 2 * log(Pg);

            if d_cheetah >= d_grass
                tmp = [tmp 0];
            else
                tmp = [tmp 255];
            end
        end
        res = [res; tmp];
    end
    figure(k)
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
    output = error_rate;
end