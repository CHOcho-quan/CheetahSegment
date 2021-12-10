% EM.m
load('/Users/quan/Documents/MATLAB/SL/5/TrainingSamplesDCT_8_new.mat')

[rg, cg] = size(TrainsampleDCT_BG);
[rc, cc] = size(TrainsampleDCT_FG);

% Reading Image
img = imread('cheetah.bmp');
img = im2double(img);
[rows,cols] = size(img);
img_mask = imread('cheetah_mask.bmp');
img_mask = img_mask / 255;

% Using Matrix Representation
zigzags = zeros(rows,cols,64);
for row = 1:rows-7
    for col = 1:cols-7
        DCT = (dct2(img(row:row+7,col:col+7)));
        zigzag_matrix = zigzag(DCT);
        zigzag_matrix = zigzag_matrix';
        zigzags(row,col,:) = zigzag_matrix;
    end
end

BG_m = containers.Map();
BG_s = containers.Map();
BG_p = containers.Map();
FG_m = containers.Map();
FG_s = containers.Map();
FG_p = containers.Map();
C = 8;

% EM Training Phase
tic;
for i = 1:5
    [mu_grass, sigma_grass, pi_grass] = EM_learn(TrainsampleDCT_BG, C);
    [mu_cheetah, sigma_cheetah, pi_cheetah] = EM_learn(TrainsampleDCT_FG, C);
    BG_m(num2str(i)) = mu_grass;
    BG_s(num2str(i)) = sigma_grass;
    BG_p(num2str(i)) = pi_grass;
    FG_m(num2str(i)) = mu_cheetah;
    FG_s(num2str(i)) = sigma_cheetah;
    FG_p(num2str(i)) = pi_cheetah;
end
toc

% Compute Error
dim = [1, 2, 4, 8, 16, 24, 32, 40, 48, 56, 64];
errors = zeros(5, 5, 11);
for i = 1:5
    for j = 1:5
        for di = 1:size(dim, 2)
            res = zeros(rows, cols);
            d = dim(di);
            mu_BG = BG_m(num2str(i));
            mu_BG = mu_BG(:, 1:d);
            mu_FG = FG_m(num2str(j));
            mu_FG = mu_FG(:, 1:d);
            sigma_BG = BG_s(num2str(i));
            sigma_BG = sigma_BG(1:d, 1:d, :);
            sigma_FG = FG_s(num2str(j));
            sigma_FG = sigma_FG(1:d, 1:d, :);
            pi_BG = BG_p(num2str(i));
            pi_FG = FG_p(num2str(j));

            % Prior
            prior_BG = rg / (rc + rg);
            prior_FG = rc / (rc + rg);
            
            inv_BG = zeros(d, d, C);
            inv_FG = zeros(d, d, C);
            alpha_BG = zeros(1, C);
            alpha_FG = zeros(1, C);
            
            % Calculate Inverse Ahead
            for c = 1:C
                inv_BG(:, :, c) = pinv(sigma_BG(:, :, c));
                inv_FG(:, :, c) = pinv(sigma_FG(:, :, c));
                alpha_BG(1, c) = sum(log(diag(sigma_BG(:, :, c))));
                alpha_FG(1, c) = sum(log(diag(sigma_FG(:, :, c))));
            end
            
            % Calculate Error
            for row = 1:rows-7
                for col = 1:cols-7
                    zigzag_matrix = zigzags(row, col, :);
                    zigzag_matrix = zigzag_matrix(:);
                    zigzag_matrix = zigzag_matrix(1:d);
                    pBG = 0;
                    pFG = 0;
                    for c = 1:C
                        cur_mu_FG = mu_FG(c, :)';
                        cur_mu_BG = mu_BG(c, :)';
                        X_FG = zigzag_matrix - cur_mu_FG;
                        X_BG = zigzag_matrix - cur_mu_BG;
                        pBG = pBG + exp(log(pi_BG(c)) - 0.5 * X_BG' * inv_BG(:, :, c) * X_BG - 0.5 * alpha_BG(c));
                        pFG = pFG + exp(log(pi_FG(c)) - 0.5 * X_FG' * inv_FG(:, :, c) * X_FG - 0.5 * alpha_FG(c));
                    end
                    pFG = log(pFG) + log(prior_FG);
                    pBG = log(pBG) + log(prior_BG);
                    if pFG >= pBG
                        res(row, col) = 1;
                    end
                end
            end
            err = sum(sum(res ~= img_mask)) / (rows*cols)
            errors(i, j, di) = err;
        end
    end
end