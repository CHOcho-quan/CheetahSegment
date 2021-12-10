function [mu, sigma, pi] = EM_learn(dataset, C)
    [N, d] = size(dataset);
    
    % Randomized Initialization
    mu0 = zeros(C, d);
    sigma0 = zeros(d, d, C);
    index = randi(N, 1, C);
    for i = 1:C
        sigma0(:, :, i) = diag(ones(1, d) .* rand(1, d) + 1e-5);
        mu0(i, :) = dataset(index(i), :);
    end
    
    % EM Procedure
    pi0 = ones(1, C) / C;
    mu = mu0;
    sigma = sigma0;
    pi = pi0;
    for i = 1:5000
        % E-step solve H matrix
        H = zeros(N, C);
        for c = 1:C
            mu_c = mu0(c, :);
            sigma_c = sigma0(:, :, c);
            H(:, c) = mvnpdf(dataset, mu_c, sigma_c);
            H(:, c) = H(:, c) * pi0(c);
        end
        H = H ./ sum(H, 2);
        
        % M-step getting new parameters
        pi = sum(H) / N;
        for c = 1:C
            mu(c, :) = sum(H(:, c) .* dataset) / (sum(H(:, c)));
            x_c = dataset - mu0(c, :);
            s = sum((x_c .* x_c) .* H(:, c) / (sum(H(:, c)))) + 1e-5;
            sigma(:, :, c) = diag(s);
        end
        
        % Early Stopping
        if sum(abs(mu - mu0 ./ mu0), 'all') < 1e-2
            break;
        end

        mu0 = mu;
        sigma0 = sigma;
        pi0 = pi;
    end
end