function [delta_alphas, shadows, intersection, delta_alpha2, loc_est] = loc_algo_3D_v2(E, C, WW, LL, threshold, sigmas, stepsize, LargestCConly)
% Localize a 3D object by intersecting shadows from all light sources
    Ns = size(E, 1);
    Nf = size(E, 2);
    delta_alphas = zeros(WW * LL, Nf);
    shadows = zeros(WW * LL, Nf);
    loc_est = [-100, -100];
    norms = zeros(1, Nf);
    for j = 1:Nf
        norms(j) = norm(E(:, j));
    end
    max_norm = max(norms);
    for j = 1:Nf
        E_sub = E(:, j);
        C_sub = C((j - 1) * Ns + 1 : j * Ns, :);
        if norms(j) >= 0.01 * max_norm
            delta_alphas(:, j) = (C_sub' * C_sub + sigmas(1) * stepsize^2 * eye(size(C_sub, 2))) \ (C_sub' * E_sub);
            d_alpha_max = max(abs(delta_alphas(:, j)));
            shadows(:, j) = double(delta_alphas(:, j) <= -d_alpha_max * threshold);
        else % if norm of the column of E is too small, define shadow as the whole floor
            shadows(:, j) = ones(WW * LL, 1);
        end
    end
    if LargestCConly == 1 % Keep only the connected component with the maximum sum in the estimated shadow
        for j = 1:Nf
            shadow2D = vec2mat(shadows(:, j), LL);
            delta_alpha_2D = abs(vec2mat(delta_alphas(:, j), LL));
            CC = bwconncomp(shadow2D, 4);
            if length(CC.PixelIdxList) > 1 % If more than one connected component
                sum_abs_delta_alpha = zeros(1, length(CC.PixelIdxList));
                for k = 1:length(CC.PixelIdxList)
                    sum_abs_delta_alpha(k) = sum(delta_alpha_2D(CC.PixelIdxList{k}));
                end
                [~, index] = max(sum_abs_delta_alpha);
                for k = 1:length(CC.PixelIdxList)
                    if k ~= index
                        shadow2D(CC.PixelIdxList{k}) = 0;
                    end
                end
            end
            shadow2D = shadow2D';
            shadows(:, j) = shadow2D(:);
        end
    end
    intersection = double(sum(shadows, 2) >= Nf);
    delta_alpha2 = zeros(WW * LL, 1);
    col_index = find(intersection == 1);
    if ~isempty(col_index) && max(max(abs(E))) ~= 0
        C_sub = C(:, col_index);
        delta_alpha_local = (C_sub' * C_sub + sigmas(2) * stepsize ^ 2 * eye(size(C_sub, 2))) \ (C_sub' * E(:));
        delta_alpha_local = delta_alpha_local .* double(delta_alpha_local <= 0);
        delta_alpha2(col_index) = delta_alpha_local;
        delta_alpha2_2D = vec2mat(delta_alpha2, LL);
        %intersection2D = vec2mat(intersection, LL);
        loc_est(1) = sum(stepsize * ((1 : WW) - 0.5) .* (sum(abs(delta_alpha2_2D), 2))') /...
                     sum(sum(abs(delta_alpha2_2D), 2));
        loc_est(2) = sum(stepsize * ((1 : LL) - 0.5) .* sum(abs(delta_alpha2_2D), 1)) /...
                     sum(sum(abs(delta_alpha2_2D), 1));
    end
end
