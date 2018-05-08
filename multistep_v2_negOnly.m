function [delta_alpha, locs, widths, lengths, BoundingBox, alpha_calc] = multistep_v2_negOnly(E_vec, C, WW, LL, threshold, sigma, stepsize, iterations)
% The Localized ridge regression algorithm, which can be extended to any 
% number of iterations
    mask = ones(WW, LL);
    locs = zeros(2, iterations);
    widths = zeros(1, iterations);
    lengths = zeros(1, iterations);
    
    mask_y = sum(mask, 1); mask_x = sum(mask, 2);
    y1 = find(mask_y ~= 0, 1, 'first'); y2 = find(mask_y ~= 0, 1, 'last');
    x1 = find(mask_x ~= 0, 1, 'first'); x2 = find(mask_x ~= 0, 1, 'last');
    
    for i = 1 : iterations
        coordinates = zeros(length(find(mask == 1)), 2);
        k = 1;
        for x = 1 : WW
            for y = 1 : LL
                if mask(x, y) == 1
                    coordinates(k, :)=[x, y];
                    k = k + 1;
                end
            end
        end
        %disp(coordinates);
        col_index = pickcol(coordinates, WW, LL);
        C_new = C(:, col_index);
        delta_alpha_local = (C_new' * C_new + sigma(i) * stepsize ^ 2 * eye(size(C_new, 2))) \ (C_new' * E_vec);
        delta_alpha = zeros(WW * LL, 1);
        delta_alpha(col_index) = delta_alpha_local;
        [~, ~, alpha_calc] = plot_alpha(delta_alpha, WW, LL, stepsize);
        mask = double(alpha_calc <= threshold * min(min(alpha_calc)));
        if i >= 1 % Preserve the most intense connected component only
            CC = bwconncomp(mask, 4);
            if length(CC.PixelIdxList) > 1 % If more than one connected component
                sum_abs_delta_alpha = zeros(1, length(CC.PixelIdxList));
                for k = 1:length(CC.PixelIdxList)
                    sum_abs_delta_alpha(k) = sum(abs(alpha_calc(CC.PixelIdxList{k})));
                end
                [~, index] = max(sum_abs_delta_alpha);
                for k = 1:length(CC.PixelIdxList)
                    if k ~= index
                        mask(CC.PixelIdxList{k}) = 0;
                    end
                end
            end
        end
        mask_y = sum(mask, 1); mask_x = sum(mask, 2);
        y1 = min(find(mask_y ~= 0)); y2 = max(find(mask_y ~= 0));
        x1 = min(find(mask_x ~= 0)); x2 = max(find(mask_x ~= 0));
        %y1=find(mask_y~=0,1,'first');y2=find(mask_y~=0,1,'last');
        %x1=find(mask_x~=0,1,'first');x2=find(mask_x~=0,1,'last');
        loc_x = sum(stepsize * ((1 : WW) - 0.5) .* (sum(abs(alpha_calc), 2))') /...
                sum(sum(abs(alpha_calc), 2));
        loc_y = sum(stepsize * ((1 : LL) - 0.5) .* sum(abs(alpha_calc), 1)) /...
                sum(sum(abs(alpha_calc), 1));
        locs(1, i) = loc_x; locs(2, i) = loc_y;
        %disp(x1);disp(x2);
        %save('tmp.mat','delta_alpha','alpha_calc');
        %widths(i)=x2(1)-x1(1)+1;
        %lengths(i)=y2(1)-y1(1)+1;
    end    
    BoundingBox = [y1 - 1, x1 - 1, y2 - y1 + 1, x2 - x1 + 1];
end

