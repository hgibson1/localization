function [x] = DrawRectangle2(alpha_2d, WW, LL, stepsize)
    x = [];
    for i = 1 : WW
        for j = 1 : LL
            if abs(alpha_2d(i, j)) >=0% 1e-6
                x = [x; [i - 1, j - 1, 1, 1, alpha_2d(i,j)]];
            end
        end
    end
    x(:, 1:4) = x(:, 1:4) * stepsize;
    x(x(:, 5) > 1, 5) = 1;
    x(x(:, 5) < -1, 5) = -1;
end

