function [x, y, alpha] = plot_alpha(alpha, W, L, stepsize)
    alpha = vec2mat(alpha, L);
    if max(max(abs(alpha))) ~= 0
        alpha = alpha / max(max(abs(alpha)));
    end
    %alpha=alpha/max(max(alpha));
    [x, y] = meshgrid(1:L, 1:W);
    x = x * stepsize;
    y = y * stepsize;
    %figure;
    %contourf(x,y,alpha);colorbar;
end

