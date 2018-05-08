function col_index = pickcol(coordinates,W,L)
    %disp(max(coordinates(:,1))>W);disp(max(coordinates(:,2))>L);
    if max(coordinates(:,1))>W | max(coordinates(:,2))>L
        col_index=0;
    else
        col_index=(coordinates(:,1)-1)*L+coordinates(:,2);
        col_index=col_index';
    end
end

