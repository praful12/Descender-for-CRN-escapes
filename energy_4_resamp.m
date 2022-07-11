%function for energy constraint in optimal resampling

function tot_en = energy_4_resamp(c,b,f_c,ft_delta_x_arr)

seg_4_resamp = size(ft_delta_x_arr,1);
tot_en = 0;
for i = 1:seg_4_resamp
    ftdx = ft_delta_x_arr{i,1};
    sum_ftdx = sum(ftdx,2);
    %tot_en = tot_en + c(i)^2*sum_ftdx(ceil(b(i)/c(i)*f_c)+1);
    %tot_en = tot_en + 1/c(i)^2*sum_ftdx(ceil(b(i)/c(i)*f_c));
    sum_idx = min(ceil(c(i)/b(i)*f_c),size(sum_ftdx,1));
    tot_en = tot_en + sum_ftdx(sum_idx);
end


end