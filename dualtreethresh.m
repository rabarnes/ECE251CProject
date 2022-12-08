function im_sp_th = dualtreethresh(im_sp, level, threshold_weight)
    [a, d] = dualtree2(abs(im_sp), 'Level', level);
    for l = 1:numel(d)
        thresh = abs(threshold_weight*max(d{l}(:), [], 'all'));
        d{l} = threshold2(d{l}, thresh);
    end
    im_sp_th = idualtree2(a, d);
end