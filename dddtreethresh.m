function im_sp_th = dddtreethresh(im_sp, level, fname, threshold_weight)
    wt = dddtree2('cplxdddt', abs(im_sp), level, fname); % self1, self2, or dddtf1
    dddLevel = numel(wt.cfs{1}(:)) / numel(wt.cfs{1}(:,:,1,1,1));
    for l = 1:dddLevel
        thresh = abs(threshold_weight*max(wt.cfs{1}(:,:,l), [], 'all'));
        wt.cfs{1}(:,:,l) = threshold2(wt.cfs{1}(:,:,l), thresh);
    end
    im_sp_th = idddtree2(wt);
end