function im_sp_th = dwtthresh(im_sp, n_dwt, fname, threshold_weight)
        cA_hist = struct;
        cH_hist = struct;
        cV_hist = struct;
        cD_hist = struct;

        cA = im_sp;
        for i = 1:n_dwt
            [cA,cH,cV,cD] = dwt2(cA, fname);
            cA_hist(i).data = cA;
            cH_hist(i).data = cH;
            cV_hist(i).data = cV;
            cD_hist(i).data = cD;
        end
        
        % only use threshold from first DWT
        thresh = abs(threshold_weight*max([cA_hist(1).data(:) ; cH_hist(1).data(:) ; cV_hist(1).data(:) ; cD_hist(1).data(:)], [], 'all'));
    
        for i = n_dwt:-1:1
            % thresholding
            % keep top thresh_coef values in wavelet transform, else set to 0
            cA_denoise = threshold2(cA_hist(i).data, thresh);
            cD_denoise = threshold2(cD_hist(i).data, thresh);
            cH_denoise = threshold2(cH_hist(i).data, thresh);
            cV_denoise = threshold2(cV_hist(i).data, thresh);
             
            % inv W trans
            % get denoised image
    %         im_sp_th = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise, LoR, HiR);
            im_sp_th = idwt2(cA_denoise, cH_denoise, cV_denoise, cD_denoise,  fname);
    
            if (i > 1) 
                cA_hist(i-1).data = im_sp_th;
            end
        end
end