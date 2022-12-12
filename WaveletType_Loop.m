function [mse, peak_snr] = WaveletType_Loop(i, wname,iter_length,imdata,threshold_weight)


% brain = imread('brain.jpg');
% imdata = im2double(brain(:,:,1)); %black and white image, all layers same, conv to double [0 1]
% imdata = phantom('Modified Shepp-Logan', 256);

% size_dif = size(imdata,1)-size(imdata,2);
% imdata = [zeros(size(imdata,1), size_dif/2 +1) imdata zeros(size(imdata,1), size_dif/2 +1)]; %evenly pad zeros to square image, make divisble by 4
% imdata = [imdata;zeros(2,size(imdata,2))];
rows = size(imdata,1);
cols = size(imdata,2);


ft_weight = 1/sqrt(size(imdata,1)*size(imdata,2));
F_imdata = fftshift(fft2(imdata).*ft_weight);

%Use gaussian mask for testing
percent_sp = 0.2;

nn = 1;
mse_iters = zeros(1,nn);
psnr_iters = zeros(1,nn);



for big_iter = 1:nn %run 10 sets of random masks and collect average data

%     mask = make_gauss_mask(rows,percent_sp);
    mask = make_mask(rows, 2);
    %figure; imshow(abs(mask));title("Mask Image");

    F_imdata_sp = F_imdata.*mask;
    im_sp = ifft2(ifftshift(F_imdata_sp))./ft_weight;
    im_sp = im_sp./(abs(max(im_sp,[],'all')));
    im_og = im_sp;
    if i == 1
        figure;
        imshow(abs(imdata), []);
        title("Orignial Image");
        figure;
        imshow(abs(im_og), []);
        title("Sparse Image");
    end
    im_final = zeros(rows,cols);

    % iter_length = 55;
    mean_squared_error = ones(1,iter_length);
    snr = zeros(1,iter_length);
    % dwtmode('per');
    for n = 1:iter_length
        %wavelet transform
    %     im_sp_th = dualtreethresh(im_sp, 2, threshold_weight);
        if (i == 1)
            im_sp_th = dwtthresh(im_sp, 2, 'bior1.1', threshold_weight);
        elseif (i == 2)
            im_sp_th = dualtreethresh(im_sp, 2, threshold_weight);
        elseif (i == 3)
            % self1, self2, or dddtf1
            im_sp_th = dddtreethresh(im_sp, 2, 'self2', threshold_weight);
        end
        %go to k space and downsample
        F_sp_th = fftshift(fft2(im_sp_th).*ft_weight);
        F_sp_th_masked = F_sp_th.*(1-mask); %find k-space info to fill in gaps


        %fill in gaps, find iteration's "difference image"
        F_err = F_sp_th_masked + F_imdata_sp;
        delta_im = ifft2(ifftshift(F_err))./ft_weight;

        %add diff im to original in image space
        new_im = im_sp + delta_im;
        im_sp = new_im./abs(max(new_im, [],'all'));
        im_final = im_final + delta_im;
        im_final = im_final./abs(max(im_final,[],'all'));
        im_final(im_final < 10^(-10)) = 0;
        mean_squared_error(n) = immse(imdata, abs(im_final));
        snr(n) = psnr(abs(im_final),imdata);
    end

    minV = min(min(abs(imdata)));
    maxV = max(max(abs(imdata)));
%     figure;
% %     subplot(1,3,1); 
%     imshow(abs(imdata), [minV maxV]);
%     title("Orignial Image");
% %     subplot(1,3,2); 
%     imshow(abs(im_og), [minV maxV]);
%     title("Sparse Image");
%     subplot(1,3,3); 
    figure;
    imshow(abs(im_final), [minV maxV]);
    if (i == 1)
       title("DWT");
    elseif (i == 2)
       title("Dual-Tree CDWT");
    elseif (i == 3)
       title("Double Density Dual-Tree CDWT");
    end
    %
    %
    % figure; plot(mean_squared_error)
    % title("Mean Square Error vs Iterations")
    %
    % figure; plot(snr)
    % title("Peak SNR vs iterations")
    %
    % disp("Final MSE: "+mean_squared_error(end));
    % disp("Final PSNR: "+snr(end));
    mse_iters(big_iter) = mean_squared_error(end)
    psnr_iters(big_iter) = snr(end)-snr(1);

end
mse = mean(mse_iters);
peak_snr = mean(psnr_iters);


end




function im_sp_th = dddtreethresh(im_sp, level, fname, threshold_weight)
    wt = dddtree2('cplxdddt', abs(im_sp), level, fname); % self1, self2, or dddtf1
    dddLevel = numel(wt.cfs{1}(:)) / numel(wt.cfs{1}(:,:,1,1,1));
    for l = 1:dddLevel
        thresh = abs(threshold_weight*max(wt.cfs{1}(:,:,l), [], 'all'));
        wt.cfs{1}(:,:,l) = threshold2(wt.cfs{1}(:,:,l), thresh);
    end
    im_sp_th = idddtree2(wt);
end


function im_sp_th = dualtreethresh2(im_sp, level, threshold_weight)
    wt = dddtree2('cplxdt', abs(im_sp), level, 'dtf4'); % self1, self2, or dddtf1
    dddLevel = numel(wt.cfs{1}(:)) / numel(wt.cfs{1}(:,:,1,1,1));
    for l = 1:dddLevel
        thresh = abs(threshold_weight*max(wt.cfs{1}(:,:,l), [], 'all'));
        wt.cfs{1}(:,:,l) = threshold2(wt.cfs{1}(:,:,l), thresh);
    end
    im_sp_th = idddtree2(wt);
end

function im_sp_th = dualtreethresh(im_sp, level, threshold_weight)
    [a, d] = dualtree2(abs(im_sp), 'Level', level, 'FilterLength', 10);
    for l = 1:numel(d)
        thresh = abs(threshold_weight*max(d{l}(:), [], 'all'));
        d{l} = threshold2(d{l}, thresh);
    end
    im_sp_th = idualtree2(a, d);
end

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