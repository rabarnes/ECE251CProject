close all
brain = imread('brain.jpg');
imdata = im2double(brain(:,:,1)); %black and white image, all layers same, conv to double [0 1]
size_dif = size(imdata,1)-size(imdata,2);
imdata = [zeros(size(imdata,1), size_dif/2) imdata zeros(size(imdata,1), size_dif/2)]; %evenly pad zeros to square image
rows = size(imdata,1);
cols = size(imdata,2);

ft_weight = 1/(norm(imdata,"fro"));
% ft_weight = 1/(max(imdata,[],'all'));
F_imdata = fftshift(fft2(imdata)*ft_weight);


%downsample, remove every mth row, keep middle third
m= 4;
mask = make_mask(rows,4);
F_imdata_sp = F_imdata.*mask;
im_sp = ifft2(ifftshift(F_imdata_sp))/ft_weight;


%wavelet transform
[C,S] = wavedec2(im_sp,2,'haar');
[H1,V1,D1] = detcoef2('all',C,S,1); %take wavelet transform w/ haar filter, h,v,d are detail matrices
A1 = appcoef2(C,S,'haar',1);

% thresh_pct = 0.8; %set bottom x% to zero, use 5% value as threshold
% C_sorted = sort(C,'ascend');
% thresh = abs(C_sorted(round(thresh_pct*size(C,2))));
thresh = abs(0.025*max(C));

for i = 1:size(C,2)
    if abs(C(i))<thresh
        C(i) = 0;
    else 
        C(i) = C(i) - thresh*exp(1j*angle(C(i))); %soft thresholding, subtract threshold val from all other pixels
    end

end

im_sp_th = waverec2(C,S,'haar'); %reconstructed image after thresholding

%go to k space and downsample
F_sp_th = fftshift(fft2(im_sp_th)*ft_weight);
F_sp_th_masked = F_sp_th.*mask;


%find err in k space and ifft to get difference image
F_err = F_imdata_sp - F_sp_th_masked; %y - Ax error between 2D FFT recon and wavelet recon in k space
delta_im = ifft2(ifftshift(F_err))/ft_weight;

%add diff im to original in image space
new_im = im_sp - delta_im;
new_F_imdata = fftshift(fft2(new_im)*ft_weight);

%show images

figure();
imshow(imdata);
title("Original Brain MRI image, full k-space")

figure();
imshow(abs(F_imdata))
title("Original Brain MRI image in k-space")

figure();
imshow(abs(F_imdata_sp))
title("Sampled k-space")

figure();
imshow(abs(im_sp))
title("2D IFFT reconstructed Sampled Image")

figure();
subplot(2,2,1)
imshow(abs(A1))
title("Coarse Approximation Wavelet Coefficients")
subplot(2,2,2)
imshow(abs(H1))
title("Horizontal Detail Wavelet Coefficients")
subplot(2,2,3)
imshow(abs(V1))
title("Vertical Detail Wavelet Coefficients")
subplot(2,2,4)
imshow(abs(D1))
title("Diagonal Detail Wavelet Coefficients")

figure();
imshow(abs(im_sp_th))
title("reconstructed image after thresholding wavelet decompositions")

figure();
imshow(abs(F_sp_th))
title("k space of thresholding");


figure();
imshow(abs(delta_im))
title("Difference Image")

figure();
imshow(abs(new_im))
title("Final Reconstructed Image after 1 loop")


