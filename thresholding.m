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


%downsample, remove every mth row
m= 5;
F_imdata_sp = F_imdata;
F_imdata_sp(1:m:round(2*end/5),:) = zeros(length(1:m:round(2*rows/5)),cols);
F_imdata_sp(round(3*end/5):m:end,:) = zeros(length(round(3*rows/5):m:rows),cols);
im_sp = ifft2(ifftshift(F_imdata_sp))/ft_weight;


%wavelet transform
[C,S] = wavedec2(im_sp,2,'haar');
[H1,V1,D1] = detcoef2('all',C,S,1); %take wavelet transform w/ haar filter, h,v,d are detail matrices
A1 = appcoef2(C,S,'haar',1);

thresh = 0.05*max(C);
for i = 1:size(C,2)
    if abs(C(i))<thresh
        C(i) = 0;
    else 
        C(i) = C(i) - thresh*exp(1j*angle(C(i))); %soft thresholding, subtract threshold val from all other pixels
    end

end

im_sp_th = waverec2(C,S,'haar');

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


