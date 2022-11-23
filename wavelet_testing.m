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



