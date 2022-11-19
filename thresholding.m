brain = imread('brain.jpg');
imdata = im2double(brain(:,:,1)); %black and white image, all layers same, conv to double [0 1]
size_dif = size(imdata,1)-size(imdata,2);
imdata = [zeros(size(imdata,1), size_dif/2) imdata zeros(size(imdata,1), size_dif/2)]; %evenly pad zeros to square image
rows = size(imdata,1);
cols = size(imdata,2);

ft_weight = 1/(norm(imdata,"fro"));
F_imdata = fft2(imdata)*ft_weight;


%downsample
m= 20;
F_imdata_sp = F_imdata;
F_imdata_sp(1:m:end,:) = zeros(length(1:m:rows),cols);

%wavelet transform


new_im = ifft2(F_imdata_sp)/ft_weight;
imshow(new_im);
