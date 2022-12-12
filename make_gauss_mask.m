function [gaussian_mask, percent_of_points] = make_gauss_mask(dim_im, percent_sample)

% dim_im = size(im);
num_pixels = dim_im * dim_im;
% percent_sample = 0.3;
num_subsample = round(num_pixels*percent_sample);
% random gaussian
u = [0 0];
sigma = [1 0; 0 1];
% generate subsamples of multivariate gaussian
R = mvnrnd(u, sigma, round(num_subsample));
% scale all values of R to indices of image
a = 1;
b = dim_im; % TODO fix this for rectangular images
min_R = min(R(:));
max_R = max(R(:));
R = ((b-a)*(R - min_R)) / (max_R-min_R) + a;
R = round(R);
% create mask by setting mask indices to 1
gaussian_mask = zeros([dim_im dim_im]);
for i=1:num_subsample
    gaussian_mask(R(i,1), R(i,2)) = 1;
end

num_of_orig_data = dim_im*dim_im;
num_of_gauss = sum(sum(gaussian_mask));
percent_of_points = num_of_gauss/num_of_orig_data *100;

end