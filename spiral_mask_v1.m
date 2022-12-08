%%% MRI Image
imdata = phantom('Modified Shepp-Logan', 256);

%%% Create Probability Density Function - PDF
PDF = create_PDF(imdata);
% figure; imshow(PDF); title("PDF");

rows = size(imdata,1);
cols = size(imdata,2);

% sz0 = 54;
% imdata = padarray(imdata,sz0,0, 'both');
% imdata = padarray(imdata',sz0,0, 'both');
% imdata = imdata';

ft_weight = 1/sqrt(size(imdata,1)*size(imdata,2));
F_imdata = fftshift(fft2(imdata).*ft_weight);


figure; 
subplot(1,2,1); imshow(abs(imdata)); title('Shepp-Logan Image');
subplot(1,2,2); imshow(abs(F_imdata)); title('FT Shepp-Logan');




%%
%sz1 = 128;
sz1 = 182;

%Nitlv = 32;
Nitlv = 45.5;
fov = [sz1,sz1];
radius = [0,1];

kmax = 1;

if length(radius)<2
	error(' radius must be at least length=2');
end

dr = 1/1500/max(fov/Nitlv)*4;
r = 0:dr:kmax;
kmax = max(r);

fov =  interp1(radius*kmax, fov, r, 'spline');
dtheta = 2*pi*dr.*fov/Nitlv;
theta = cumsum(dtheta);

i = sqrt(-1);
C = r.*exp(i*theta);


figure; plot(C*sz1, '.');

% sz2 = 128;
sz2 = 182;

% figure; hold on;
% plot(C*128, '.');
% plot(C*exp(-i*pi/2)*sz2, '.');
% plot(C*exp(-i*pi)*sz2, '.');
% plot(C*exp(-i*3*pi/2)*sz2, '.');
% 
% plot(C*exp(-i*1*pi/4)*sz2, '.');
% plot(C*exp(-i*3*pi/4)*sz2, '.');
% plot(C*exp(-i*5*pi/4)*sz2, '.');
% plot(C*exp(-i*7*pi/4)*sz2, '.');
% 
% plot(C*exp(-i*1*pi/6)*sz2, '.');
% plot(C*exp(-i*5*pi/6)*sz2, '.');
% plot(C*exp(-i*7*pi/6)*sz2, '.');
% plot(C*exp(-i*11*pi/6)*sz2, '.');
% 
% plot(C*exp(-i*1*pi/3)*sz2, '.');
% plot(C*exp(-i*2*pi/3)*sz2, '.');
% plot(C*exp(-i*4*pi/3)*sz2, '.');
% plot(C*exp(-i*5*pi/3)*sz2, '.');
% 
% hold off;


%%
scale = 182; %128

C_phase = zeros(length(C), 4);
C_phase(:,1) = C*scale;
C_phase(:,2) = C*exp(-i*pi/2)*scale;
C_phase(:,3) = C*exp(-i*pi)*scale;
C_phase(:,4) = C*exp(-i*3*pi/2)*scale;

C_phase(:,5) = C*exp(-i*1*pi/4)*scale;
C_phase(:,6) = C*exp(-i*3*pi/4)*scale;
C_phase(:,7) = C*exp(-i*5*pi/4)*scale;
C_phase(:,8) = C*exp(-i*7*pi/4)*scale;

C_phase(:,9) = C*exp(-i*1*pi/6)*scale;
C_phase(:,10) = C*exp(-i*5*pi/6)*scale;
C_phase(:,11) = C*exp(-i*7*pi/6)*scale;
C_phase(:,12) = C*exp(-i*11*pi/6)*scale;

C_phase(:,13) = C*exp(-i*1*pi/3)*scale;
C_phase(:,14) = C*exp(-i*2*pi/3)*scale;
C_phase(:,15) = C*exp(-i*4*pi/3)*scale;
C_phase(:,16) = C*exp(-i*5*pi/3)*scale;

C_phase(:,17) = C*exp(-i*1*pi/12)*scale;
C_phase(:,18) = C*exp(-i*11*pi/12)*scale;
C_phase(:,19) = C*exp(-i*13*pi/12)*scale;
C_phase(:,20) = C*exp(-i*23*pi/12)*scale;

C_phase(:,21) = C*exp(-i*3*pi/7)*scale;
C_phase(:,22) = C*exp(-i*4*pi/7)*scale;
C_phase(:,23) = C*exp(-i*10*pi/7)*scale;
C_phase(:,24) = C*exp(-i*11*pi/7)*scale;

C_phase(:,25) = C*exp(-i*3*pi/8)*scale;
C_phase(:,26) = C*exp(-i*5*pi/8)*scale;
C_phase(:,27) = C*exp(-i*11*pi/8)*scale;
C_phase(:,28) = C*exp(-i*13*pi/8)*scale;

C_phase(:,29) = C*exp(-i*1*pi/7)*scale;
C_phase(:,30) = C*exp(-i*6*pi/7)*scale;
C_phase(:,31) = C*exp(-i*8*pi/7)*scale;
C_phase(:,32) = C*exp(-i*13*pi/7)*scale;



rC_phase = round(real(C_phase));
iC_phase = -round(imag(C_phase));

mask1 = zeros(scale*2,scale*2);

counter = 0;
for n = 1:size(C_phase, 2)
    for m = 1:size(C_phase,1)
        mask1(iC_phase(m,n)+1+scale, rC_phase(m,n)+1+scale) = 1;
        counter = counter + 1;
    end
end

%figure; imagesc(mask1); colormap('gray');


%
mask2 = mask1(55:310, 55:310);
figure; imagesc(mask2); colormap('gray');



[N1, N2] = size(F_imdata);
num_of_orig_data = N1*N2;
num_of_spiral = sum(sum(mask2));

percent_of_points = num_of_spiral/num_of_orig_data *100

%%

F_imdata_masked = F_imdata .* mask2;

figure; imshow(abs(F_imdata_masked)); title('FT Shepp-Logan');


imdata_masked = ifft2(ifftshift(F_imdata_masked))./ft_weight;

figure; imagesc(abs(imdata_masked)); title('Masked Shepp-Logan');


