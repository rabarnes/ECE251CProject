function [mask_out, percent_of_points] = make_spiral_mask(dim_im, density_knob)
    % Example: make_spiral_mask(128, 4)
    % density_knob - lower means higher density
    
    half_size = dim_im./2;
    circ_size = ceil(half_size*sqrt(2));
    
    Nitlv = circ_size/4;
    fov = [circ_size,circ_size];
    radius = [0,1];
    
    kmax = 1;
    dr = 1/1500/max(fov/Nitlv)*density_knob;
    r = 0:dr:kmax;
    kmax = max(r);
    
    fov =  interp1(radius*kmax, fov, r, 'spline');
    dtheta = 2*pi*dr.*fov/Nitlv;
    theta = cumsum(dtheta);
    
    i = sqrt(-1);
    C = r.*exp(i*theta);
    
    C_phase = zeros(length(C), 4);
    C_phase(:,1) = C*circ_size;
    C_phase(:,2) = C*exp(-i*pi/2)*circ_size;
    C_phase(:,3) = C*exp(-i*pi)*circ_size;
    C_phase(:,4) = C*exp(-i*3*pi/2)*circ_size;
    
    C_phase(:,5) = C*exp(-i*1*pi/4)*circ_size;
    C_phase(:,6) = C*exp(-i*3*pi/4)*circ_size;
    C_phase(:,7) = C*exp(-i*5*pi/4)*circ_size;
    C_phase(:,8) = C*exp(-i*7*pi/4)*circ_size;
    
    C_phase(:,9) = C*exp(-i*1*pi/6)*circ_size;
    C_phase(:,10) = C*exp(-i*5*pi/6)*circ_size;
    C_phase(:,11) = C*exp(-i*7*pi/6)*circ_size;
    C_phase(:,12) = C*exp(-i*11*pi/6)*circ_size;
    
    C_phase(:,13) = C*exp(-i*1*pi/3)*circ_size;
    C_phase(:,14) = C*exp(-i*2*pi/3)*circ_size;
    C_phase(:,15) = C*exp(-i*4*pi/3)*circ_size;
    C_phase(:,16) = C*exp(-i*5*pi/3)*circ_size;
    
    C_phase(:,17) = C*exp(-i*1*pi/12)*circ_size;
    C_phase(:,18) = C*exp(-i*11*pi/12)*circ_size;
    C_phase(:,19) = C*exp(-i*13*pi/12)*circ_size;
    C_phase(:,20) = C*exp(-i*23*pi/12)*circ_size;
    
    C_phase(:,21) = C*exp(-i*3*pi/7)*circ_size;
    C_phase(:,22) = C*exp(-i*4*pi/7)*circ_size;
    C_phase(:,23) = C*exp(-i*10*pi/7)*circ_size;
    C_phase(:,24) = C*exp(-i*11*pi/7)*circ_size;
    
    C_phase(:,25) = C*exp(-i*3*pi/8)*circ_size;
    C_phase(:,26) = C*exp(-i*5*pi/8)*circ_size;
    C_phase(:,27) = C*exp(-i*11*pi/8)*circ_size;
    C_phase(:,28) = C*exp(-i*13*pi/8)*circ_size;
    
    C_phase(:,29) = C*exp(-i*1*pi/7)*circ_size;
    C_phase(:,30) = C*exp(-i*6*pi/7)*circ_size;
    C_phase(:,31) = C*exp(-i*8*pi/7)*circ_size;
    C_phase(:,32) = C*exp(-i*13*pi/7)*circ_size;
    
    rC_phase = round(real(C_phase));
    iC_phase = -round(imag(C_phase));
    
    mask_temp = zeros(circ_size*2,circ_size*2);
    
    for n = 1:size(C_phase, 2)
        for m = 1:size(C_phase,1)
            mask_temp(iC_phase(m,n)+1+circ_size, rC_phase(m,n)+1+circ_size) = 1;
        end
    end
    
    mask_out = mask_temp(55:310, 55:310);
    
    num_of_orig_data = dim_im*dim_im;
    num_of_spiral = sum(sum(mask_out));
    
    percent_of_points = num_of_spiral/num_of_orig_data *100;
end