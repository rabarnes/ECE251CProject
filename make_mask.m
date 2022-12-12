function  [mask, percent_of_points] = make_mask(size, m)
    
    mask = ones(size,size);
    mask(1:m:round(4*size/9),:) = zeros(length(1:m:round(4*size/9)),size);
    mask(round(5*size/9):m:end,:) = zeros(length(round(5*size/9):m:size),size);

    num_of_orig_data = size*size;
    num_of_mask = sum(sum(mask));
    
    percent_of_points = num_of_mask/num_of_orig_data *100;
end