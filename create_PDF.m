function pdf = create_PDF(input_image)
    rows = size(input_image, 1);
    cols = size(input_image, 2);
    [x, y] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
    partial_sampling_factor = 0.5;
    
    r = sqrt(x.^2 + y.^2);
    r = r/max(abs(r(:)));
    pdf = (1-r).^2;
    
    minval = 0;
    maxval = 1;
    val = 0.5;
    PCTG = floor(partial_sampling_factor*rows*cols);
    
    while true
	    val = minval/2 + maxval/2;
	    pdf = (1-r).^2 + val; 
        pdf(pdf>1) = 1; 
    
	    N = floor(sum(pdf(:)));
	    if N > PCTG % infeasible
		    maxval = val;
	    end
	    if N < PCTG % not optimal
		    minval = val;
	    end
	    if N == PCTG % optimal
		    break;
	    end
    end

end