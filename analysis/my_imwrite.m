function my_imwrite(A,cmap,clim,savname)
%MY_IMWRITE save png image with fixed color scale
%   Detailed explanation goes here

    Amax = clim(2);
    Amin = clim(1);
    A = (A - Amin)/(Amax-Amin);
    A(A<0) = 0;
    A(A>1) = 1;
    Ncolor = length(cmap);
    A = uint8(A*Ncolor);
    imwrite(A, cmap, savname, 'png');

end

