function H = myMean1D(X,W,y)
%% Weighted 1D histogram
% function H = myMean1D(X,W,y,z)
%
% Inputs: X - Nx1 data
%         W - Nx1 weights
%         y - vector of M(1) bin centers along 1st dim of X
% Output: H - M(1)xM(2) bar heights
% 
% Grid vectors (y, z) must have linear spacing


%# bin centers (integers)
xNumBins = numel(y); 

%# map X/Y values to bin indices
Xi = round( interp1(y, 1:xNumBins, X, 'linear', 'extrap') );

%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);

%# count number of elements in each bin
H = accumarray(Xi(:), W(:), [xNumBins 1]);
norm = 1; %accumarray(Xi(:), 1, [xNumBins 1]);

H=(H./norm)';

end



