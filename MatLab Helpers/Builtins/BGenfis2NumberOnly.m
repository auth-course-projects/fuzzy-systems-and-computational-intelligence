function numRule = BGenfis2NumberOnly(Xin,Xout,radii,xBounds,options,user_centers)
%BGENFIS2 Generates a Sugeno-type FIS using subtractive clustering.
%
%   Given separate sets of input and output data, GENFIS2 generates a fuzzy
%   inference system (FIS) using fuzzy subtractive clustering. GENFIS2 can be
%   used to generate an initial FIS for ANFIS training by first applying
%   subtractive clustering on the data. GENFIS2 accomplishes this by extracting
%   a set of rules that models the data behavior. The rule extraction method
%   first uses the SUBCLUST function to determine the number of rules and
%   antecedent membership functions and then uses linear least squares
%   estimation to determine each rule's consequent equations.
%
%   FIS = GENFIS2(XIN,XOUT,RADII) returns a Sugeno-type FIS given input data XIN
%   and output data XOUT. The matrices XIN and XOUT have one column per FIS 
%   input and output, respectively. RADII specifies the range of influence of
%   the cluster center for each input and output dimension, assuming the data 
%   falls within a unit hyperbox (range [0 1]). Specifying a smaller cluster
%   radius will usually yield more, smaller clusters in the data, and hence more
%   rules. When RADII is a scalar it is applied to all input and output
%   dimensions. When RADII is a vector, it has one entry for each input and
%   output dimension.
%
%   FIS = GENFIS2(...,XBOUNDS) also specifies a matrix XBOUNDS used to normalize
%   the data XIN and XOUT into a unit hyperbox (range [0 1]). XBOUNDS is size 
%   2-by-N, where N is the total number of inputs and outputs. Each column of
%   XBOUNDS provides the minimum and maximum values for the corresponding input
%   or output data set. If XBOUNDS is an empty matrix or not provided, the
%   minimum and maximum data values found in XIN and XOUT, are used as defaults.
%
%   FIS = GENFIS2(...,XBOUNDS,OPTIONS) specifies options for changing the
%   default algorithm parameters, type HELP SUBCLUST for details.
%
%   FIS = GENFIS2(...,XBOUNDS,OPTIONS,USER_CENTERS) accepts user-supplied
%   cluster centers.  USER_CENTERS has size J-by-N where J is the number of
%   clusters and N is the total number of inputs and outputs.
%
%   Examples
%       Xin1 =  7*rand(50,1);
%       Xin2 = 20*rand(50,1)-10;
%       Xin  = [Xin1 Xin2];
%       Xout =  5*rand(50,1);
%       fis = genfis2(Xin,Xout,0.5) specifies a range of influence of 0.5 for
%       all data dimensions.
%
%       fis = genfis2(Xin,Xout,[0.5 0.25 0.3]) specifies the ranges of influence
%       in the first, second, and third data dimensions are 0.5, 0.25, and 0.3
%       times the width of the data space, respectively.
%
%       fis = genfis2(Xin,Xout,0.5,[0 -10 0; 7 10 5]) specifies the data in
%       the first column of Xin are scaled from [0 7], the data in the
%       second column of Xin are scaled from [-10 10], and the data in Xout are
%       scaled from [0 5].
%
%   See also SUBCLUST, GENFIS1, ANFIS.

%   Copyright 1994-2018 The MathWorks, Inc.

%   Reference
%   S. Chiu, "Fuzzy Model Identification Based on Cluster Estimation," J. of
%   Intelligent & Fuzzy Systems, Vol. 2, No. 3, 1994.

    warning(message('fuzzy:general:warnGenfis2_Deprecation'))

    [numData,numInp] = size(Xin);
    [numData2,numOutp] = size(Xout);

    if numData ~= numData2
        % There's a mismatch in the input and output data matrix dimensions
        if numData == numOutp
            % The output data matrix should have been transposed, we'll fix it
            Xout = Xout';
            numOutp = numData2;
        else
            error(message("fuzzy:general:errGenfis2_XinXoutSizeMismatch"))
        end
    end

    if nargin < 5
        verbose = 0;
        if nargin < 4
            xBounds = [];
        end
        [centers,sigmas] = subclust([Xin Xout],radii,xBounds);
    elseif nargin < 6
        verbose = options(4);
        [centers,sigmas] = subclust([Xin Xout],radii,xBounds,options);
    else
        if isempty(options)
            verbose = 0;
        else
            verbose = options(4);
        end
        numParams = numInp + numOutp;
        if numParams~=size(user_centers,2)
            error(message("fuzzy:general:errGenfis2_InvalidClusterCenterDataSize", ...
                numParams))
        end
        centers = user_centers;
        % if only one value is given as the range of influence, then apply that
        % value to all data dimensions
        if length(radii) == 1 && numParams ~= 1
            radii = radii * ones(1,numParams);
        end
        % distance multipliers for accumulating and squashing cluster potentials
        X = [Xin Xout];
        if isempty(xBounds)
            % No data scaling range values are specified, use the actual minimum and maximum values
            minX = min(X);
            maxX = max(X);
            % If the actual min and max values have a range of zero, calculate a small range
            % relative to the data, to allow a sigma value to be calculated.
            index = find(maxX == minX);
            minX(index) = minX(index) - 0.0001*(1 + abs(minX(index)));
            maxX(index) = maxX(index) + 0.0001*(1 + abs(maxX(index)));
        else
            % Use the user supplied data range values in xBounds
            minX = xBounds(1,:);
            maxX = xBounds(2,:);
            % Verify correct dimensions and values for xBounds were supplied
            if length(minX) ~=  size(X,2)
                error(message("fuzzy:general:errGenfis2_InvalidXBoundSize"))
            elseif any(maxX == minX)
                error(message("fuzzy:general:errGenfis2_ZeroXBoundDataRange"))
            end
        end
        sigmas = (radii .* (maxX - minX)) / sqrt(8.0);
    end

    if verbose
        disp('Setting up matrix for linear least squares estimation...');
    end

    % Discard the clusters' output dimensions
    centers = centers(:,1:numInp);
    sigmas = sigmas(:,1:numInp);

    % Distance multipliers
    distMultp = (1.0 / sqrt(2.0)) ./ sigmas;

    numRule = size(centers,1);
end

%% Helper functions -------------------------------------------------------
function range = getRange(data,bounds,n1,n2)
if isempty(bounds)
    % No data scaling range values were specified, use the actual minimum
    % and maximum values of the data.
    minData = min(data);
    maxData = max(data);
else
    minData = bounds(1,n1:n2);
    maxData = bounds(2,n1:n2);
end
range = [minData ; maxData]';
end