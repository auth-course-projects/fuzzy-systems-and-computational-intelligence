clear, clc

% Specifications
split = [0.6 0.2 0.2];
NR = [5,8,12,16,20];        % Nuber of Fuzzy Rules per model variant
RAD = [1, 0.853125, 0.7, 0.5875, 0.525];

nf = -1;    % include all features

%% Load & Pre-process dataset
load( 'avila.txt', 'avila' )
dataset = unique( avila, 'rows' );
%   - apply smoothing
smoothing_factors = [.5, .75, .5, .25, .5, .5, .5, .5, .5, .5];
for f_i = 1 : size( dataset, 2 ) - 1
    dataset( :, f_i ) = smoothdata( dataset( :, f_i ), ...
        'SmoothingFactor', 1 );
end
%   - apply normalization
% dataset = [normalize( dataset( :, 1:end-1 ) ) dataset( :, end )];

%% Split Dataset ( 60-20-20 split )
[training, validation, testing, ~] = ...
    AnfisWrapper.partition_cl( dataset, split );
clear avila

% disp( probs )
% disp( boxdist( [probs.Training, probs.Validation, probs.Testing] ) )

%% Get Rads from NR ( for Substractive Clustering )
scw = SubstractiveClusteringWrapper( training, validation, testing, true );
if 0 == exist( 'RAD', 'var' )
    
    RAD = zeros( size( NR ) );
    for nr_i = 1 : length( NR )

        nr = NR( nr_i );
        [nr_hat, rad] = scw.getFeasibleNfNrCombination( nf, nr );

        if nr_hat == nr
            RAD( nr_i ) = rad;
        else
            warning( ['Could not converge to nr=' num2str( nr ) ...
                '. Closest nr_hat = ' num2str( nr_hat )] )
            return
        end

    end
    
end

%% Train models
metrics = Metrics.empty( length( RAD ), 0 ); 
for rad_i = 1 : length( RAD )
    
    initial_fis = AnfisWrapper.initial_fis_sc( RAD( rad_i ), training( :, [scw.nf2indices(nf) end] ) );    
    model = AnfisWrapper( initial_fis, validation( :, [scw.nf2indices(nf) end] ), 100, true );
    model = model.train( training( :, [scw.nf2indices(nf) end] ) );     
    
    %% Test trained model and get metrics
    [testing_output, metrics( rad_i )] = model.test( testing( :, [scw.nf2indices(nf) end] ) );

    %% Learning Curves
    %   - trained mfs for each input
    figure
    model.plot_initial_mfs();
    suptitle( ['MODEL ' num2str(rad_i) ' | Initial Input MFs'] )
    figure
    model.plot_trained_mfs();
    suptitle( ['MODEL ' num2str(rad_i) ' | Trained Input MFs'] )

    %   - learning curves ( training/validation error )
    figure
    model.plot_learning_curves();

    %   - prediction error
    figure
    plot( testing_output - testing( :, end ) );
    title( ['MODEL ' num2str(rad_i) ' | Prediction Errors'] )
    xlabel( 'sample index' )
    ylabel( 'error' )
    
end
    
    
