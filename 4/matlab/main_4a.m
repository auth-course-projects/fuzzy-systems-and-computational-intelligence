clearvars -except RAD
clc

% Specifications
split = [0.6 0.2 0.2];
NR = [4,8,12,16,20];        % Nuber of Fuzzy Rules per model variant
% RAD = [0.1875, 0.10546875, 0.08125, 0.0625, 0.0546875];

NR = 8;
% RAD = 0.10546875;
nf = -1;

%% Load & Pre-process dataset
load( 'avila.txt', 'avila' )
dataset = unique( avila, 'rows' );
% dataset = [smoothdata( dataset( :, 1:end-1 ), 'SmoothingFactor', 1 ) ...
%     dataset( :, end )];

% smoothing_factors = [.5, .75, .5, .3, .1, .2, .3, .4, .5, .5];
% for f_i = 1 : size( dataset, 2 ) - 1
%    
%     dataset( :, f_i ) = smoothdata( dataset( :, f_i ), ...
%         'SmoothingFactor', 1 );
%     
% end

% dataset = [smoothdata( dataset( :, 1:end-1 ) ) dataset(:,end)];
% dataset = [normalize( dataset( :, 1:end-1 ) ) dataset(:,end)];
% load( 'wifi-localization.dat', 'wifi_localization' )
% dataset = unique( wifi_localization, 'rows' );

%% Split Dataset ( 60-20-20 split )
[training, validation, testing, probs] = ...
    AnfisWrapper.partition_cl( dataset, split );
clear avila CCPP wifi_localization

disp( probs )
disp( boxdist( [probs.Training, probs.Validation, probs.Testing] ) )

%% Get Rads from NR ( for Substractive Clustering )
if 0 == exist( 'RAD', 'var' )
    
    RAD = zeros( size( NR ) );
    scw = SubstractiveClusteringWrapper( training, validation, testing );
    for nr_i = 1 : length( NR )

        nr = NR( nr_i );
        [nr_hat, rad] = scw.getFeasibleNfNrCombination( nf, nr, 0.2 );

        if nr_hat == nr
            RAD( nr_i ) = rad;
        else
            warning( ['Could not converge to nr=' num2str( nr ) ...
                '. Closest nr_hat = ' num2str( nr_hat )] )
            return
        end

    end

%     clear scw
    
end

%% Train models
metrics = Metrics.empty( length( RAD ), 0 ); 
for rad_i = 1 : length( RAD )
   
    initial_fis = AnfisWrapper.initial_fis_sc( RAD( rad_i ), training( :, [scw.nf2indices(nf) end] ) );    
    model = AnfisWrapper( initial_fis, validation( :, [scw.nf2indices(nf) end] ), 100, true );
    model = model.train( training( :, [scw.nf2indices(nf) end] ) );     
    
    %% Test trained model and get metrics
    [testing_output, metrics( rad_i )] = model.test( testing );

    %% Learning Curves
%     %   - trained mfs for each input
%     figure
%     model.plot_initial_mfs();
%     suptitle( ['MODEL ' num2str(rad_i) ' | Initial Input MFs'] )
%     figure
%     model.plot_trained_mfs();
%     suptitle( ['MODEL ' num2str(rad_i) ' | Trained Input MFs'] )
% 
    %   - learning curves ( training/validation error )
    figure
    model.plot_learning_curves();
% 
%     %   - prediction error
%     figure
%     plot( testing_output - testing( :, end ) );
%     title( ['MODEL ' num2str(rad_i) ' | Prediction Errors'] )
%     xlabel( 'sample index' )
%     ylabel( 'error' )
    
    break
    
end
    
    
