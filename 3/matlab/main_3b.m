clearvars -except NFxNRxRAD
clc

% Specifications
split = [0.6 0.2 0.2];
NF = [3,9,15,21];           % Feature sub-set length
NR = [4,8,12,16,20];        % Nuber of Fuzzy Rules
NFxNR_length = length( NF ) * length( NR );
K_FOLD = 5;

%% Load & Pre-process dataset
load( 'superconduct.csv', 'superconduct' )
dataset = unique( superconduct, 'rows' );

%% Split Dataset ( 60-20-20 split )
[training, validation, testing] = AnfisWrapper.partition( dataset, split );
clear dataset superconduct

%% Get optimal (nf,nr) combination using Grid Search
cv = cvpartition( length( training ), 'k', K_FOLD );
if 0 == exist( 'NFxNRxRAD', 'var' )
    NFxNRxRAD = zeros( NFxNR_length, 3, cv.NumTestSets );
end
checkRADSFirst = zeros( NFxNR_length, 1 );
validation_errors_cv = zeros( cv.NumTestSets, NFxNR_length );

% Cross Validation
for cv_test_i = 1 : cv.NumTestSets
    
    trn_index = cv.training( cv_test_i );
    val_index = cv.test( cv_test_i );
    
    % Get training/validation sets for this cross validation test
    cv_training = training( trn_index == 1, : );
    cv_validation = training( val_index == 1, : );
    
    % Get NFxNRxRAD for this test
    scw_cv_i = SubstractiveClusteringWrapper( cv_training, cv_validation, [] );
    
    if cv_test_i == 1 && 0 == exist( 'NFxNRxRAD', 'var' )
        NFxNRxRAD( :, :, cv_test_i ) = ...
            scw_cv_i.getFeasibleNFxNRCombinations( NF, NR );
    else
        if cv_test_i == 1
            checkRADSFirst( : ) = NFxNRxRAD( :, 3, 1 );
        else
            checkRADSFirst( : ) = NFxNRxRAD( :, 3, cv_test_i - 1 );
        end
        
        NFxNRxRAD( :, :, cv_test_i ) = ...
            scw_cv_i.getFeasibleNFxNRCombinations( NF, NR, checkRADSFirst );
    end
    
    for i = 1:NFxNR_length
        
        grid_point( : ) = NFxNRxRAD( i, :, cv_test_i );
        nf = grid_point(1);
        rad = grid_point(3);
       
        % Train model using SC
        initial_fis = AnfisWrapper.initial_fis_sc( rad, cv_training( :, [scw_cv_i.nf2indices(nf) end] ) );
        model = AnfisWrapper( initial_fis, cv_validation( : , [scw_cv_i.nf2indices(nf) end] ), 20 );
        model = model.disableDisplay();
        model = model.train( cv_training( :, [scw_cv_i.nf2indices(nf) end] ) );

        % Mean validation error
        validation_errors_cv( cv_test_i, i ) = mean( model.validation_error );
        
    end
    
end

clear scw_cv_i cv_training cv_validation

% Get grid point that lead to minimum validation error
validation_errors = mean( validation_errors_cv );
[~, opt_grid_point_index] = min( validation_errors );

% Get optimum parameters
nf_opt = NFxNRxRAD( opt_grid_point_index, 1, 1 );
nr_opt = NFxNRxRAD( opt_grid_point_index, 2, 1 );
rads_opt( : ) = NFxNRxRAD( opt_grid_point_index, 3, : );
rad_opt = trimmean( rads_opt( rads_opt > 0 ), 50  );

% 3D Stem
NFxNR = cartprod( NF, NR );
validation_errors_ij = zeros( length( NF ), length( NR ) );
for i = 1:NFxNR_length
    j = ceil( i / length( NF ) );
    validation_errors_ij( i - (j-1)*length( NF ), j ) = validation_errors( i );
end
figure
stem3( NR, NF, validation_errors_ij,':*r' )
hold on
stem3( nr_opt, nf_opt, validation_errors_ij( find( NF == nf_opt ), find( NR == nr_opt ) ), ':*g' )
hold off
 
% Keep only the "best" features of training/validation/testing set
scw = SubstractiveClusteringWrapper( training, validation, testing );

%% Train model using optimum params
initial_fis = AnfisWrapper.initial_fis_sc( rad_opt, training( :, [scw.nf2indices(nf_opt) end] ) );
model = AnfisWrapper( initial_fis, validation( :, [scw.nf2indices(nf_opt) end] ), 60 );
model = model.train( training( :, [scw.nf2indices(nf_opt) end] ) );

%% Test trained model and get metrics
[testing_output, metrics] = model.test( testing( :, [scw.nf2indices(nf_opt) end] ) );

%% Learning Curves
%   - trained mfs for each input
figure
model.plot_initial_mfs();
suptitle( 'Initial Input MFs' )
figure
model.plot_trained_mfs();
suptitle( 'Trained Input MFs' )

%   - learning curves ( training/validation error )
figure
model.plot_learning_curves();

%   - prediction error
figure
plot( testing_output - testing( :, end ) );
title( 'Prediction Errors' )
xlabel( 'sample index' )
ylabel( 'error' )