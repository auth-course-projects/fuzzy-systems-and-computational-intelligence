clear, clc

% Specifications
split = [0.6 0.2 0.2];
NF = [5,10,15,20];           % Feature sub-set length
NR = [5,10,15,20,25];        % Nuber of Fuzzy Rules
K_FOLD = 5;

NFxNR = cartprod( NF, NR );
NFxNR_length = length( NF ) * length( NR );

%% Load & Pre-process dataset
load( 'isolet.dat', 'isolet' )
dataset = unique( isolet, 'rows' );

%% Split Dataset ( 60-20-20 split )
[training, validation, testing, probs] = ...
    AnfisWrapper.partition_cl( dataset, split );
clear dataset isolet

% disp( probs )
% disp( boxdist( [probs.Training, probs.Validation, probs.Testing] ) )

%% Get optimal (nf,nr) combination using Grid Search
cv = cvpartition( length( training ), 'k', K_FOLD );
validation_errors_cv = zeros( cv.NumTestSets, NFxNR_length );

% Cross Validation
for cv_test_i = 1 : cv.NumTestSets
    
    trn_index = cv.training( cv_test_i );
    val_index = cv.test( cv_test_i );
    
    % Get training/validation sets for this cross validation test
    cv_training = training( trn_index == 1, : );
    cv_validation = training( val_index == 1, : );
    
    % Will use SubstractiveClusteringWrapper to perform feature importance
    scw_cv_i = SubstractiveClusteringWrapper( cv_training, cv_validation, [], true );
    
    for i = 1:NFxNR_length
        
        grid_point( : ) = NFxNR( i, : );
        nf = grid_point(1);
        nr = grid_point(2);
       
        % Train model using SC
        initial_fis = AnfisWrapper.initial_fis_fcm( nr, cv_training( :, [scw_cv_i.nf2indices(nf) end] ) );
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
nf_opt = NFxNR( opt_grid_point_index, 1 );
nr_opt = NFxNR( opt_grid_point_index, 2 );

% 3D Stem
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
scw = SubstractiveClusteringWrapper( training, validation, testing, true );

%% Train model using optimum params
initial_fis = AnfisWrapper.initial_fis_fcm( nr_opt, training( :, [scw.nf2indices(nf_opt) end] ) );
model = AnfisWrapper( initial_fis, validation( :, [scw.nf2indices(nf_opt) end] ), 100 );
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