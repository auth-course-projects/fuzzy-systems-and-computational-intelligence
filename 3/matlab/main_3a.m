clear, clc

% Specifications
split = [0.6 0.2 0.2];
tsk_param_nmfs = [2 3 2 3];
tsk_param_types = ["constant" "constant" "linear" "linear"];
mf_types = ["gbellmf", "gbellmf", "gbellmf", "gbellmf"];
n_models = length( tsk_param_nmfs );

%% Load & Pre-process dataset
load('CCPP.dat', 'CCPP')
dataset = unique( CCPP, 'rows' );
n_inputs = size( dataset, 2 ) - 1;

%% Split Dataset ( 60-20-20 split )
[training, validation, testing] = AnfisWrapper.partition( dataset, split );
clear dataset CCPP

%% Training
metrics = Metrics.empty( n_models, 0 );
for i = n_models
    
    initial_fis = AnfisWrapper.initial_fis_gp( tsk_param_nmfs(i), ...
        mf_types(i), tsk_param_types(i), training );
    
    model = AnfisWrapper( initial_fis, validation, 80 );
    model = model.train( training );
    
    %% Test trained model and get metrics
    [testing_output, metrics] = model.test( testing );
    
    %% Learning Curves
    %   - trained mfs for each input
    figure
    model.plot_initial_mfs();
    suptitle( ['MODEL ' num2str( i ) ' | Initial Input MFs'] )
    figure
    model.plot_trained_mfs();
    suptitle( ['MODEL ' num2str( i ) ' | Trained Input MFs'] )
    
    %   - learning curves ( training/validation error )
    figure
    model.plot_learning_curves( i );

    %   - prediction error
    figure
    plot( testing_output - testing( :, end ) );
    title( ['MODEL ' num2str(i) ' | Prediction Errors'] )
    xlabel( 'sample index' )
    ylabel( 'error' )
    
end