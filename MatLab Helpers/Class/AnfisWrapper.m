classdef AnfisWrapper
    %ANFISWRAPPER Wrapper and useful methods for anfis() related
    %operations.
    
    properties
        
        anfis_options
        initial_fis
        classification
        
        trained_ninputs
        trained_inputs_indices = NaN
        trained_fis
        
        train_error
        validation_error
        
    end
    
    methods
        function obj = AnfisWrapper( initial_fis, validation, epoch, classification )
           
            obj.anfis_options = anfisOptions( 'InitialFIS', initial_fis );
            obj.anfis_options.ValidationData = validation;
            obj.anfis_options.EpochNumber = epoch;
            
            obj.initial_fis = initial_fis;
            if nargin == 4 && classification
                
                obj.classification = true;
                
                % Fix ( for classification )
                for i = 1 : length( obj.initial_fis.output.mf )
                    
                    obj.initial_fis.output.mf( i ).type = 'constant';
                    obj.initial_fis.output.mf( i ).params = ...
                        obj.initial_fis.output.mf( i ).params( end );
                    
                end
                
            else
                
                obj.classification = false;
                
            end
            
        end
        
        function obj = disableDisplay(obj)
            
            obj.anfis_options.DisplayANFISInformation = 0;
            obj.anfis_options.DisplayErrorValues = 0;
            obj.anfis_options.DisplayStepSize = 0;
            obj.anfis_options.DisplayFinalResults = 0;
            
        end
        
        function obj = train( obj, training )
            
            [~, trainError, ~, chkFIS, chkError] = ...
                anfis( training, obj.anfis_options );
            
            obj.trained_ninputs = size( training, 2 ) - 1;
            
            obj.trained_fis = chkFIS;            
%             if obj.classification
%                 
%                 % Fix ( for classification )
%                 for i = 1 : length( obj.trained_fis.output.mf )
%                     obj.trained_fis.output.mf( i ).type = 'constant';
%                     obj.trained_fis.output.mf( i ).params = ...
%                         obj.trained_fis.output.mf( i ).params( end );
%                 end
%                 
%             end
            
            obj.train_error = trainError;
            obj.validation_error = chkError;
            
            if obj.trained_ninputs > 9
                obj.trained_inputs_indices = ...
                    sort( randperm( obj.trained_ninputs, 9 ) );
            end
            
        end
        
        function [testing_output, metrics] = test( obj, testing )
            
            % Eval FIS
            testing_output = evalfis( testing(:, 1:end-1),...
                obj.trained_fis );
            
            % If classification task, round output to closest integer index
            % ( index of closest class )
            if obj.classification
                testing_output = round( testing_output );
            end
            
            % Get Performance Metrics
            metrics = Metrics( testing(:,end), testing_output, ...
                obj.classification );
            
        end
        
        function plot_initial_mfs( obj )
            
            AnfisWrapper.plot_mfs( obj.initial_fis, ...
                obj.trained_ninputs, obj.trained_inputs_indices );
            
        end
        
        function plot_trained_mfs( obj )
            
            AnfisWrapper.plot_mfs( obj.trained_fis, ...
                obj.trained_ninputs, obj.trained_inputs_indices );
            
        end
        
        function plot_learning_curves( obj, modelIndex )
            
            i = 1 : length( obj.train_error );
            plot( i, obj.train_error, i, obj.validation_error );
            
            if nargin > 1 && modelIndex > 0
                title( ['MODEL ' num2str(modelIndex) ' | Learning Curves'] )
            else
                title( 'Learning Curves' )  
            end
            xlabel( '# of Iterations' )
            ylabel( 'MSE' )
            legend( 'Training','Validation' )
            
        end
    end
    
    methods ( Static )
        
        function plot_mfs( fis, n_inputs, input_indices )
            
            if nargin < 3
                input_indices = NaN;
            end
            
            if round( sqrt( n_inputs) )^2 ~= n_inputs
                if n_inputs < 4
                    n_rows = 2;
                    n_cols = 2;
                else
                    % max no. of sub plots is 9 ( 3x3 grid )
                    if ( nargin < 3 || any( isnan( input_indices ) ) )
                        input_indices = sort( randperm( n_inputs, 9 ) );
                    end
                    n_rows = 3;
                    n_cols = 3;
                end
            else
                n_rows = sqrt( n_inputs );
                n_cols = sqrt( n_inputs );
            end
            
            if any( isnan( input_indices ) )
                input_indices = 1 : n_inputs;
            end
            
            for input_i = 1 : ( n_rows * n_cols )
        
                input_label_i = input_indices( input_i );
                
                subplot( n_rows, n_cols, input_i )
                plotmf( fis, 'input', input_label_i )
                title( ['Input ' num2str(input_label_i)] )

            end
            
        end
       
        function fis = initial_fis_gp(nmfs, input_mf_type, output_mf_type, training)
            %ANFISWRAPPER.INITIAL_FIS_GP Get initial FIS for anfis() using 
            %grid partitioning as input space splitting method.
            %
            genOpt = genfisOptions( 'GridPartition' );
            genOpt.NumMembershipFunctions = nmfs;
            genOpt.InputMembershipFunctionType = input_mf_type;
            genOpt.OutputMembershipFunctionType = output_mf_type;
            
            fis = genfis( training(:, 1:end-1), training(:, end), genOpt );
            
        end
        
        function fis = initial_fis_sc( rad, training )
            %ANFISWRAPPER.INITIAL_FIS_SC Get initial FIS for anfis() using 
            %substractive clustering as input space splitting method.
            %
            
            genOpt = genfisOptions( 'SubtractiveClustering' );
            genOpt.ClusterInfluenceRange = rad;
%             genOpt.DataScale = 'auto';

            fis = genfis( training(:,1:end-1), training(:,end), genOpt );
       
        end
        
        
        function fis = initial_fis_fcm( nr, training )
            %ANFISWRAPPER.INITIAL_FIS_FCM Get initial FIS for anfis() using 
            %Fuzzy C-means as input space splitting method.
            %
            
            genOpt = genfisOptions( 'FCMClustering' );
            genOpt.FISType = 'sugeno';
            genOpt.NumClusters = nr;

            fis = genfis( training(:,1:end-1), training(:,end), genOpt );
       
        end
        
        function [training, validation, testing, probs] = partition_cl( dataset, split )
            %ANFISWRAPPER.PARTITION Partitions dataset in training, 
            %validation and testing datasets based on SPLIT. Each subset
            %should contain the same frequencies of classes across its
            %points, so we sort by class before calling actual partition()
            %method.
            %
            
            ncols = size( dataset, 2 );
            
            arr = grpstats( table( dataset(:,1), dataset(:,end) ), 'Var2' );
            dataset_freqs = arr.GroupCount ./ length( dataset );
            
            [training, validation, testing] = AnfisWrapper.partition( ...
                sortrows( dataset, ncols ), split );
            
            arr = grpstats( table( training(:,1), training(:,end) ), 'Var2' );
            training_freqs = arr.GroupCount ./ length( training );
            
            arr = grpstats( table( validation(:,1), validation(:,end) ), 'Var2' );
            validation_freqs = arr.GroupCount ./ length( validation );

            arr = grpstats( table( testing(:,1), testing(:,end) ), 'Var2' );
            testing_freqs = arr.GroupCount ./ length( testing );
            
            % Format probabilities table
            probs = table( categorical( ...
                1:length(dataset_freqs))', dataset_freqs, ...
                training_freqs, validation_freqs, testing_freqs, ...
                'VariableNames', {'Class', 'Dataset' 'Training' ...
                'Validation' 'Testing'} ...
            );
            
        end
        
        function [training, validation, testing] = partition(dataset, split)
            %ANFISWRAPPER.PARTITION Partitions dataset in training, 
            %validation and testing datasets based on SPLIT.
            %
            %   ANFIS.PARTITION(dataset, [0.6 0.2 0.2]) splits dataset in
            %   training set ( 60% ), validation and testing sets each
            %   having appr. 20% of the rows. The assignment of elements 
            %   in each sub-set is done in round robin manner, assigning 
            %   6-2-2 number of elements sequentially for every 10 elements 
            %   of dataset.
            %

            features_length = size( dataset, 2 );
            dataset_length = length( dataset );
            split_lengths = floor( dataset_length * split );

            training = zeros( split_lengths(1), features_length );
            validation = zeros( split_lengths(2), features_length );
            testing = zeros( split_lengths(3), features_length );

            tr = 1;
            vl = 1;
            ts = 1;
            i = 1;
            split10 = 10 * split;
            while i <= dataset_length

                training( tr : tr + split10(1), : ) = dataset( i : i + split10(1), : );
                tr = tr + split10(1);
                i = i + split10(1);

                validation( vl : vl + split10(2), : ) = dataset( i : i + split10(2), : );
                vl = vl + split10(2);
                i = i + split10(2);

                testing( ts : ts + split10(3), : ) = dataset( i : i + split10(3), : );
                ts = ts + split10(3);
                i = i + split10(3);

                if ( i + sum( split10 ) > dataset_length )
                    while( i <= dataset_length )

                        i_init = i;

                        if ( tr <= split_lengths(1) )
                            training( tr, : ) = dataset( i, : );
                            tr = tr + 1;
                            i = i + 1;
                        end

                        if ( vl <= split_lengths(2) )
                            validation( vl, : ) = dataset( i, : );
                            vl = vl + 1;
                            i = i + 1;
                        end

                        if ( ts <= split_lengths(3) )
                            testing( ts, : ) = dataset( i, : );
                            ts = ts + 1;
                            i = i + 1;
                        end

                        if ( i == i_init )
                            break
                        end

                    end

                    break
                end

            end

        end
        
    end
    
end

