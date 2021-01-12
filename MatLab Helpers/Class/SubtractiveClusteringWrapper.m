classdef SubtractiveClusteringWrapper
    %SUBTRACTIVECLUSTERINGWRAPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
      
        FSS_ALGO = 'relief';    % 'relief', 'nca', 'ensemble'
        FSS_ALGO_RELIEF_K = 100;
        
        BISECTION_MAX_ITERS = 20;
        BISECTION_TOLERATE_DISTANCE = 2;    % +- 2 clusters
        
        % -----------------
        % Bisection Policy
        % -----------------
        %   - "tolerate": if max_iters reached and the distance of number
        %   of cluster to the asked number of clusters is less than
        %   tolerate distance, then it is assumed that the asked number of
        %   clusters has been found for the last computed radius.
        %   - "abort": throws exception and stops execution (the respective
        %   (nf,nr) combination is assumed not feasible
        ON_BISECTION_MAX_REACHED_POLICY_TOLERATE = 1;
        ON_BISECTION_MAX_REACHED_POLICY_ABORT = 2;
        ON_BISECTION_MAX_REACHED_POLICY = 1;    
        
    end
    
    properties
        
        training
        validation
        testing
        
        feature_indices
        NFMap
        use_builtins
        
    end
    
    methods
        function obj = SubtractiveClusteringWrapper( training, validation, testing, getFI, FI, use_builtins )
            %SUBSTRACTIVECLUSTERINGWRAPPER Construct an instance of this class
            %   Detailed explanation goes here
            
            % Save parts
            obj.training = training;
            obj.validation = validation;
            obj.testing = testing;
            
            % Get feature indices ordered by importance (descending order)
            if nargin == 6
                obj.use_builtins = use_builtins;
            elseif nargin >= 4 && getFI == true
                if ( nargin == 5 && isvector(FI) )
                    obj.feature_indices = FI;
                else
                    switch( obj.FSS_ALGO )
                        case 'relief'
                            obj.feature_indices = obj.getFeatureImportanceRL;
                        case 'nca'
                            obj.feature_indices = obj.getFeatureImportanceNCA;
                        otherwise 
                            obj.feature_indices = obj.getFeatureImportanceFE;
                    end
                end
            else
                obj.feature_indices = 1 : size( obj.training, 2 ) - 1;
            end
            
        end
        
        function feature_indices = getFeatureImportanceFE( obj )

            ens = fitrensemble( ...
                obj.training( :, 1:end-1 ), obj.training( :, end ), ...
                'Learners', templateTree( 'Reproducible', true ), ...
                'OptimizeHyperparameters', 'auto', ...
                'HyperparameterOptimizationOptions', struct('UseParallel',true) ...
            );
        
            feature_importance = predictorImportance( ens );
            [~, feature_indices] = sort( feature_importance, 2, 'descend' );

        end
        
        function feature_indices = getFeatureImportanceNCA( obj )
             
            mdl = fsrnca( ...
                obj.training( :, 1:end-1 ), obj.training( :, end ), ...
                'FitMethod', 'exact', 'Solver', 'lbfgs' ...
            );
        
            feature_importance = mdl.FeatureWeights;
            [~, feature_indices] = sort( feature_importance, 2, 'descend' );

        end
        
        function feature_indices = getFeatureImportanceRL( obj )
           
            [feature_indices, ~] = relieff( obj.training(:, 1:end-1), ...
                obj.training(:, end), obj.FSS_ALGO_RELIEF_K );

        end
        
        function feature_indices = nf2indices( obj, nf )
            %SubstractiveClusteringWrapper.N_FEATURES Uses ReliefF to get
            % NF most important features of dataset. Returns column indices
            
            if nargin == 2 && nf > 0
                feature_indices = obj.feature_indices( 1:nf );
            else
                feature_indices = obj.feature_indices;
            end
        end
        
        function [obj, NFxNRxRAD] = getFeasibleNFxNRCombinations( obj, NF, NR )
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % Setup NFMap (speedup boost in radii computations)
            obj = obj.initNFMap(NF);
            
            NFxNR = cartprod( NF, NR );
            NFxNRxRAD = zeros( length( NFxNR ), 3 );
            for i = 1 : length( NFxNR )
               
                nf = NFxNR( i, 1 );
                nr = NFxNR( i, 2 );
                
                [obj, nr_hat, rad] = obj.bisectionForNfNr( nf, nr, 0.1, 1, -1, -1, 0 );
                if ( nr_hat ~= nr )
                    error( "(nf,nr) = (%d,%d) is NOT feasible. nr_hat = %d\n", nf, nr, nr_hat )
                end
                
                NFxNRxRAD( i, : ) = [nf, nr_hat, rad];         
                fprintf( "(nf,nr) = (%d,%d) IS feasible! radii = %f\n", nf, nr, rad )
                    
            end
        end
        
        function obj = initNFMap(obj, NF)
           
            if ( ~isstruct( obj.NFMap ) || numel(fieldnames(( obj.NFMap ))) ~= numel( NF ) )
                obj.NFMap = struct();
                for i = 1 : length( NF )
                    obj.NFMap.(['nf_' num2str(abs(NF(i)))]) = containers.Map('KeyType', 'double', 'ValueType', 'uint32');
                end
            end
            
        end
        
        function [obj, nr_hat, radius] = bisectionForNfNr(obj, nf, nr, radiusLeft, radiusRight, clustersAtLeft, clustersAtRight, bisectionIterCount)
            
            fprintf('radius BT [%f %f]\n', radiusLeft, radiusRight);
            
            nr_hat = -1;
            radius = -1;
            bisectionIterCount = bisectionIterCount + 1;
            if ( bisectionIterCount > obj.BISECTION_MAX_ITERS )
               
                if ( ...
                    obj.ON_BISECTION_MAX_REACHED_POLICY == obj.ON_BISECTION_MAX_REACHED_POLICY_TOLERATE && ...
                    abs( clustersAtLeft - nr ) <= obj.BISECTION_TOLERATE_DISTANCE ...
                )
                   
                    warning( 'BISECTION_MAX_ITERS (=%d) REACHED. TOLERATE Policy has been successfully deployed!\n', obj.BISECTION_MAX_ITERS )
                    nr_hat = nr;
                    radius = radiusLeft;
                    
                end
                
                return;
                
            end
            
            %% Value @ left & right edge
            if ( nargin < 6 || clustersAtLeft <= -1 )
                [obj, clustersAtLeft] = obj.computeNrGivenNfAndRad( nf, radiusLeft );
            end
            
            if ( nargin < 7 || clustersAtRight <= 0 )
                [obj, clustersAtRight] = obj.computeNrGivenNfAndRad( nf, radiusRight );
            end
        
            %% Decide and next step
            if ( nr > clustersAtLeft )
                warning( 'nr > clustersAtLeft (%d > %d)', nr, clustersAtLeft )
                radiusLeftDiv2 = radiusLeft/2;
                [obj, clustersAtLeftDiv2] = obj.computeNrGivenNfAndRad( nf, radiusLeftDiv2 );
                [obj, nr_hat, radius] = obj.bisectionForNfNr( nf, nr, radiusLeftDiv2, radiusLeft, clustersAtLeftDiv2, clustersAtLeft, bisectionIterCount );
            elseif ( nr < clustersAtRight )                
                if abs( clustersAtRight - nr ) <= obj.BISECTION_TOLERATE_DISTANCE
                    warning( 'nr < clustersAtRight (%d < %d)\n\tTOLERATE Policy has been successfully deployed!\n', nr, clustersAtRight )
                    nr_hat = nr;
                    radius = radiusRight;
                else
                    warning( 'nr < clustersAtRight (%d < %d)', nr, clustersAtRight )
                    nr_hat = clustersAtRight;
                    radius = radiusRight;
                end
            elseif ( nr == clustersAtLeft )
                nr_hat = nr;
                radius = radiusLeft;
            elseif ( nr == clustersAtRight )
                nr_hat = nr;
                radius = radiusRight;
            else
                % Get middle point
                %   - Naive: 
                radiusMiddle = ( radiusLeft + radiusRight ) / 2.0;
                %   - Linear Interpolation
                % radiusMiddle = (( clustersAtRight - nr ) * radiusLeft + ( nr - clustersAtLeft ) * radiusRight ) / ( clustersAtRight - clustersAtLeft );
                % Compute for middle point:
                [obj, clustersAtMiddle] = obj.computeNrGivenNfAndRad( nf, radiusMiddle );
                % Decide Direction
                if ( nr > clustersAtMiddle )
                    [obj, nr_hat, radius] = obj.bisectionForNfNr( nf, nr, radiusLeft, radiusMiddle, clustersAtLeft, clustersAtMiddle, bisectionIterCount );
                elseif ( nr < clustersAtMiddle )
                    [obj, nr_hat, radius] = obj.bisectionForNfNr( nf, nr, radiusMiddle, radiusRight, clustersAtMiddle, clustersAtRight, bisectionIterCount );
                else
                    nr_hat = nr;
                    radius = radiusMiddle;
                end
            end 
            
        end
        
        function [obj, clusters] = computeNrGivenNfAndRad( obj, nf, radius )
            
            % Check if combination has been computed
            nf_map_key = ['nf_' num2str(abs(nf))];
            if ( obj.NFMap.(nf_map_key).isKey(radius) )
                clusters = obj.NFMap.(nf_map_key)(radius);
%                 fprintf('\n-->value was previously SAVED! (nf=%02d | radius=%.10f => nr=%02d)\n', nf, radius, clusters)
                return
            end
            
            genOpt = SubtractiveClusteringWrapper.genfisOptionsSC;
            genOpt.ClusterInfluenceRange = radius;
            
            nf_indices = obj.nf2indices(nf);
            clusters = BGenfisNumberOnly( obj.training( :, nf_indices ), obj.training( :, end), genOpt );
            
            obj.NFMap.(nf_map_key)(radius) = clusters;
            fprintf('\tradius = %f --> clusters = %d\n', radius, clusters);
            
        end
        
    end
    
    
    methods (Static)
        
        function genOpt = genfisOptionsSC
            genOpt = genfisOptions( 'SubtractiveClustering' );
%             genOpt.SquashFactor = 1.1;
%             genOpt.RejectRatio = 0.3;
            genOpt.Verbose = false; 
        end
       
        function scw = fromDatasetSplit( dataset, split )
            
            [training, validation, testing] = ...
                AnfisWrapper.partition( unique( dataset, 'rows' ), split );
            clear dataset
            
            scw = SubstractiveClusteringWrapper( training, validation, testing );
            
        end
        
    end
    
end

