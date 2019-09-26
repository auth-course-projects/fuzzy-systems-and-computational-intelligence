classdef SubstractiveClusteringWrapper
    %SUBSTRACTIVECLUSTERINGWRAPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        training
        validation
        testing
        
        feature_indices
        
    end
    
    methods
        function obj = SubstractiveClusteringWrapper( training, validation, testing, getFI )
            %SUBSTRACTIVECLUSTERINGWRAPPER Construct an instance of this class
            %   Detailed explanation goes here
            
            % Save parts
            obj.training = training;
            obj.validation = validation;
            obj.testing = testing;
            
            % Get feature indices ordered by importance (descending order)
            if nargin == 4 && getFI == true
                obj.feature_indices = obj.getFeatureImportanceRL;
            else
                obj.feature_indices = 1 : size( obj.training, 2 ) - 1;
            end
        end
        
        function feature_indices = getFeatureImportanceFE( obj )
           
            ens = fitrensemble( obj.training( :, 1:end-1 ), ...
                obj.training( :, end ), 'Method', 'LSBoost', 'Learners',...
                templateTree( 'MaxNumSplits', 1 ) );
            
            feature_importance = predictorImportance( ens );
            [~, feature_indices] = sort( feature_importance, 2, 'descend' );

        end
        
        function feature_indices = getFeatureImportanceRL( obj )
           
            [ranks, ~] = relieff( obj.training(:, 1:end-1), ...
                obj.training(:, end), 10 );
            [~, feature_indices] = sort( ranks, 2, 'descend' );

        end
        
        function feature_indices = nf2indices( obj, nf )
            %SubstractiveClusteringWrapper.N_FEATURES Uses ReliefF to get
            %$nf most important features of dataset. Returns column indices
            
            if nf > 0
                feature_indices = obj.feature_indices( 1:nf );
            else
                feature_indices = obj.feature_indices;
            end
        end
        
        function NFxNRxRAD = getFeasibleNFxNRCombinations( obj, NF, NR, checkRADSFirst )
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            NFxNR = cartprod( NF, NR );
            NFxNRxRAD = zeros( length( NFxNR ), 3 );
            for i = 1 : length( NFxNR )
               
                nf = NFxNR( i, 1 );
                nr = NFxNR( i, 2 );
                
                if ( nargin == 4 )
                    checkRADFirst = checkRADSFirst( i );
                else
                    checkRADFirst = NaN;
                end
                
                [nr_hat, rad] = obj.getFeasibleNfNrCombination( nf, nr, checkRADFirst );
                if ( nr ~= nr_hat )
                
                    fprintf( "(%d, %d) is NOT feasible. Closest nr = %d\n",...
                        nf, nr, nr_hat )
                    
                    NR( NR == nr ) = nr_hat;
                    NFxNRxRAD = obj.getFeasibleNFxNRCombinations( NF, NR );
                    
                    return
                    
                else
                    
%                     fprintf( "(%d, %d) IS feasible!\n", nf, nr )
                    NFxNRxRAD( i, : ) = [nf, nr_hat, rad];
                    
                end
            end
            
        end
        
        function [nr_hat, rad] = getFeasibleNfNrCombination( obj, nf, nr, checkRADFirst )
            
            rad = 0;
            rad_step = 0.2;
            max_iter = 30;
            genOpt = genfisOptions( 'SubtractiveClustering' );
            genOpt.DataScale = 'auto';
            
            if nargin == 4 && ~isnan( checkRADFirst )
                
                genOpt.ClusterInfluenceRange = checkRADFirst;                
                fis = genfis( obj.training( :, obj.nf2indices(nf)), ...
                    obj.training( :, end), genOpt );
                
                nr_hat = length( fis.rule );
                if nr_hat == nr
                   
                    rad = checkRADFirst;
                    return
                    
                end
                    
            end
            
            while( 1 )
                
                rad = rad + rad_step;
                
                genOpt.ClusterInfluenceRange = rad;                
                fis = genfis( obj.training( :, obj.nf2indices(nf)), ...
                    obj.training( :, end), genOpt );
                
                nr_hat = length( fis.rule );
                max_iter = max_iter - 1;
                if ( nr_hat == nr || 0 == max_iter )
                    
                    break
                    
                elseif ( nr_hat < nr )

                    rad = rad - rad_step;
                    rad_step = rad_step / 2;
                    
                end
                
            end

        end
    end
    
    methods (Static)
       
        function scw = fromDatasetSplit( dataset, split )
            
            [training, validation, testing] = ...
                AnfisWrapper.partition( unique( dataset, 'rows' ), split );
            clear dataset
            
            scw = SubstractiveClusteringWrapper( training, validation, testing );
            
        end
        
    end
    
end

