classdef Metrics
    
    properties
        
        N
        classification = false
        
        % Regression Metrics
        MSE
        RMSE
        NMSE
        R2
        NDEI
        
        % Classification Metrics
        OA
        PU
        UA
        K_HAT
        
    end
    
    methods
        
        function obj = Metrics( output, true_output, classification )
            
            obj.N = length( output );
            if obj.N == 0
                return
            end
            
            if nargin == 3 && classification
                obj.classification = classification;
                [obj.OA, obj.PU, obj.UA, obj.K_HAT] = ...
                    Metrics.classification_metrics( output, true_output );
            else
                [obj.MSE, obj.RMSE, obj.NMSE, obj.R2, obj.NDEI] = ...
                    Metrics.regression_metrics( output, true_output );
            end
            
        end
        
    end
    
    methods ( Static )
    
        function [mse, rmse, nmse, r2, ndei] = regression_metrics( output, true_output )
           
            %   - MSE
            mse = sumsqr( true_output - output ) / length( output );
            %   - RMSE
            rmse = sqrt( mse );
            %   - NMSE
            nmse = ( mse * length( output ) ) / ...
                sumsqr( true_output - mean( true_output ) );
            %   - R^2
            r2 = 1 - nmse;
            %   - NDEI
            ndei = sqrt( nmse );
            
        end
        
        function [OA, PU, UA, K_HAT] = classification_metrics( true_output, output )
           
            % Normalize
            classes_range = [min(true_output), max(true_output)];
            output( output < classes_range( 1 ) ) = classes_range( 1 );
            output( output > classes_range( 2 ) ) = classes_range( 2 );
            
            % Get confusion matrix
            C = confusionmat( true_output, output );
            N = length( output );
            
            x_ir = sum( C )';
            x_jc = sum( C, 2 );
            diag_C = diag( C );
                      
            %   - OVERALL ACCURACY
            OA = sum( diag_C ) / N;
            %   - PRODUCER ACCURACY            
            PU = diag_C ./ x_jc;
            %   - USER ACCURACY
            UA = diag_C ./ x_ir;
            %   - K_HAT
            p_o = OA;
            p_e = ( 1/N^2 ) * sum( x_ir .* x_jc );
            K_HAT = (p_o - p_e) / (1 - p_e);
            
        end
        
    end
    
end

