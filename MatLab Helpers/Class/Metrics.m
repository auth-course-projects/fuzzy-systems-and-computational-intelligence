classdef Metrics
    
    properties
        
        N
        
        MSE
        RMSE
        NMSE
        R2
        NDEI
        
    end
    
    methods
        
        function obj = Metrics(output, true_output)
            
            obj.N = length( output );
            
            %   - MSE
            obj.MSE = sumsqr( true_output - output ) / obj.N;
            %   - RMSE
            obj.RMSE = sqrt( obj.MSE );
            %   - NMSE
            obj.NMSE = ( obj.MSE * obj.N ) / ...
                sumsqr( true_output - mean( true_output ) );
            %   - R^2
            obj.R2 = 1 - obj.NMSE;
            %   - NDEI
            obj.NDEI = sqrt( obj.NMSE );
            
        end
        
    end
end

