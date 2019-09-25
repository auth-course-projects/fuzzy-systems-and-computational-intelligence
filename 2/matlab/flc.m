function delta_theta = flc( dh, dv, theta, which )

    if ( nargin < 4 )
        
        which = 'modified';
        
    end

    persistent fis
    if isempty( fis )
        
        fis = readfis( ['car_flc_' which '.fis'] );
        
    end
    
    [delta_theta, ~, ~, ~, ~] = evalfis( [dh theta dv], fis );

end