function rlen = route_length(route, route_i_end, print_interp)
    % Fit Curcve to Data
    pchip_struct = pchip( route(1:route_i_end-1, 1), route(1:route_i_end-1, 2) );
    
    % Print Interpolation Curve
    if ( print_interp )
        plot( route(1:route_i_end-1, 1), ppval( pchip_struct, route(1:route_i_end-1, 1) ), '-r' );
    end

    % Caclulate Distance
    rlen = integral( @(xx) ppval( pchip_struct, xx), route(1, 1), route(route_i_end-1, 1) );
end