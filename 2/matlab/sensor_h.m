function dh = sensor_h( pos )

    if pos( 2 ) <= 1
        
        dh = 5 - pos( 1 );
        
    elseif pos( 2 ) <= 2
        
        dh = 6 - pos( 1 );
        
    elseif pos( 2 ) <= 3
        
        dh = 7 - pos( 1 );
        
    else
        
        dh = 1;     % max
        
    end
       
end

