clear, clear flc; clc
drawArrow = @(x, y, varargin) quiver( x(1), y(1), x(2)-x(1), y(2)-y(1), ...
    0, varargin{:} );    % Source: https://stackoverflow.com/a/25730013/4718869

%% Data
target = [10; 3.2];
speed = 0.05;   % m / sec

% Initial Pose
pos_initial = [3.8; 0.5];
theta_initial = 0;

which_fis = 'modified';  % initial_no_dv, initial, modified

%% Control Loop
route = zeros( 200, 2 );
route_i = 1;

pos = pos_initial;
theta = theta_initial;

disp( "Initial Target: (x, y) = (" + num2str( target(1) ) + ", " + ...
            num2str( target(2) ) + ")" );

% FIX: Track step when one coordinate becomes the same as target point's
% To track step, will retain previous point with the same property
ref_point = NaN;
conv_step_n = 0;
n = 1;
        
while( 1 )
   
    % Get Sensor Measurements
    dh = sensor_h( pos );
    dv = sensor_v( pos );
    
    % Run FLC
    delta_theta = flc( dh, dv, theta, which_fis );
    
    % Apply
    theta = theta + delta_theta;
    
    % Calc next position
    if (route_i > 0)
        pos = pos + speed * [cosd(theta); sind(theta)];
    end
    route( route_i, : ) = pos';
    
    % FIX: Check if any of the coords became equal
    % For given problem, check if y coord has been approx. to millimeter
    if ( abs( pos( 2 ) - target( 2 ) ) < 1e-3 )
        
        
        if ( isnan( ref_point ) )
            ref_point = pos( 1 );
        else
            % Get convergence step
            current_step = abs( pos(1) - ref_point );
            conv_step_n = conv_step_n + ( current_step - conv_step_n ) / n;
            
            % Predict on given step
            % If going the next step will increase the error from target
            % then, stop on this step
            if ( abs( target(1) - ( pos(1) + conv_step_n ) ) > ...
                    abs( target(1) - pos(1) ) )
                disp( "Final Position: (x, y) = (" + num2str( pos(1) ) + ", " + ...
                    num2str( pos(2) ) + ")" );
                break;
            end
            
            % Increment for next step
            ref_point = pos(1);
            n = n + 1;
        end
        
    end
    
    % If pos_y > target_y: STOP
    if ( pos( 1 ) >= target( 1 ) )
        
        disp( "Final Position: (x, y) = (" + num2str( pos(1) ) + ", " + ...
            num2str( pos(2) ) + ")" );
        break;
        
    else
        
        route_i = route_i + 1;
    
    end
    
end

%% Results
distance = sqrt( sumsqr( [ pos(1) - target(1); pos(2) - target(2) ] ) );
disp( "Euclidean Distance: " + num2str( 1000 * distance ) + " mm" )

% Draw environment > Borders
figure
line( [5 5], [0 1],  'Color','black', 'LineWidth', 2 ), hold on
line( [6 6], [1 2],  'Color','black', 'LineWidth', 2 )
line( [7 7], [2 3],  'Color','black', 'LineWidth', 2 )
line( [5 6], [1 1],  'Color','black', 'LineWidth', 2 )
line( [6 7], [2 2],  'Color','black', 'LineWidth', 2 )
line( [7 10], [3 3], 'Color','black', 'LineWidth', 2 )
line( [10 10], [0 3], 'Color','black', 'LineWidth', 2 )

% Draw Environment > Face
r = rectangle( 'Position', [5 0 5 1] );
set( r, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor', 'none' );
r = rectangle( 'Position', [6 1 4 1] );
set( r, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none' );
r = rectangle( 'Position',[7 2 3 1] );
set( r, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none' );

% Draw initial pose
plot( pos_initial(1), pos_initial(2), 'og' )

arrow_end = pos_initial' + 10 * speed * [cosd( theta_initial ) sind( theta_initial )];
drawArrow( [pos_initial(1) arrow_end(1)], [pos_initial(2) arrow_end(2)], 'linewidth', 2, 'color', 'g' );

% Draw target
plot( target(1), target(2), '+r', 'MarkerSize', 20 )

% Draw result route
plot( route(1:route_i - 1, 1), route(1:route_i - 1, 2), 'ob' )
% Get result route length & Draw Interpolating Curve
rlen = route_length(route, route_i, true);
disp( "Route Distance: " + num2str( rlen ) + " m" )

% Draw final pose
plot( pos(1), pos(2), 'og' )

arrow_end = pos' + 10 * speed * [cosd( theta ) sind( theta )];
drawArrow( [pos(1) arrow_end(1)], [pos(2) arrow_end(2)], 'linewidth', 2, 'color', 'g' );
hold off;
