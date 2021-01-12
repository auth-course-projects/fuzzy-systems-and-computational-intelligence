clear, clc
drawArrow = @(x, y, varargin) quiver( x(1), y(1), x(2)-x(1), y(2)-y(1), ...
    0, varargin{:} );    % Source: https://stackoverflow.com/a/25730013/4718869
linspace_N = 10;

%% Calculate
x = linspace(0,10, linspace_N);
y = linspace(0,4, linspace_N);
pos_all = cartprod(x,y);
dh = zeros(length(pos_all), 1);
dv = zeros(length(pos_all), 1);

for i = 1 : length(pos_all)
    pos = pos_all(i, 1:2);
    if (...
        (pos(1) >= 5 && pos(2) <=1) ||...
        (pos(1) >= 6 && pos(2) <=2) ||...
        (pos(1) >= 7 && pos(2) <=3)...
    )
        continue;
    end
    
    dh(i) = sensor_h(pos);
    dv(i) = sensor_v(pos);
end

%% Plot
hold on;
for i = 1 : length(pos_all)
    pos = [pos_all(i,1), pos_all(i,2)];
    if (...
        (pos(1) >= 5 && pos(2) <=1) ||...
        (pos(1) >= 6 && pos(2) <=2) ||...
        (pos(1) >= 7 && pos(2) <=3)...
    )
        continue;
    end
    
    plot(pos(1), pos(2),'*b');
    drawArrow([pos_all(i,1), pos_all(i,1) + dh(i)], [pos_all(i,2), pos_all(i,2)], 'color', [0 1 0])
    drawArrow([pos_all(i,1), pos_all(i,1)], [pos_all(i,2), pos_all(i,2) - dv(i)], 'color', [0.1216 0.7255 0.1255])
end

% Draw environment > Borders
line( [5 5], [0 1],  'Color','black', 'LineWidth', 2 ), hold on
line( [6 6], [1 2],  'Color','black', 'LineWidth', 2 )
line( [7 7], [2 3],  'Color','black', 'LineWidth', 2 )
line( [5 6], [1 1],  'Color','black', 'LineWidth', 2 )
line( [6 7], [2 2],  'Color','black', 'LineWidth', 2 )
line( [7 10], [3 3], 'Color','black', 'LineWidth', 2 )
line( [10 10], [0 3], 'Color','black', 'LineWidth', 2 )

% Draw Environment > Face
% r = rectangle( 'Position', [5 0 5 1] );
% set( r, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor', 'none' );
% r = rectangle( 'Position', [6 1 4 1] );
% set( r, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none' );
% r = rectangle( 'Position',[7 2 3 1] );
% set( r, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none' );

hold off