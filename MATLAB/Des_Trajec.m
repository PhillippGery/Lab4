% --- 1. Define Parameters ---
% These values are taken directly from your Desmos screenshot
A = 0.1485;
B = 0.1136;
a = 2;
b = 3;

% --- 2. Define Time Vector ---
% The trajectory runs from t = 0 to T = 2*pi
T_end = 2 * pi;
t = linspace(0, T_end, 1000); % 1000 points for a smooth curve

% --- 3. Calculate Trajectory ---
x = A * sin(a * t);
y = B * sin(b * t);

% --- 4. Create the Plot ---
figure;
hold on; % Keep the plot active to add the boundary box

% Plot the Lissajous curve
plot(x, y, 'b', 'LineWidth', 2);

% --- 5. Add Workspace Boundary ---
% The problem specifies a workspace boundary of |x| <= 0.16 and |y| <= 0.16
workspace_lim = 0.16;

% Draw the boundary box as a dashed red rectangle
rectangle('Position', [-workspace_lim, -workspace_lim, 2*workspace_lim, 2*workspace_lim], ...
          'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% --- 6. Format the Plot ---
title('Lissajous Curve Trajectory');
xlabel('x_d (m)');
ylabel('y_d (m)');
grid on;
axis equal; % Ensures x and y axes have the same scale

% Set axis limits to be slightly larger than the workspace
axis_padding = 0.02;
axis_lim = workspace_lim + axis_padding;
axis([-axis_lim, axis_lim, -axis_lim, axis_lim]);

% Add a legend
legend('Trajectory', 'Workspace Boundary', 'Location', 'northwest');

hold off; % Release the plot
 
