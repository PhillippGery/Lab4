%[text] Step 1: Create Trapezoidal Velocity Trajectories
% clear all; % Auskommentiert, damit 'final_scaled_paths' nicht gelöscht wird
close all;
clc;

% --- 1. Pfad-Daten laden ---
% ANNAHME: Die Variable 'final_scaled_paths' existiert aus deinem Logo-Skript
if ~exist('final_scaled_paths', 'var')
    try
        load('PG_daten.mat'); % Versucht, die Daten aus einer .mat-Datei zu laden
    catch
        error('Die Pfad-Variable ''final_scaled_paths'' wurde nicht gefunden. Bitte lade zuerst dein Logo.');
    end
end

% --- 2. Parameter definieren ---
v_max = 0.25;    % Maximale Geschwindigkeit (m/s)
t_accel = 1.0;  % Beschleunigungs-/Verzögerungszeit (in Sekunden).
dt = 0.002;   % Zeitschritt (s)
num_transition_points = 100; % Punkte für An-/Abfahrt (für glatte Interpolation)

% NEU: 3D-Offset für die gesamte Zeichnung (X, Y, Z) in Metern
drawing_offset = [0, -0.2, -0.2]; % z.B. 10cm in X, 0cm in Y, 2cm in Z
fprintf('Zeichnung wird auf (%.3f, %.3f, %.3f) verschoben.\n', ...
         drawing_offset(1), drawing_offset(2), drawing_offset(3));

% --- 3. Master-Pfad & Distanzen für ALLE Segmente berechnen ---
p_origin = [0, 0, 0];    % Start- und Endpunkt ist 3D-Ursprung
current_pos = p_origin;   % Startet am Ursprung

G_master_cell = {};       % Cell array, um alle Pfad-Matrizen zu sammeln
led_distance_markers = [];% Array, um die Länge jedes Segments zu speichern
num_paths = length(final_scaled_paths);

fprintf('Starte 3D-Pfad-Zusammenführung für %d Logo-Segmente...\n', num_paths);

for k = 1:num_paths
    % Holt den nächsten 2D-Logo-Pfad
    path_k_2D = final_scaled_paths{k};
    
    % --- 3D-Verschiebung des Pfades ---
    x_shifted = path_k_2D(:, 1) + drawing_offset(1);
    y_shifted = path_k_2D(:, 2) + drawing_offset(2);
    z_shifted = ones(size(path_k_2D, 1), 1) * drawing_offset(3);
    
    % path_k ist jetzt ein 3D-Pfad
    path_k = [x_shifted, y_shifted, z_shifted];
    
    p_start = path_k(1, :);  % 3D-Startpunkt des Logo-Segments k
    p_end = path_k(end, :);    % 3D-Endpunkt des Logo-Segments k
    
    % --- Segment A (Transition, LED AUS) ---
    % Fährt von 'current_pos' zum 3D-Startpunkt 'p_start'
    path_A_x = linspace(current_pos(1), p_start(1), num_transition_points)';
    path_A_y = linspace(current_pos(2), p_start(2), num_transition_points)';
    path_A_z = linspace(current_pos(3), p_start(3), num_transition_points)'; % 3D
    path_A = [path_A_x, path_A_y, path_A_z];
    
    % Zum Master-Pfad hinzufügen
    if k == 1
        G_master_cell{end+1} = path_A; % Erste Anfahrt vom Ursprung
    else
        G_master_cell{end+1} = path_A(2:end, :); % Übergang zwischen Segmenten
    end
    
    % 3D-Distanz berechnen
    d_A = norm(p_start - current_pos);
    led_distance_markers(end+1) = d_A;
    
    % --- Segment B (Logo-Pfad, LED AN) ---
    G_master_cell{end+1} = path_k(2:end, :); % (vermeidet doppelten Startpunkt)
    
    % 3D-Distanz des Logo-Pfades berechnen
    dx_B = diff(path_k(:,1)); 
    dy_B = diff(path_k(:,2)); 
    dz_B = diff(path_k(:,3)); % 3D
    d_B = sum(sqrt(dx_B.^2 + dy_B.^2 + dz_B.^2));
    led_distance_markers(end+1) = d_B;
    
    % Aktualisiere die 'current_pos' für die nächste Iteration
    current_pos = p_end;
    fprintf('Pfad %d (Länge %.3fm) und Übergang (Länge %.3fm) hinzugefügt.\n', k, d_B, d_A);
end

% --- Segment C (Finale Rückfahrt, LED AUS) ---
% Fährt vom letzten 3D-Punkt (current_pos) zurück zum 3D-Ursprung
p_target = p_origin; 
path_C_x = linspace(current_pos(1), p_target(1), num_transition_points)';
path_C_y = linspace(current_pos(2), p_target(2), num_transition_points)';
path_C_z = linspace(current_pos(3), p_target(3), num_transition_points)'; % 3D
path_C = [path_C_x, path_C_y, path_C_z];

G_master_cell{end+1} = path_C(2:end, :);
d_C = norm(p_target - current_pos); % 3D-Distanz
led_distance_markers(end+1) = d_C;
fprintf('Rückfahrt zum Ursprung (Länge %.3fm) hinzugefügt.\n', d_C);

% --- 4. Finalen Master-Pfad und Gesamtdistanz erstellen ---
G_master = vertcat(G_master_cell{:});
x_master = G_master(:, 1);
y_master = G_master(:, 2);
z_master = G_master(:, 3); % 3D

d = sum(led_distance_markers); % Gesamtdistanz der gesamten Operation
fprintf('GESAMTDISTANZ (3D): %.4f m\n', d);

% Kumulativen 3D-Distanzvektor für G_master erstellen
dx_master = diff(x_master); 
dy_master = diff(y_master); 
dz_master = diff(z_master); % 3D
s_master = [0; cumsum(sqrt(dx_master.^2 + dy_master.^2 + dz_master.^2))];

% --- 5. Trapez-/Dreieck-Profil für die GESAMTE Distanz 'd' ---
% (Dieser Block ist identisch)
d_ramp_full = v_max * t_accel;
if d >= d_ramp_full
    d_const = d - d_ramp_full; t_const = d_const / v_max;
    tfinal = 2 * t_accel + t_const; v_profile_max = v_max;
    fprintf('Profil: Trapez | v_max: %.2f m/s | t_final: %.2f s\n', v_profile_max, tfinal);
else
    t_const = 0; tfinal = 2 * t_accel; v_profile_max = d / t_accel;
    fprintf('Profil: Dreieck | v_max (reduziert): %.2f m/s | t_final: %.2f s\n', v_profile_max, tfinal);
end

% --- 6. Geschwindigkeitsprofil v_desired(t) erstellen ---
% (Dieser Block ist identisch)
t = (0:dt:tfinal)'; 
N = length(t);
v_desired = zeros(N, 1);
for i = 1:N
    if t(i) <= t_accel
        v_desired(i) = (v_profile_max / t_accel) * t(i);
    elseif t(i) > t_accel && t(i) <= t_accel + t_const
        v_desired(i) = v_profile_max;
    else
        t_remaining = tfinal - t(i);
        v_desired(i) = (v_profile_max / t_accel) * t_remaining;
    end
end
v_desired(v_desired < 0) = 0; v_desired(end) = 0;

% --- 7. Pfad neu "timen" durch 3D-Interpolation ---
d_desired = cumtrapz(t, v_desired);
[s_master_unique, ia] = unique(s_master);
x_master_unique = x_master(ia);
y_master_unique = y_master(ia);
z_master_unique = z_master(ia); % 3D

% Interpolieren, um die finalen Trajektorien x(t), y(t) und z(t) zu erhalten
x = interp1(s_master_unique, x_master_unique, d_desired, 'linear', 'extrap');
y = interp1(s_master_unique, y_master_unique, d_desired, 'linear', 'extrap');
z = interp1(s_master_unique, z_master_unique, d_desired, 'linear', 'extrap'); % 3D

% Sicherstellen, dass der letzte Punkt exakt der 3D-Ursprung ist
x(end) = p_origin(1);
y(end) = p_origin(2);
z(end) = p_origin(3);

% --- 8. LED-Status-Vektor erstellen ---
% (Dieser Block ist identisch)
led_status = zeros(N, 1);
cumulative_dist_markers = cumsum(led_distance_markers);
current_dist_marker = 0;
for k = 1:length(led_distance_markers)
    is_led_on = (mod(k, 2) == 0); % Segmente 2, 4, 6... sind Logo-Pfade
    if is_led_on
        indices = (d_desired > current_dist_marker) & ...
                  (d_desired <= cumulative_dist_markers(k));
        led_status(indices) = 1;
    end
    current_dist_marker = cumulative_dist_markers(k);
end
led_status(end) = 0;

% --- 9. 3D-Geschwindigkeitsprofil plotten ---
v_actual = sqrt( (diff(x)/dt).^2 + (diff(y)/dt).^2 + (diff(z)/dt).^2 ); % 3D
figure;
plot(t(1:end-1), v_actual, 'k', 'MarkerSize', 4);
hold on;
plot(t, v_desired, 'g--', 'LineWidth', 1.5);
yline(v_profile_max, 'r--', 'LineWidth', 1.5);
title('3D-Geschwindigkeitsprofil (Gesamte Trajektorie)');
grid on; xlabel('time (s)'); ylabel('velocity (m/s)');
legend('Tatsächliche Geschwindigkeit (v)', 'Erwünschte (v\_desired)', 'Max Geschw.');
ylim([0 v_max * 1.1]);

% --- 10. 3D-Trajektorien-Plot ---
figure;
hold on;
% Teile des Pfades basierend auf LED-Status finden
path_off = [x(led_status==0), y(led_status==0), z(led_status==0)];
path_on  = [x(led_status==1), y(led_status==1), z(led_status==1)];
% Plotten
plot3(path_off(:,1), path_off(:,2), path_off(:,3), 'k.', 'MarkerSize', 2); % An/Abfahrt
plot3(path_on(:,1), path_on(:,2), path_on(:,3), 'b.', 'MarkerSize', 4);  % Logo
% Start- und Endpunkt hervorheben
plot3(p_origin(1), p_origin(2), p_origin(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);
axis equal; grid on;
title('Gefahrene 3D-Trajektorie (inkl. An- und Abfahrt)');
xlabel('X (Meter)'); ylabel('Y (Meter)'); zlabel('Z (Meter)');
legend('An/Abfahrt (LED AUS)', 'Logo-Pfad (LED AN)', 'Start/Ende (0,0,0)');
view(30, 20); % Guter 3D-Blickwinkel


%%
%[text] ## Step 2: Forward Kinematics
%[text] (2c) Calculate T0
% these values were obtained from the URDF directly
L1 = 0.2435;
L2 = 0.2132;
W1 = 0.1311;
W2 = 0.0921;
H1 = 0.1519;
H2 = 0.0854;

% home position of end effector
M = [-1 0 0 L1+L2;
    0 0 1 W1+W2;
    0 1 0 H1-H2;
    0 0 0 1];

% screw axes
S1 = [0 0 1 0 0 0]';
S2 = [0 1 0 -H1 0 0]';
S3 = [0 1 0 -H1 0 L1]';
S4 = [0 1 0 -H1 0 L1+L2]';
S5 = [0 0 -1 -W1 L1+L2 0]';
S6 = [0 1 0 H2-H1 0 L1+L2]';
S = [S1 S2 S3 S4 S5 S6];

% body screw axes
B1 = ECE569_Adjoint(M)\S1;
B2 = ECE569_Adjoint(M)\S2;
B3 = ECE569_Adjoint(M)\S3;
B4 = ECE569_Adjoint(M)\S4;
B5 = ECE569_Adjoint(M)\S5;
B6 = ECE569_Adjoint(M)\S6;
B = [B1 B2 B3 B4 B5 B6];

% joint angles
theta0 = [-1.6800   -1.4018   -1.8127   -2.9937   -0.8857   -0.0696]';
%% 

% calculate the 4x4 matrix representing the transition
% from end effector frame {b} to the base frame {s} at t=0: Tsb(0)

% TODO: implement ECE569_FKinSpace and ECE569_FKinBody
T0_space = ECE569_FKinSpace(M,S,theta0);
T0_body = ECE569_FKinBody(M,B,theta0);
T0_space-T0_body;
T0 = T0_body;
%[text] Calculate Tsd at every time step.
% Calculate Tsd(t) for t=0 to t=tfinal
% Tsd(t) = T0 * Td(t)
N = length(t);
Tsd = zeros(4,4,N);

for i=1:N
    pd = [x(i); y(i); 0];
    R = eye(3);
    Td = [R, pd; 0, 0, 0, 1]; 
    Tsd(:,:,i) = T0*Td;
end

%%
%[text] (2d) Plot (x,y,z) in the s frame
xs = Tsd(1,4,:);
ys = Tsd(2,4,:);
zs = Tsd(3,4,:);
plot3(xs(:), ys(:), zs(:), 'LineWidth', 1)
title('Trajectory \{s\} frame')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
hold on
plot3(xs(1),ys(1),zs(1),'go','MarkerSize',10,'LineWidth',2)
plot3(xs(end),ys(end),zs(end),'rx','MarkerSize',10,'LineWidth',2)
legend('Trajectory', 'Start', 'End')
grid on
hold off
%%
%[text] ## Step 3: Inverse Kinematics
initialguess = theta0;
Td = T0;

% you need to implement IKinBody
[thetaSol, success] = ECE569_IKinBody(B,M,Td,theta0,1e-6,1e-6);
if (~success)
    close(f);
    error('Error. \nCould not perform IK at index %d.',1)
end
%%
%[text] (3c) Perform IK at each time step
thetaAll = zeros(6,N);
thetaAll(:,1) = theta0;

% you can comment out the waitbar functions if they aren't working
% (sometimes they don't work with .mlx files)
% If the code gets stuck here, you will need to restart MATLAB
f = waitbar(0,['Inverse Kinematics (1/',num2str(N),') complete.']);



for i=2:N
    % TODO: use previous solution as current guess
    initialguess = thetaAll(:,i-1);

    % TODO: calculate thetaSol for Tsd(:,:,i) with initial guess
    [thetaSol, success] = ECE569_IKinBody(B,M,Tsd(:,:,i),initialguess,1e-6,1e-6);
    if (~success)
        close(f);
        error('Error. \nCould not perform IK at index %d.',i)
    end
    thetaAll(:,i) = thetaSol;
    waitbar(i/N,f,['Inverse Kinematics (',num2str(i),'/',num2str(N),') complete.']);
end
close(f);
%%
%[text] (3c) Verify that the joint angles don't change very much
dj = diff(thetaAll');
plot(t(1:end-1), dj)
title('First Order Difference in Joint Angles')
legend('J1','J2','J3','J4','J5','J6','Location','northeastoutside')
grid on
xlabel('time (s)')
ylabel('first order difference')
%%
%[text] (3d) Verify that the joints we found actually trace out our trajectory (forward kinematics)
actualTsd = zeros(4,4,N);

for i=1:N
    % TODO: use forward kinematics to calculate Tsd from our thetaAll
    actualTsd(:,:,i) = ECE569_FKinBody(M, B, thetaAll(:,i));
end

xs = actualTsd(1,4,:);
ys = actualTsd(2,4,:);
zs = actualTsd(3,4,:);
plot3(xs(:), ys(:), zs(:), 'LineWidth', 1)
title('Verified Trajectory \{s\} frame')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
hold on
plot3(xs(1),ys(1),zs(1),'go','MarkerSize',10,'LineWidth',2)
plot3(xs(end),ys(end),zs(end),'rx','MarkerSize',10,'LineWidth',2)
legend('Trajectory', 'Start', 'End')
grid on
hold off

%%
%[text] (3e) Verify that the end effector does not enter a kinematic singularity, by plotting the determinant of your body jacobian
body_dets = zeros(N,1);
for i=1:N
    Jb = ECE569_JacobianBody(B, thetaAll(:,i)); 
    body_dets(i) = det(Jb);
end
plot(t, body_dets)
title('Manipulability')
grid on
xlabel('time (s)')
ylabel('det of J_B')

%%
led = led_status; % <--- NEUE ZEILE
% save to the CSV file
data = [t thetaAll' led];

writematrix(data, 'Sign_06.csv')