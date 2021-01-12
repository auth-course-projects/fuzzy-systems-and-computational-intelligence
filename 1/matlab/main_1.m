clear, clc
mi = 19;     % mi = Ki / Kp

%% Design Requirements
Va_max = 200;
Omega_d_max = 150;

%% Linear PI Controller
Kp = Va_max / Omega_d_max;
Ki = mi * Kp;

%% System Definition
% Closed Loop TF (T_L = 0)
H_k = tf( 18.69 * [Kp Ki],[1 12.064+18.69*Kp 18.69*Ki]);
% step(H_k)
linearStepInfo = stepinfo(H_k);

%% Fuzzy PI Controller
% Initial Gains
FZ_alpha_initial = 1 / mi;
FZ_Ke_initial = 1;
FZ_K1_initial = Kp / (FZ_alpha_initial * FZ_Ke_initial);

% Final Gains
FZ_alpha_final = FZ_alpha_initial / 2.5;    % FZ_alpha_initial / 1.75; 
                                            % for 2nd use scenario
FZ_Ke_final = FZ_Ke_initial * 1.5;
FZ_K1_final = FZ_K1_initial * 1.2;

%% Simulation
initialSimOut = sim( 'initial', 'SimulationMode','normal','AbsTol','1e-5',...
    'SaveState','on','StateSaveName','xout',...
    'SaveOutput','on','OutputSaveName','yout',...
    'SaveFormat', 'Dataset'...
);
finalSimOut = sim( 'final', 'SimulationMode','normal','AbsTol','1e-5',...
    'SaveState','on','StateSaveName','xout',...
    'SaveOutput','on','OutputSaveName','yout',...
    'SaveFormat', 'Dataset'...
);

% Results
figure
initialStepInfo = stepinfo( initialSimOut.FPI.signals.values(:,2), ...
    initialSimOut.FPI.time);
plot( initialSimOut.FPI.time, initialSimOut.FPI.signals.values(:,2) );
hold on

finalStepInfo = stepinfo( finalSimOut.FPI.signals.values(:,2), ...
    finalSimOut.FPI.time);
plot( finalSimOut.FPI.time, finalSimOut.FPI.signals.values(:,2) );
hold off

% Compare
%   - rise time
RiseTime = struct('linear', [num2str(1000 * linearStepInfo.RiseTime) ' ms'], ...
    'fz_initial', [num2str(1000 * initialStepInfo.RiseTime) ' ms'], ...
    'fz_final', [num2str(1000 * finalStepInfo.RiseTime) ' ms'] ...
);
%   - settling time
SettlingTime = struct('linear', [num2str(1000 * linearStepInfo.SettlingTime) ' ms'], ...
    'fz_initial', [num2str(1000 * initialStepInfo.SettlingTime) ' ms'], ...
    'fz_final', [num2str(1000 * finalStepInfo.SettlingTime) ' ms'] ...
);
%   - overshoot
Overshoot = struct('linear', [num2str(linearStepInfo.Overshoot) ' %'], ...
    'fz_initial', [num2str(initialStepInfo.Overshoot) ' %'], ...
    'fz_final', [num2str(finalStepInfo.Overshoot) ' %'] ...
);
