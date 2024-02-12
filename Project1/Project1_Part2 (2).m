% Project 1 - Part 2

% Values To Loop Through
J1 = 100; % Rotational Inertia 1 [kg-m^2]
b1 = 1;% Damping Coefficient 1 [N-m-s/rad]
J2 = 1; % Rotational Inertia 2 [kg-m^2]
b2 = 1; % Damping Coefficient 2 [N-m-s/rad]
T_app_arr = [1,100]; % Constant Applied Torque [N-m]
solver_arr = {'ode1','ode4','ode45'}; % Solvers To Use
dT_arr = [0.1, 1]; % Time Step [s]
% The Following is An Initial Conditon Only For System Option 1
k_arr = [10,100,1000]; % Spring Constant [N/m]
simulation_arr = {'P1_Part2_Option1.slx','P1_Part2_Option2.slx'};

% Preallocate the output cell arrays for efficiency
SimulationNumberData = repmat({''}, 1, 300);
J1Data = repmat({''}, 1, 300);
b1Data = repmat({''}, 1, 300);
J2Data = repmat({''}, 1, 300);
b2Data = repmat({''}, 1, 300);
T_app_arrData = repmat({''}, 1, 300);
solver_arrData = repmat({''}, 1, 300);
dT_arrData = repmat({''}, 1, 300);
k_arrData = repmat({''}, 1, 300);
SimulationData = repmat({''}, 1, 300);
CPUTimeData = repmat({''}, 1, 300);

% Set the titles for the output CSV
SimulationNumberData{1} = 'Simulation #';
J1Data{1} = 'J1 (kg-m^2)';
b1Data{1} = 'b1 (N-m-s/rad)';
J2Data{1} = 'J2 (kg-m^2)';
b2Data{1} = 'b2 (N-m-s/rad)';
T_app_arrData{1} = 'T_app [N-m]';
solver_arrData{1} = 'Solver';
dT_arrData{1} = 'dT (s)';
k_arrData{1} = 'k (N/m)';
SimulationData{1} = 'Simulation';
CPUTimeData{1} = 'CPU Time For Simulation (s)';

% Instantiate counter variables
SimulatedNumber = 0;
IndexVal = 1;
IndexCompile = 0;
IndexCompile1 = 0;

h = waitbar(0, 'Simulation Progress');
for index1 = T_app_arr
    T_app = index1;
    for index2 = 1:length(simulation_arr)
        simulation = simulation_arr{index2};
        for index3 = 1:length(solver_arr)
            solver = solver_arr{index3};
            if strcmp(solver,'ode1') || strcmp(solver,'ode4') % Fixed Time Step
                for index4 = dT_arr % consider timestep for fixed time step case
                    dT = index4;
                    if strcmp(simulation,'P1_Part2_Option1.slx') % if option 1, consider k
                        for index5 = k_arr
                            k = index5;
                            % ---------------------------------------------------------------------------------------- (Fixed Time Step, Option 1)
                            % Tell User Progress on Simulations
                            if IndexCompile == 0
                                clc;
                                disp('Please wait 15 sec while I compile the Simulink model');
                                disp('');
                                disp('Compiling...');
                                simout = sim("P1_Part2_Option1.slx","Solver",solver,"SolverType", 'Fixed-step');
                                IndexCompile= IndexCompile+1;
                            end
                            
                            % Tell User Progress on Simulations
                            clc; 
                            disp(['Please wait 2 min while I compute each simulation case']);
                            disp('');
                            waitbar((SimulatedNumber+1)/40,h,'Simulation Progress');
                            fprintf('Computing Simulation %d/40\r', SimulatedNumber+1);

                            % Simulated Rotational Speed Solution
                            tic;
                            simout = sim(simulation,"Solver",solver,"SolverType", 'Fixed-step');
                            simtime = toc;
                            % Output data
                            if k==10
                                W2a10 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                                Ta10 = simout.tout;
                            elseif k==100
                                W2a100 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                                Ta100 = simout.tout;
                            else %k==1000
                                W2a1000 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                                Ta1000 = simout.tout;
                            end
                              
                            % Increment Indexes
                            SimulatedNumber = SimulatedNumber+1;
                            IndexVal = IndexVal+1;
                            % Update Output With Values
                            SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                            J1Data{IndexVal} = num2str(J1);
                            b1Data{IndexVal} = num2str(b1);
                            J2Data{IndexVal} = num2str(J2);
                            b2Data{IndexVal} = num2str(b2);
                            T_app_arrData{IndexVal} = num2str(T_app);
                            solver_arrData{IndexVal} = num2str(solver);
                            dT_arrData{IndexVal} = num2str(dT);
                            k_arrData{IndexVal} = num2str(k);
                            input_string = simulation;
                            str_length = length(input_string);
                            substring = input_string(str_length - 10 : str_length - 4);
                            SimulationData{IndexVal} = num2str(substring);
                            CPUTimeData{IndexVal} = num2str(simtime);
                        end
                    else % do not consider k for option 2
                        % ---------------------------------------------------------------------------------------- (Fixed Time Step, Option 2)
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 2 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/40,h,'Simulation Progress');
                        fprintf('Computing Simulation %d/40\r', SimulatedNumber+1);

                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim(simulation,"Solver",solver,"SolverType", 'Fixed-step');
                        simtime = toc;
                        % Output data
                        W2b = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                        Tb = simout.tout;
    
                        % Increment Indexes
                        SimulatedNumber = SimulatedNumber+1;
                        IndexVal = IndexVal+1;
                        % Update Output With Values
                        SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                        J1Data{IndexVal} = num2str(J1);
                        b1Data{IndexVal} = num2str(b1);
                        J2Data{IndexVal} = num2str(J2);
                        b2Data{IndexVal} = num2str(b2);
                        T_app_arrData{IndexVal} = num2str(T_app);
                        solver_arrData{IndexVal} = num2str(solver);
                        dT_arrData{IndexVal} = num2str(dT);
                        k_arrData{IndexVal} = 'N/A-No Spring';
                        input_string = simulation;
                        str_length = length(input_string);
                        substring = input_string(str_length - 10 : str_length - 4);
                        SimulationData{IndexVal} = num2str(substring);
                        CPUTimeData{IndexVal} = num2str(simtime);
                    end
                end
            else % Variable Time Step
                if strcmp(simulation,'P1_Part2_Option1.slx') % if option 1, consider k
                    for index5 = k_arr
                        k = index5;
                        % ---------------------------------------------------------------------------------------- (Variable Time Step, Option 1)
                        if IndexCompile1 == 0
                            simout = sim("P1_Part2_Option1.slx","Solver",solver,"SolverType", 'Fixed-step');
                            IndexCompile1= IndexCompile1+1;
                        end
                        
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 2 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/40,h,'Simulation Progress');
                        fprintf('Computing Simulation %d/40\r', SimulatedNumber+1);

                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim(simulation,"Solver",solver,"SolverType", 'Variable-step');
                        simtime = toc;
                        % Output data
                        if k==10
                            W2c10 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                            Tc10 = simout.tout;
                        elseif k==100
                            W2c100 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                            Tc100 = simout.tout;
                        else %k==1000
                            W2c1000 = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                            Tc1000 = simout.tout;
                        end
    
                        % Increment Indexes
                        SimulatedNumber = SimulatedNumber+1;
                        IndexVal = IndexVal+1;
                        % Update Output With Values
                        SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                        J1Data{IndexVal} = num2str(J1);
                        b1Data{IndexVal} = num2str(b1);
                        J2Data{IndexVal} = num2str(J2);
                        b2Data{IndexVal} = num2str(b2);
                        T_app_arrData{IndexVal} = num2str(T_app);
                        solver_arrData{IndexVal} = num2str(solver);
                        dT_arrData{IndexVal} = 'N/A-Variable';
                        k_arrData{IndexVal} = num2str(k);
                        input_string = simulation;
                        str_length = length(input_string);
                        substring = input_string(str_length - 10 : str_length - 4);
                        SimulationData{IndexVal} = num2str(substring);
                        CPUTimeData{IndexVal} = num2str(simtime);
                    end
                else % do not consider k for option 2
                    % ---------------------------------------------------------------------------------------- (Variable Time Step, Option 2)
                    % Tell User Progress on Simulations
                    clc; 
                    disp(['Please wait 2 min while I compute each simulation case']);
                    disp('');
                    waitbar((SimulatedNumber+1)/40,h,'Simulation Progress');
                    fprintf('Computing Simulation %d/40\r', SimulatedNumber+1);

                    % Simulated Rotational Speed Solution
                    tic;
                    simout = sim(simulation,"Solver",solver,"SolverType", 'Variable-step');
                    simtime = toc;
                    % Output data
                    W2d = simout.w2.Data; % Note w2 is shaft speed (rad/s)
                    Td = simout.tout;
    
                    % Increment Indexes
                    SimulatedNumber = SimulatedNumber+1;
                    IndexVal = IndexVal+1;
                    % Update Output With Values
                    SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                    J1Data{IndexVal} = num2str(J1);
                    b1Data{IndexVal} = num2str(b1);
                    J2Data{IndexVal} = num2str(J2);
                    b2Data{IndexVal} = num2str(b2);
                    T_app_arrData{IndexVal} = num2str(T_app);
                    solver_arrData{IndexVal} = num2str(solver);
                    dT_arrData{IndexVal} = 'N/A-Variable';
                    k_arrData{IndexVal} = 'N/A-No Spring';
                    input_string = simulation;
                    str_length = length(input_string);
                    substring = input_string(str_length - 10 : str_length - 4);
                    SimulationData{IndexVal} = num2str(substring);
                    CPUTimeData{IndexVal} = num2str(simtime);
                end
            end
        end
    end
end

close(h);
disp('');
disp('Simulation complete. Outputting data now.');

filename = 'Project1_Part2_OutputData.csv';
dataTable = table(SimulationNumberData',J1Data',b1Data',J2Data',b2Data',T_app_arrData',solver_arrData',solver_arrData',dT_arrData',k_arrData',SimulationData',CPUTimeData');
writetable(dataTable, filename, 'WriteVariableNames', false);
winopen(filename);

% Plot the data
figure;  % Create a new figure
scatter(Ta10, W2a10, 'blue', 'o', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Fixed Option 1, k=10
hold on;
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Ta100, W2a100, 'blue',  's', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Fixed Option 1, k=100
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Ta1000, W2a1000, 'blue',  'd', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Fixed Option 1, k=1000
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Tb, W2b, 'red', 'o', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Fixed Option 2
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Tc10, W2c10, 'green', 'o', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Variable Option 1, k=10
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Tc100, W2c100, 'green', 's', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Variable Option 1, k=100
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Tc1000, W2c1000, 'green', 'd', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Variable Option 1, k=1000
ylim([-2.5, 27.5]);
xlim([0, 25]);
scatter(Td, W2d, 'o', 'magenta', 'filled', 'SizeData', 10, 'MarkerEdgeColor', 'black'); % Variable Option 2
ylim([-2.5, 27.5]);
xlim([0, 25]);
hold off;
xlabel('Time (s)');
ylabel('Shaft Speed (rad/s)');
title('Shaft Speed vs Time (for Various Combination Options and Stiffnesses [k])');
% Add a legend
legend('Fixed, Option 1 (k=10)', 'Fixed, Option 1 (k=100)', 'Fixed, Option 1 (k=1000)', 'Fixed, Option 2', 'Variable, Option 1 (k=10)', 'Variable, Option 1 (k=100)', 'Variable, Option 1 (k=1000)', 'Variable, Option 2');
set(gcf, 'WindowState', 'maximized');

% OLD DATA ---------------------------------------------------------

%{
plot_ts_ode1 = zeros(1, 3);
plot_err_ode1 = zeros(1, 3);
plot_CPU_ode1 = zeros(1, 3);
plot_ts_ode4 = zeros(1, 3);
plot_err_ode4 = zeros(1, 3);
plot_CPU_ode4 = zeros(1, 3);
plot_ts_ode45 = zeros(1, 1);
plot_err_ode45 = zeros(1, 1);
plot_CPU_ode45 = zeros(1, 1);
plot_ts_ode23tb = zeros(1, 1);
plot_err_ode23tb = zeros(1, 1);
plot_CPU_ode23tb = zeros(1, 1);

ode1_0001_idx=0;
ode1_01_idx=0;
ode1_1_idx=0;
ode4_0001_idx=0;
ode4_01_idx=0;
ode4_1_idx=0;
ode45_1_idx=0;
ode23tb_1_idx=0;

indexmain_ode1=1;
indexmain_ode4=1;
indexmain_ode45=1;
indexmain_ode23tb=1;
for index8 = 1:length(solver_Data)
    if (solver_Data(index8) == "ode1") && (dT_Data(index8) == "0.001") && (ode1_0001_idx==0)
        plot_ts_ode1(indexmain_ode1)=str2double(dT_Data{index8});
        plot_err_ode1(indexmain_ode1)=str2double(MaxErrorData{index8});
        plot_CPU_ode1(indexmain_ode1)=str2double(CPUTimeData{index8});
        indexmain_ode1=indexmain_ode1+1;
        ode1_0001_idx=ode1_0001_idx+1;
    end
    if (solver_Data(index8) == "ode1") && (dT_Data(index8) == "0.1") && (ode1_01_idx==0)
        plot_ts_ode1(indexmain_ode1)=str2double(dT_Data{index8});
        plot_err_ode1(indexmain_ode1)=str2double(MaxErrorData{index8});
        plot_CPU_ode1(indexmain_ode1)=str2double(CPUTimeData{index8});
        indexmain_ode1=indexmain_ode1+1;
        ode1_01_idx=ode1_01_idx+1;
    end
    if (solver_Data(index8) == "ode1") && (dT_Data(index8) == "1") && (ode1_1_idx==0)
        plot_ts_ode1(indexmain_ode1)=str2double(dT_Data{index8});
        plot_err_ode1(indexmain_ode1)=str2double(MaxErrorData{index8});
        plot_CPU_ode1(indexmain_ode1)=str2double(CPUTimeData{index8});
        indexmain_ode1=indexmain_ode1+1;
        ode1_1_idx=ode1_1_idx+1;
    end
    if (solver_Data(index8) == "ode4") && (dT_Data(index8) == "0.001") && (ode4_0001_idx==0)
        plot_ts_ode4(indexmain_ode4)=str2double(dT_Data{index8});
        plot_err_ode4(indexmain_ode4)=str2double(MaxErrorData{index8});
        plot_CPU_ode4(indexmain_ode4)=str2double(CPUTimeData{index8});
        indexmain_ode4=indexmain_ode4+1;
        ode4_0001_idx=ode4_0001_idx+1;
    end
    if (solver_Data(index8) == "ode4") && (dT_Data(index8) == "0.1") && (ode4_01_idx==0)
        plot_ts_ode4(indexmain_ode4)=str2double(dT_Data{index8});
        plot_err_ode4(indexmain_ode4)=str2double(MaxErrorData{index8});
        plot_CPU_ode4(indexmain_ode4)=str2double(CPUTimeData{index8});
        indexmain_ode4=indexmain_ode4+1;
        ode4_01_idx=ode4_01_idx+1;
    end
    if (solver_Data(index8) == "ode4") && (dT_Data(index8) == "1") && (ode4_1_idx==0)
        plot_ts_ode4(indexmain_ode4)=str2double(dT_Data{index8});
        plot_err_ode4(indexmain_ode4)=str2double(MaxErrorData{index8});
        plot_CPU_ode4(indexmain_ode4)=str2double(CPUTimeData{index8});
        indexmain_ode4=indexmain_ode4+1;
        ode4_1_idx=ode4_1_idx+1;
    end
    if (solver_Data(index8) == "ode45") && (dT_Data(index8) == "1") && (ode45_1_idx==0)
        plot_ts_ode45(indexmain_ode45)=str2double(dT_Data{index8});
        plot_err_ode45(indexmain_ode45)=str2double(MaxErrorData{index8});
        plot_CPU_ode45(indexmain_ode45)=str2double(CPUTimeData{index8});
        indexmain_ode45=indexmain_ode45+1;
        ode45_1_idx=ode45_1_idx+1;
    end
    if (solver_Data(index8) == "ode23tb") && (dT_Data(index8) == "1") && (ode23tb_1_idx==0)
        plot_ts_ode23tb(indexmain_ode23tb)=str2double(dT_Data{index8});
        plot_err_ode23tb(indexmain_ode23tb)=str2double(MaxErrorData{index8});
        plot_CPU_ode23tb(indexmain_ode23tb)=str2double(CPUTimeData{index8});
        indexmain_ode23tb=indexmain_ode23tb+1;
        ode23tb_1_idx=ode23tb_1_idx+1;
    end 
end

% Plot the data
figure;  % Create a new figure
scatter(plot_ts_ode1, plot_err_ode1, 'o', 'filled', 'MarkerEdgeColor', 'b','SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(plot_ts_ode4, plot_err_ode4, 'o', 'filled', 'MarkerEdgeColor', 'r','SizeData', 10);
scatter(plot_ts_ode45, plot_err_ode45, 'o', 'filled', 'MarkerEdgeColor', 'y','SizeData', 10);
scatter(plot_ts_ode23tb, plot_err_ode23tb, 'o', 'filled', 'MarkerEdgeColor', 'g','SizeData', 10);
hold off;  % Release the current plot hold
% Add labels and title
xlabel('Time Step (s)');
ylabel('Max Simulation Error (w)');
title('Max Simulation Error vs Time Step');
% Add a legend
legend('ode1', 'ode4','ode45','ode23tb');

% Plot the data
figure;  % Create a new figure
scatter(plot_ts_ode1, plot_CPU_ode1, 'o', 'filled', 'MarkerEdgeColor', 'b','SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(plot_ts_ode4, plot_CPU_ode4, 'o', 'filled', 'MarkerEdgeColor', 'r','SizeData', 10);
scatter(plot_ts_ode45, plot_CPU_ode45, 'o', 'filled', 'MarkerEdgeColor', 'y','SizeData', 10);
scatter(plot_ts_ode23tb, plot_CPU_ode23tb, 'o', 'filled', 'MarkerEdgeColor', 'g','SizeData', 10);
hold off;  % Release the current plot hold
% Add labels and title
xlabel('Time Step (s)');
ylabel('CPU Time Taken (s)');
title('CPU Time Taken vs Time Step');
% Add a legend
legend('ode1', 'ode4','ode45','ode23tb');

% Plot the data
figure;  % Create a new figure
scatter(plot_CPU_ode1, plot_err_ode1, 'o', 'filled', 'MarkerEdgeColor', 'b','SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(plot_CPU_ode4, plot_err_ode4, 'o', 'filled', 'MarkerEdgeColor', 'r','SizeData', 10);
scatter(plot_CPU_ode45, plot_err_ode45, 'o', 'filled', 'MarkerEdgeColor', 'y','SizeData', 10);
scatter(plot_CPU_ode23tb, plot_err_ode23tb, 'o', 'filled', 'MarkerEdgeColor', 'g','SizeData', 10);
hold off;  % Release the current plot hold
% Add labels and title
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
title('Max Simulation Error vs CPU Time Taken');
% Add a legend
legend('ode1', 'ode4','ode45','ode23tb');




%dT = [0.001, 0.1, 1]; % Time Step [s]
%solver = ["ode1", "ode4"]; % Fixed Time Step Solver [Euler]
%simout = sim("P1_demo.slx","Solver",solver,"FixedStep",string(dT));
%}

