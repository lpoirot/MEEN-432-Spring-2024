% Project 1 Demo/Scaffolding
% this code is to be used a start-off point, do not expect this code to be
% perfect already.

% Initial Conditions To Loop Through
w_0_arr = [10,0]; % Initial Angular Velocity [rad/s]
J_arr = [100,0.01]; % Rotational Inertia [kg-m^2]
b_arr = [10,0.1]; % Damping Coefficient [N-m-s/rad]
A_arr = [0,100]; % Constant Applied Torque [N-m]
solver_arr = {'ode1','ode4','ode45','ode23tb'};
dT_arr = [0.001, 0.1, 1]; % Time Step [s]
f_arr = [0.1,100]; % Sinusoidal Toque - Frequency
torquetype = {'Constant','Sinusoidal'};
phi = [0,deg2rad(90)];

% Instantiate variables
SimulatedNumber = 0;
IndexVal = 1;

% Preallocate the output cell arrays for efficiency
SimulationNumberData = repmat({''}, 1, 300);
w_0Data = repmat({''}, 1, 300);
J_Data = repmat({''}, 1, 300);
b_Data = repmat({''}, 1, 300);
solver_Data = repmat({''}, 1, 300);
dT_Data = repmat({''}, 1, 300);
A_Data = repmat({''}, 1, 300);
f_Data = repmat({''}, 1, 300);
torquetype_Data = repmat({''}, 1, 300);
CPUTimeData = repmat({''}, 1, 300);
MaxErrorData = repmat({''}, 1, 300);
% Set the titles for the output CSV
SimulationNumberData{1} = 'Simulation #';
w_0Data{1} = 'w0';
J_Data{1} = 'J';
b_Data{1} = 'b';
solver_Data{1} = 'Solver';
dT_Data{1} = 'dT';
A_Data{1} = 'A';
f_Data{1} = 'f';
torquetype_Data{1} = 'Torque Type';
CPUTimeData{1} = 'CPU Time For Simulation (sec)';
MaxErrorData{1} = 'Max Error w';

indexmain=1;
h = waitbar(0, 'Simulating...','Name', 'Simulation Progress');
for index1 = w_0_arr
    w_0 = index1; % Initial Angular Velocity [rad/s]
    for index2 = J_arr
        J = index2; % Rotational Inertia [kg-m^2]
        for index3 = b_arr
            b = index3; % Damping Coefficient [N-m-s/rad]
            for index4 = 1:length(solver_arr)
                solver = solver_arr{index4};
                % check if solver is fixed or variable application
                if strcmp(solver,'ode1') || strcmp(solver,'ode4') % Fixed Application
                    disp(solver);
                    for index5 = dT_arr
                        dT = index5; % Time Step [s]
                        for index6 = A_arr % Constant Applied Torque [N-m]
                            % Tell User Progress on Simulations
                            clc; 
                            disp(['Please wait 10 min while I compute each simulation case']);
                            disp('');
                            waitbar((SimulatedNumber+1)/256,h,fprintf('Computing Simulation %d/256\r', SimulatedNumber+1));

                            A = index6;
                            f = 0;
                            phi = deg2rad(90);
    
                            % Run simulation and track CPU time
                            %for index4 = A_arr
                            %A = index4; % Constant Applied Torque [N-m]
                        
                            
                            % Simulated Rotational Speed Solution
                            tic;
                            simout = sim("P1_demo.slx","Solver",solver,"SolverType", 'Fixed-step');
                            simtime = toc;
                            % Output data
                            W = simout.w.Data;
                            W_DOT = simout.w_dot.Data;
                            T = simout.tout;

                            % Theoretical Rotational Speed Solution Constant Applied Torque (Solution to first order ODE)
                            %w_theor = (A/b)-((A/b)-w_0)*exp(t/J);
                            % Create a vector of every theoretical value
                            t=0:dT:25;
                            W_max_difference = zeros(size(t));
                            for i=1:length(t)
                                curr_t=t(i);
                                w_theor = (A/b)-((A/b)-w_0)*exp(curr_t/J);
                                W_max_difference(i)=abs(W(i)-w_theor);
                            end
                            overallMaxDifference = max(W_max_difference);
                            
    
                            % Increment Indexes
                            SimulatedNumber = SimulatedNumber+1;
                            IndexVal = IndexVal+1;
                            % Update Output With Values
                            SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                            w_0Data{IndexVal} = num2str(w_0);
                            J_Data{IndexVal} = num2str(J);
                            b_Data{IndexVal} = num2str(b);
                            solver_Data{IndexVal} = num2str(solver);
                            dT_Data{IndexVal} = num2str(dT);
                            A_Data{IndexVal} = num2str(A);
                            f_Data{IndexVal} = num2str(f);
                            torquetype_Data{IndexVal} = 'Constant';
                            CPUTimeData{IndexVal} = num2str(simtime);
                            MaxErrorData{IndexVal} = num2str(overallMaxDifference);
                            
    
    
                            %subplot(4,4,indexmain);
                            %plot(W,T);
                            %hold on;
                            %plot(W_DOT,T);
                            %hold on;
                        end
                        for index7 = f_arr  % Sinusoidal Toque - Frequency
                            % Tell User Progress on Simulations
                            clc; 
                            disp(['Please wait 10 min while I compute each simulation case']);
                            disp('');
                            waitbar((SimulatedNumber+1)/256,h,fprintf('Computing Simulation %d/256\r', SimulatedNumber+1));
                            
                            f = index7;
                            phi = 0;
                            A = 1;
    
                            % Run simulation and track CPU time
                            %for index4 = A_arr
                            %A = index4; % Constant Applied Torque [N-m]
                        
                            % Simulated Rotational Speed Solution
                            tic;
                            simout = sim("P1_demo.slx","Solver",solver,"SolverType", 'Fixed-step');
                            simtime = toc;
                            % Output data
                            W = simout.w.Data;
                            W_DOT = simout.w_dot.Data;
                            T = simout.tout;
    
                            % Theoretical Rotational Speed Solution - Sinusoidal Torque (runge kutta 4th order)
                            simout1 = sim("P1_demo.slx","Solver",'ode4',"SolverType", 'Fixed-step');
                            w_theor = simout1.w.Data;
                            % Create a vector of every theoretical value
                            t=0:dT:25;
                            W_max_difference = zeros(size(t));
                            for i=1:length(t)
                                W_max_difference(i)=abs(W(i)-w_theor(i));
                            end
                            overallMaxDifference = max(W_max_difference);
                            
                            
                            % Increment Indexes
                            SimulatedNumber = SimulatedNumber+1;
                            IndexVal = IndexVal+1;
                            % Update Output With Values
                            SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                            w_0Data{IndexVal} = num2str(w_0);
                            J_Data{IndexVal} = num2str(J);
                            b_Data{IndexVal} = num2str(b);
                            solver_Data{IndexVal} = num2str(solver);
                            dT_Data{IndexVal} = num2str(dT);
                            A_Data{IndexVal} = num2str(A);
                            f_Data{IndexVal} = num2str(f);
                            torquetype_Data{IndexVal} = 'Sinusoidal';
                            CPUTimeData{IndexVal} = num2str(simtime);
                            MaxErrorData{IndexVal} = num2str(overallMaxDifference);
                            
    
    
                            %subplot(4,4,indexmain);
                            %plot(W,T);
                            %hold on;
                            %plot(W_DOT,T);
                            %hold on;
                        end
                    end
                else % Variable Application
                    for index6 = A_arr % Constant Applied Torque [N-m]
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 10 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/256,h,fprintf('Computing Simulation %d/256\r', SimulatedNumber+1));

                        A = index6;
                        f = 0;
                        phi = deg2rad(90);

                        % Run simulation and track CPU time
                        %for index4 = A_arr
                        %A = index4; % Constant Applied Torque [N-m]
                    
                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim("P1_demo.slx","Solver",solver,"SolverType", 'Variable-step');
                        simtime = toc;
                        % Output data
                        W = simout.w.Data;
                        W_DOT = simout.w_dot.Data;
                        T = simout.tout;

                        % Theoretical Rotational Speed Solution Constant Applied Torque (Solution to first order ODE)
                        %w_theor = (A/b)-((A/b)-w_0)*exp(t/J);
                        % Create a vector of every theoretical value
                        t=0:dT:25;
                        W_max_difference = zeros(size(t));
                        for i=1:length(t)
                            curr_t=t(i);
                            w_theor = (A/b)-((A/b)-w_0)*exp(curr_t/J);
                            W_max_difference(i)=abs(W(i)-w_theor);
                        end
                        overallMaxDifference = max(W_max_difference);
                        
                        
                        % Increment Indexes
                        SimulatedNumber = SimulatedNumber+1;
                        IndexVal = IndexVal+1;
                        % Update Output With Values
                        SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                        w_0Data{IndexVal} = num2str(w_0);
                        J_Data{IndexVal} = num2str(J);
                        b_Data{IndexVal} = num2str(b);
                        solver_Data{IndexVal} = num2str(solver);
                        dT_Data{IndexVal} = num2str(dT);
                        A_Data{IndexVal} = num2str(A);
                        f_Data{IndexVal} = num2str(f);
                        torquetype_Data{IndexVal} = 'Constant';
                        CPUTimeData{IndexVal} = num2str(simtime);
                        MaxErrorData{IndexVal} = num2str(overallMaxDifference);
                        


                        %subplot(4,4,indexmain);
                        %plot(W,T);
                        %hold on;
                        %plot(W_DOT,T);
                        %hold on;
                    end
                    for index7 = f_arr  % Sinusoidal Toque - Frequency
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 10 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/256,h,fprintf('Computing Simulation %d/256\r', SimulatedNumber+1));
                        
                        f = index7;
                        phi = 0;
                        A = 1;

                        % Run simulation and track CPU time
                        %for index4 = A_arr
                        %A = index4; % Constant Applied Torque [N-m]
                    
                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim("P1_demo.slx","Solver",solver,"SolverType", 'Variable-step');
                        simtime = toc;
                        % Output data
                        W = simout.w.Data;
                        W_DOT = simout.w_dot.Data;
                        T = simout.tout;

                        % Theoretical Rotational Speed Solution - Sinusoidal Torque (runge kutta 4th order)
                        simout1 = sim("P1_demo.slx","Solver",'ode4',"SolverType", 'Fixed-step'); %Note: need to do fixed-step to do runge kutta 4th order
                        w_theor = simout1.w.Data;
                        % Create a vector of every theoretical value
                        t=0:dT:25;
                        W_max_difference = zeros(size(t));
                        for i=1:length(t)
                            W_max_difference(i)=abs(W(i)-w_theor(i));
                        end
                        overallMaxDifference = max(W_max_difference);
                        
                        
                        % Increment Indexes
                        SimulatedNumber = SimulatedNumber+1;
                        IndexVal = IndexVal+1;
                        % Update Output With Values
                        SimulationNumberData{IndexVal} = num2str(SimulatedNumber);
                        w_0Data{IndexVal} = num2str(w_0);
                        J_Data{IndexVal} = num2str(J);
                        b_Data{IndexVal} = num2str(b);
                        solver_Data{IndexVal} = num2str(solver);
                        dT_Data{IndexVal} = num2str(dT);
                        A_Data{IndexVal} = num2str(A);
                        f_Data{IndexVal} = num2str(f);
                        torquetype_Data{IndexVal} = 'Sinusoidal';
                        CPUTimeData{IndexVal} = num2str(simtime);
                        MaxErrorData{IndexVal} = num2str(overallMaxDifference);

                        %subplot(4,4,indexmain);
                        %plot(W,T);
                        %hold on;
                        %plot(W_DOT,T);
                        %hold on;
                    end
                end
            end
            %title(['w_0=',num2str(w_0),', J=',num2str(J),', b=',num2str(b),', A=',num2str(A)]);
            %legend('dT=0.001','dT=0.1','dT=1');
            %xlabel('Angular Velocity [rad/s]');
            %ylabel('Time [s]');
            %indexmain=indexmain+1;
        end
    end
end
%hold off;
%hFig = gcf;  % Get the current figure handle
%set(hFig, 'WindowState', 'maximized');

close(h);
disp('');
disp('Simulation complete. Outputting data now.');


filename = 'Project1_OutputData.csv';
dataTable = table(SimulationNumberData',w_0Data',J_Data',b_Data',solver_Data',dT_Data',A_Data',f_Data',torquetype_Data',CPUTimeData',MaxErrorData');
writetable(dataTable, filename, 'WriteVariableNames', false);
winopen(filename);


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

