% Project 1 - Part 1


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
IndexValMaxSimOde1 = 0;
IndexValMaxSimOde4 = 0;
IndexValMaxSimOde45 = 0;
IndexValMaxSimOde23tb = 0;
IndexCompile = 0;

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
MaxSimMaxSimErrorOde1 = ones(1, 300)*99999;
CPUTimeTakenOde1 = ones(1, 300)*99999;
MaxSimTimeStepOde1 = ones(1, 300)*99999;
MaxSimMaxSimErrorOde4 = ones(1, 300)*99999;
CPUTimeTakenOde4 = ones(1, 300)*99999;
MaxSimTimeStepOde4 = ones(1, 300)*99999;
MaxSimMaxSimErrorOde45 = ones(1, 300)*99999;
CPUTimeTakenOde45 = ones(1, 300)*99999;
MaxSimTimeStepOde45 = ones(1, 300)*99999;
MaxSimMaxSimErrorOde23tb = ones(1, 300)*99999;
CPUTimeTakenOde23tb = ones(1, 300)*99999;
MaxSimTimeStepOde23tb = ones(1, 300)*99999;
EigenValuesOde1 = ones(1, 300)*99999;
EigenValuesOde4 = ones(1, 300)*99999;
EigenValuesOde45= ones(1, 300)*99999;
EigenValuesOde23tb = ones(1, 300)*99999;
FreqValuesOde1 = ones(1, 300)*99999;
FreqValuesOde4 = ones(1, 300)*99999;
FreqValuesOde45 = ones(1, 300)*99999;
FreqValuesOde23tb = ones(1, 300)*99999;

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
h = waitbar(0, 'Simulation Progress');
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
                            % ---------------------------------------------------------------------------------------- (Fixed Time Step, Constant Torque)
                            % Tell User Progress on Simulations
                            if IndexCompile == 0
                                clc;
                                disp('Please wait 15 sec while I compile the Simulink model');
                                disp('');
                                disp('Compiling...');
                                simout = sim("P1_Part1.slx","Solver",solver,"SolverType", 'Fixed-step');
                                IndexCompile= IndexCompile+1;
                            end
                            
                            clc; 
                            disp(['Please wait 10 min while I compute each simulation case']);
                            disp('');
                            waitbar((SimulatedNumber+1)/256,h,'Simulation Progress');
                            fprintf('Computing Simulation %d/256\r', SimulatedNumber+1);

                            A = index6;
                            f = 0;
                            phi = deg2rad(90);
    
                            % Run simulation and track CPU time
                            %for index4 = A_arr
                            %A = index4; % Constant Applied Torque [N-m]
                        
                            
                            % Simulated Rotational Speed Solution
                            tic;
                            simout = sim("P1_Part1.slx","Solver",solver,"SolverType", 'Fixed-step');
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

                            % Plot Data
                            if strcmp(solver,'ode1')
                                IndexValMaxSimOde1 = IndexValMaxSimOde1+1;
                                MaxSimMaxSimErrorOde1(IndexValMaxSimOde1) = overallMaxDifference;
                                CPUTimeTakenOde1(IndexValMaxSimOde1) = simtime;
                                MaxSimTimeStepOde1(IndexValMaxSimOde1) = dT;
                                EigenValuesOde1(IndexValMaxSimOde1)=A;
                            else %ode4
                                IndexValMaxSimOde4 = IndexValMaxSimOde4+1;
                                MaxSimMaxSimErrorOde4(IndexValMaxSimOde4) = overallMaxDifference;
                                CPUTimeTakenOde4(IndexValMaxSimOde4) = simtime;
                                MaxSimTimeStepOde4(IndexValMaxSimOde4) = dT;
                                EigenValuesOde4(IndexValMaxSimOde4)=A;
                            end
    
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
                            % ---------------------------------------------------------------------------------------- (Fixed Time Step, Sinusoidal Torque)
                            % Tell User Progress on Simulations
                            clc; 
                            disp(['Please wait 10 min while I compute each simulation case']);
                            disp('');
                            waitbar((SimulatedNumber+1)/256,h,'Simulation Progress');
                            fprintf('Computing Simulation %d/256\r', SimulatedNumber+1);
                            
                            f = index7;
                            phi = 0;
                            A = 1;
    
                            % Run simulation and track CPU time
                            %for index4 = A_arr
                            %A = index4; % Constant Applied Torque [N-m]
                        
                            % Simulated Rotational Speed Solution
                            tic;
                            simout = sim("P1_Part1.slx","Solver",solver,"SolverType", 'Fixed-step');
                            simtime = toc;
                            % Output data
                            W = simout.w.Data;
                            W_DOT = simout.w_dot.Data;
                            T = simout.tout;
                            
    
                            % Theoretical Rotational Speed Solution - Sinusoidal Torque (runge kutta 4th order)
                            simout1 = sim("P1_Part1.slx","Solver",'ode4',"SolverType", 'Fixed-step');
                            w_theor = simout1.w.Data;
                            % Create a vector of every theoretical value
                            t=0:dT:25;
                            W_max_difference = zeros(size(t));
                            for i=1:length(t)
                                W_max_difference(i)=abs(W(i)-w_theor(i));
                            end
                            overallMaxDifference = max(W_max_difference);

                            % Plot Data
                            if strcmp(solver,'ode1')
                                IndexValMaxSimOde1 = IndexValMaxSimOde1+1;
                                MaxSimMaxSimErrorOde1(IndexValMaxSimOde1) = overallMaxDifference;
                                CPUTimeTakenOde1(IndexValMaxSimOde1) = simtime;
                                MaxSimTimeStepOde1(IndexValMaxSimOde1) = dT;
                                FreqValuesOde1(IndexValMaxSimOde1)=f;
                            else %ode4
                                IndexValMaxSimOde4 = IndexValMaxSimOde4+1;
                                MaxSimMaxSimErrorOde4(IndexValMaxSimOde4) = overallMaxDifference;
                                CPUTimeTakenOde4(IndexValMaxSimOde4) = simtime;
                                MaxSimTimeStepOde4(IndexValMaxSimOde4) = dT;
                                FreqValuesOde4(IndexValMaxSimOde4)=f;
                            end
                            
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
                        % ---------------------------------------------------------------------------------------- (Variable Time Step, Constant Torque)
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 10 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/256,h,'Simulation Progress');
                        fprintf('Computing Simulation %d/256\r', SimulatedNumber+1);

                        A = index6;
                        f = 0;
                        phi = deg2rad(90);

                        % Run simulation and track CPU time
                        %for index4 = A_arr
                        %A = index4; % Constant Applied Torque [N-m]
                    
                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim("P1_Part1.slx","Solver",solver,"SolverType", 'Variable-step');
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

                        % Plot Data
                        if strcmp(solver,'ode45')
                            IndexValMaxSimOde45 = IndexValMaxSimOde45+1;
                            MaxSimMaxSimErrorOde45(IndexValMaxSimOde45) = overallMaxDifference;
                            CPUTimeTakenOde45(IndexValMaxSimOde45) = simtime;
                            MaxSimTimeStepOde45(IndexValMaxSimOde45) = dT;
                            EigenValuesOde45(IndexValMaxSimOde45)=A;
                        else %ode23tb
                            IndexValMaxSimOde23tb = IndexValMaxSimOde23tb+1;
                            MaxSimMaxSimErrorOde23tb(IndexValMaxSimOde23tb) = overallMaxDifference;
                            CPUTimeTakenOde23tb(IndexValMaxSimOde23tb) = simtime;
                            MaxSimTimeStepOde23tb(IndexValMaxSimOde23tb) = dT;
                            EigenValuesOde23tb(IndexValMaxSimOde23tb)=A;
                        end
                        
                        
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
                        % ---------------------------------------------------------------------------------------- (Variable Time Step, Sinusoidal Torque)
                        % Tell User Progress on Simulations
                        clc; 
                        disp(['Please wait 10 min while I compute each simulation case']);
                        disp('');
                        waitbar((SimulatedNumber+1)/256,h,'Simulation Progress');
                        fprintf('Computing Simulation %d/256\r', SimulatedNumber+1);
                        
                        f = index7;
                        phi = 0;
                        A = 1;

                        % Run simulation and track CPU time
                        %for index4 = A_arr
                        %A = index4; % Constant Applied Torque [N-m]
                    
                        % Simulated Rotational Speed Solution
                        tic;
                        simout = sim("P1_Part1.slx","Solver",solver,"SolverType", 'Variable-step');
                        simtime = toc;
                        % Output data
                        W = simout.w.Data;
                        W_DOT = simout.w_dot.Data;
                        T = simout.tout;

                        % Theoretical Rotational Speed Solution - Sinusoidal Torque (runge kutta 4th order)
                        simout1 = sim("P1_Part1.slx","Solver",'ode4',"SolverType", 'Fixed-step'); %Note: need to do fixed-step to do runge kutta 4th order
                        w_theor = simout1.w.Data;
                        % Create a vector of every theoretical value
                        t=0:dT:25;
                        W_max_difference = zeros(size(t));
                        for i=1:length(t)
                            W_max_difference(i)=abs(W(i)-w_theor(i));
                        end
                        overallMaxDifference = max(W_max_difference);

                        % Plot Data
                        if strcmp(solver,'ode45')
                            IndexValMaxSimOde45 = IndexValMaxSimOde45+1;
                            MaxSimMaxSimErrorOde45(IndexValMaxSimOde45) = overallMaxDifference;
                            CPUTimeTakenOde45(IndexValMaxSimOde45) = simtime;
                            MaxSimTimeStepOde45(IndexValMaxSimOde45) = dT;
                            FreqValuesOde45(IndexValMaxSimOde45)=f;
                        else %ode23tb
                            IndexValMaxSimOde23tb = IndexValMaxSimOde23tb+1;
                            MaxSimMaxSimErrorOde23tb(IndexValMaxSimOde23tb) = overallMaxDifference;
                            CPUTimeTakenOde23tb(IndexValMaxSimOde23tb) = simtime;
                            MaxSimTimeStepOde23tb(IndexValMaxSimOde23tb) = dT;
                            FreqValuesOde23tb(IndexValMaxSimOde23tb)=f;
                        end
                        
                        
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

% Clean Up Plot Data to Delete Blank Elements
MaxSimTimeStepOde1 = MaxSimTimeStepOde1(MaxSimTimeStepOde1 ~= 99999);
MaxSimTimeStepOde4 = MaxSimTimeStepOde4(MaxSimTimeStepOde4 ~= 99999);
MaxSimMaxSimErrorOde1 = MaxSimMaxSimErrorOde1(MaxSimMaxSimErrorOde1 ~= 99999);
MaxSimMaxSimErrorOde4 = MaxSimMaxSimErrorOde4(MaxSimMaxSimErrorOde4 ~= 99999);
CPUTimeTakenOde1 = CPUTimeTakenOde1(CPUTimeTakenOde1 ~= 99999);
CPUTimeTakenOde4 = CPUTimeTakenOde4(CPUTimeTakenOde4 ~= 99999);
MaxSimTimeStepOde23tb = MaxSimTimeStepOde23tb(MaxSimTimeStepOde23tb ~= 99999);
MaxSimTimeStepOde45 = MaxSimTimeStepOde45(MaxSimTimeStepOde45 ~= 99999);
MaxSimMaxSimErrorOde23tb = MaxSimMaxSimErrorOde23tb(MaxSimMaxSimErrorOde23tb ~= 99999);
MaxSimMaxSimErrorOde45 = MaxSimMaxSimErrorOde45(MaxSimMaxSimErrorOde45 ~= 99999);
CPUTimeTakenOde23tb = CPUTimeTakenOde23tb(CPUTimeTakenOde23tb ~= 99999);
CPUTimeTakenOde45 = CPUTimeTakenOde45(CPUTimeTakenOde45 ~= 99999);
EigenValuesOde1 = EigenValuesOde1(EigenValuesOde1 ~= 99999);
EigenValuesOde4 = EigenValuesOde4(EigenValuesOde4 ~= 99999);
EigenValuesOde45 = EigenValuesOde45(EigenValuesOde45 ~= 99999);
EigenValuesOde23tb = EigenValuesOde23tb(EigenValuesOde23tb ~= 99999);
FreqValuesOde1 = FreqValuesOde1(FreqValuesOde1 ~= 99999);
FreqValuesOde4 = FreqValuesOde4(FreqValuesOde4 ~= 99999);
FreqValuesOde45 = FreqValuesOde45(FreqValuesOde45 ~= 99999);
FreqValuesOde23tb = FreqValuesOde23tb(FreqValuesOde23tb ~= 99999);

% Plot Max Simulation Error vs Time Step
figure;  % Create a new figure
scatter(MaxSimTimeStepOde1, MaxSimMaxSimErrorOde1, 'blue', 'o', 'filled', 'SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(MaxSimTimeStepOde4, MaxSimMaxSimErrorOde4, 'red', 'o', 'filled', 'SizeData', 10);
hold off;
% Add labels and title
xlabel('Time Step (s)');
ylabel('Max Simulation Error (w)');
title('Max Simulation Error vs Time Step');
% Add a legend
legend('ode1', 'ode4');

% CPU Time Taken vs Time Step
figure;  % Create a new figure
scatter(MaxSimTimeStepOde1, CPUTimeTakenOde1, 'blue', 'o', 'filled', 'SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(MaxSimTimeStepOde4, CPUTimeTakenOde4, 'red', 'o', 'filled', 'SizeData', 10);
hold off;
% Add labels and title
xlabel('Time Step (s)');
ylabel('CPU Time Taken (s)');
title('CPU Time Taken vs Time Step');
% Add a legend
legend('ode1', 'ode4');

% Max Simulation Error vs CPU Time Taken
figure;  % Create a new figure
scatter(CPUTimeTakenOde1, MaxSimMaxSimErrorOde1, 'blue', 'o', 'filled', 'SizeData', 10);
hold on;  % Hold the current plot so that the next plots are added to it
scatter(CPUTimeTakenOde4, MaxSimMaxSimErrorOde4, 'red', 'o', 'filled', 'SizeData', 10);
scatter(CPUTimeTakenOde45, MaxSimMaxSimErrorOde45, 'yellow', 'o', 'filled', 'SizeData', 10);
scatter(CPUTimeTakenOde23tb, MaxSimMaxSimErrorOde23tb, 'green', 'o', 'filled', 'SizeData', 10);
hold off;
% Add labels and title
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
title('Max Simulation Error vs CPU Time Taken');
% Add a legend
legend('ode1', 'ode4','ode45','ode23tb');

% Eigen Values Heat
figure;  % Create a new figure
subplot(2, 2, 1);
heatmap(CPUTimeTakenOde1, MaxSimMaxSimErrorOde1, EigenValuesOde1);
title('Ode1');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 2);
heatmap(CPUTimeTakenOde4, MaxSimMaxSimErrorOde4, EigenValuesOde4);
title('Ode4');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 3);
heatmap(CPUTimeTakenOde45, MaxSimMaxSimErrorOde45, EigenValuesOde45);
title('Ode4');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 4);
heatmap(CPUTimeTakenOde23tb, MaxSimMaxSimErrorOde23tb, EigenValuesOde23tb);
title('Ode23tb');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
sgtitle('Contour plot of constant system eigen values with Max Simulation Error vs CPU Time Taken');

% Freq Values Heat
figure;  % Create a new figure
subplot(2, 2, 1);
heatmap(CPUTimeTakenOde1, MaxSimMaxSimErrorOde1, FreqValuesOde1);
title('Ode1');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 2);
heatmap(CPUTimeTakenOde4, MaxSimMaxSimErrorOde4, FreqValuesOde4);
title('Ode4');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 3);
heatmap(CPUTimeTakenOde45, MaxSimMaxSimErrorOde45, FreqValuesOde45);
title('Ode4');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
subplot(2, 2, 4);
heatmap(CPUTimeTakenOde23tb, MaxSimMaxSimErrorOde23tb, FreqValuesOde23tb);
title('Ode23tb');
xlabel('CPU Time Taken (s)');
ylabel('Max Simulation Error (w)');
sgtitle('Contour plot of constant input frequencies with Max Simulation Error vs CPU Time Taken');

