%Plotting urban drive cycle
simout = sim("p3_car_urban.slx");
sim_vel = simout.vel.Data;
sim_time = simout.tout;

figure;
plot(sim_time, sim_vel*(1/mph2mps), 'b') % Remember, drive cycles are mph
hold on
plot(TimeUrban, DriveDataUrban, '--r') 
plot(TimeUrban, (DriveDataUrban)+3, '--k') 
plot(TimeUrban, (DriveDataUrban)-3, '--k') 
xlabel("Time (s)")
ylabel("Velocity (mph)") 
legend("Sim Velocity", "Drive Cycle Velocity") % , "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time for Urban Drive Cycle")


% For seeing how large the errors are
error = zeros(length(TimeUrban),1);
for j = 1:length(TimeUrban)
    time_dc = TimeUrban(j);
    vel_dc = DriveDataUrban(j);
    for i = 1:length(sim_time)
        time_s = sim_time(i);
        vel_s = sim_vel(i);

        if time_s == time_dc
            err = vel_dc - vel_s;
            error(j) = err;
        else
        end
    end
end


sim_torque = simout.torque.Data; %units: Nm
sim_angularvelocity = simout.angularvelocity.Data; %units: rpm

power = zeros(size(sim_time));
for i = 1:length(sim_time)
    angularvelocity_rads = (sim_angularvelocity(i) * ((2*pi) / 60)); % angular velocity in rad/s
    power(i) = sim_torque(i)*angularvelocity_rads; % power in Watts
end

%Energy consumed in Joules
energy = sum(power)*0.001; %Interval between each power measurment is 0.001 second
disp(['Total energy consumed for Urban: ' num2str(energy) ' Joules']);


%Plotting highway drive cycle
simout_highway = sim("p3_car_highway.slx");
sim_vel_highway = simout_highway.vel.Data;
sim_time_highway = simout_highway.tout;

figure;
plot(sim_time_highway, sim_vel_highway*(1/mph2mps), 'b') % Remember, drive cycles are mph
hold on
plot(TimeHighway, DriveDataHighway, '--r') 
plot(TimeHighway, (DriveDataHighway)+3, '--k') 
plot(TimeHighway, (DriveDataHighway)-3, '--k') 
xlabel("Time (s)")
ylabel("Velocity (mph)") 
legend("Sim Velocity", "Drive Cycle Velocity") % , "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time for Highway Drive Cycle")


% For seeing how large the errors are
error = zeros(length(TimeHighway),1);
for j = 1:length(TimeHighway)
    time_dc = TimeHighway(j);
    vel_dc = DriveDataHighway(j);
    for i = 1:length(sim_time_highway)
        time_s = sim_time_highway(i);
        vel_s = sim_vel_highway(i);

        if time_s == time_dc
            err = vel_dc - vel_s;
            error(j) = err;
        else
        end
    end
end



sim_torque_highway = simout_highway.torque.Data; %units: Nm
sim_angularvelocity_highway = simout_highway.angularvelocity.Data; %units: rpm

power_highway = zeros(size(sim_time_highway));
for j = 1:length(sim_time_highway)
    angularvelocity_rads_highway = (sim_angularvelocity_highway(j) * ((2*pi) / 60)); %angular velocity in rad/s
    power_highway(j) = sim_torque_highway(j)*angularvelocity_rads_highway; %power in Watts
end

%Energy consumed in Joules
energy_highway = sum(power_highway)*0.001; %Interval between each power measurment is 0.001 second
disp(['Total energy consumed for Highway: ' num2str(energy_highway) ' Joules']);
