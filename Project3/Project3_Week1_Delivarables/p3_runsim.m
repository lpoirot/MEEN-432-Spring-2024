%Plotting urban drive cycle
simout = sim("p3_w1_car.slx");
sim_vel = simout.vel.Data;
sim_time = simout.tout;

figure;
plot(sim_time, sim_vel*(1/mph2mps), 'b') % Remember, drive cycles are mph
hold on
plot(Time, DriveData, '--r') 
plot(Time, (DriveData)+3, '--k') 
plot(Time, (DriveData)-3, '--k') 
xlabel("Time (s)")
ylabel("Velocity (mph)") 
legend("Sim Velocity", "Drive Cycle Velocity") % , "3 mph Error Band")
title("Simulated Vehicle Velocity vs Time for Urban Drive Cycle")


% For seeing how large the errors are
error = zeros(length(Time),1);
for j = 1:length(Time)
    time_dc = Time(j);
    vel_dc = DriveData(j);
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


%Plotting highway drive cycle
simout_highway = sim("p3_w1_car_highway.slx");
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
