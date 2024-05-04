P4init;
p4sim;
simout = sim("P4_v4.slx","StopTime","3600")
carX = simout.X.Data;
carY = simout.Y.Data;
tout = simout.tout;

% Race Statistics
race = raceStat(carX, carY, tout, path, simout);
disp(race)