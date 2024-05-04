%importdata('raceStatOld.m');
%importdata('raceStat.m');
importdata('p4init.m');
importdata('rotate.m');


% **Add in a way to specify number of laps desired
% Request input from the user
%clc;
%disp("I will simulate a racecar driving as fast as possible around a track using real physics contraints.");
%disp("");
%response = input('Please enter the # of laps for the race: ');
%disp("Thank you. Please wait a few seconds while I generate the simulation.");

pure_pursuit_lookaheaddist = 5;


% Define parameters
radius = 200; % radius of curved sections
straight_length = 900; % length of straight sections
width = 15;
carspeed = 100; % a lower number is faster
n = 1; % Number of laps

% Define angles for curved sections
theta1 = linspace(pi/2, 3*pi/2, carspeed); % 180 degrees
theta2 = linspace(-pi/2, pi/2, carspeed); % 180 degrees

% Define the x and y coordinates for each section of the track
x_straight1 = linspace(0, straight_length, carspeed);
y_straight1 = zeros(1, carspeed);
x_curve1 = straight_length + radius * cos(theta2);
y_curve1 = radius * sin(theta2) + radius;
x_straight2 = linspace(straight_length, 0, carspeed);
y_straight2 = ones(1, carspeed) * (2 * radius);
x_curve2 = radius * cos(theta1);
y_curve2 = radius * sin(theta1) + radius;

% Combine parametric sections
xpath = [x_straight1, x_curve1, x_straight2, x_curve2];
ypath = [y_straight1, y_curve1, y_straight2, y_curve2];

% Calculate the total time for one lap
total_time_one_lap = length(xpath);
% Calculate the total time for 'n' laps
total_time_n_laps = total_time_one_lap * n;
% Generate the time vector for 'n' laps
t = linspace(0, total_time_n_laps, total_time_n_laps);

% Plot the race track
plot(xpath, ypath, 'LineWidth', width, 'Color', [0.8, 0.8, 0.8]);
axis equal;
title('Race Track');
xlabel('X (m)');
ylabel('Y (m)');
hold on;

% Initialize animated line for the track
h = animatedline;
axis([-radius*2 straight_length+radius*2 -radius radius*3])

% Initialize animated line for car path
path_line = animatedline('Color', 'black', 'LineWidth', 1);

% Define the car parameters
w = 30; % Width of the car
car = [-w/2, -w; w/2, -w; w/2, w; -w/2, w]'; % Car vertices
a = patch('XData',car(:,1),'YData',car(:,2));
a.EdgeColor = [1 0 0];
a.FaceColor = [1 0 0];

degreerhs = linspace(pi/2, -pi/2, carspeed);
degreelhs = linspace(3*pi/2, pi/2, carspeed);

%{
Preallocate arrays that will store the x-values, y-values, and t-values
for all laps around the track. The arrays are instantiated with the
value 9999999999 (a random large number so we can delete untouched
values later) for carspeed*500 elements (the amount of times a value is
plotted for a lap*the max amount of laps that we are saying is 1000).
%}
iteratorindex = 0;
xall = 9999999999*ones(1,carspeed*1000);
yall = 9999999999*ones(1,carspeed*1000);
%tall = 9999999999*ones(1,carspeed*1000);


% Animate the car driving around the track


sim("P4_v4.slx");
tic;
simout = sim("P4_v4.slx");
for lap = 1:n

    idx1 = 0;
    idx2 = 0;
    x = simout.X.Data;
    y = simout.Y.Data;
    vx = simout.vx.Data;
    vy = simout.vy.Data;
    tout = simout.tout;

    % Iterate through each element in the list
    for i = 1:numel(x)
        % Access the i-th element of the list
        current_x = x(i);
        current_y = y(i);
        current_vx = vx(i);
        current_vy = vy(i);
        
        totalv=sqrt(current_vx^2+current_vy^2);
        % Add in a way to display car speed:
        %disp(['Current Speed = ',num2str(totalv),' m/s']);

        % Update car image
        clearpoints(h);
        addpoints(h, current_x, current_y);

        %generate un-rotated car at new location x,y
        car = [current_x - w, current_y - w/2; current_x + w, current_y-w/2; current_x+w, current_y+w/2; current_x-w, current_y+w/2];

        % rotate car according to slope of curve 
        if (current_x <= 900) && (current_x >= 0) %straightaway
            car2 = rotate(car'-[current_x;current_y], 0)' + [current_x;current_y]';
        else
            if (current_x > 900) % rhs
                idx1=idx1+1;
                degreerhsnew = -atan((-(200-current_y))/(current_x-900));
                car2 = rotate(car'-[current_x;current_y],degreerhsnew+pi/2)' + [current_x;current_y]';
            else %lhs
                idx2=idx2+1;
                degreelhsnew = 2*pi/2+atan((-(200-current_y))/(-current_x));
                car2 = rotate(car'-[current_x;current_y],degreelhsnew+pi/2)' + [current_x;current_y]';
            end
        end

        %update car image
        a.XData = car2(:,1);
        a.YData = car2(:,2);

        % Update vehicle's path
        addpoints(path_line, current_x, current_y);

        %append output value array for racestat file
        iteratorindex=iteratorindex+1;
        xall(iteratorindex) = current_x;
        yall(iteratorindex) = current_y;
        %tall(iteratorindex) = (tout(end))*lap; % here we do (last tout value)*(number of laps) because each sim is one lap

        drawnow limitrate nocallbacks;
    end
end


% Usage: rs = raceStat(X,Y,t,path)
%
% Inputs:   X, Y are coordinates from your vehicle simulations. 
%           t is the set times corresponding to X and Y
%           path is a structure of with fields "width" (width of the 
%                  track), "l_st" (length of the straight away), and 
%                  "radius" (radius of the curved section)

% remove values that were not touched
xall = xall(xall ~= 9999999999);
yall = yall(yall ~= 9999999999);
%tall = tall(tall ~= 9999999999);

% return values
X = xall;
Y = yall;
%t = tall;
path = struct('width', width, 'l_st', straight_length, 'radius', radius);


hold off;