% Define parameters
radius = 200; % radius of curved sections
straight_length = 900; % length of straight sections
width = 15;
carspeed = 500; % a lower number is faster


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
idx1 = 0;
degreelhs = linspace(3*pi/2, pi/2, carspeed);
idx2 = 0;

% Animate the car driving around the track
for i = 1:length(xpath)
    x = xpath(i);
    y = ypath(i);
    % Update car image
    clearpoints(h);
    addpoints(h, x, y);

    %generate un-rotated car at new location x,y
    car = [x - w, y - w/2; x + w, y-w/2; x+w, y+w/2; x-w, y+w/2];

    % rotate car according to slope of curve 
    if y == 0 || y == radius*2
        car2 = rotate(car'-[x;y], 0)' + [x;y]';
    else
        if idx1<(carspeed-2)
            idx1=idx1+1;
            car2 = rotate(car'-[x;y],degreerhs(idx1)+pi/2)' + [x;y]';
        else
            idx2=idx2+1;
            car2 = rotate(car'-[x;y],degreelhs(idx2)+pi/2)' + [x;y]';
        end
    end

    %update car image
    a.XData = car2(:,1);
    a.YData = car2(:,2);

    % Update vehicle's path
    addpoints(path_line, x, y);

    drawnow limitrate nocallbacks;
end

hold off;
