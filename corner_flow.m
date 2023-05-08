% Parameters
U = 1; % Freestream velocity
k = 1; % Vortex strength
A = 1; % Constant related to wedge angle and vortex core size
R = 0.4; % Vortex core radius
n = 100; % Number of grid points in each direction
xmin = -1; xmax = 1; % x-limits of plot
ymin = -1/0.6; ymax = 1/0.6; % y-limits of plot

% Create grid
x = linspace(xmin, xmax, n);
y = linspace(ymin, ymax, n);
[X,Y] = meshgrid(x,y);

% Compute stream function
r = sqrt(X.^2 + Y.^2);
theta = atan2(Y,X);


%Psi = U*r.*sin(theta) - k*A*cos(theta).*sin(2*theta)/(2*pi) .* log(r/R);
Psi = 3*X.^2.*Y - Y.^3;

U = 3*X.^2-3*Y.^2;
V = -(-6*X*Y);

uu = @(x,y) 3*x^2-3*y^2;
vv = @(x,y) (-6*x*y);



% Initialize particle position and velocity
r = 0;
v = 0;

% Integrate equations of motion
t = 0;
rx = 0.99; %% horizontal placement of particle
ry = 1.65; %%vertical placement of particle 

xs = [rx];
ys = [ry];
ts = [t];

tmax = 2.5;
dt = 0.05;
while t < tmax
    vx = uu(rx,ry);
    vy = vv(rx,ry);
  
    rx = rx + vx*dt;
    ry = ry + vy*dt; 

    % Increment time
    t = t + dt;
    xs = [xs rx];
    ys = [ys ry];
end

[px,py] = gradient(Psi);

figure(1)
% Plot streamlines
hold on
contour(X,Y,Psi,100);
yline(0)
xl1 = 0:0.1:1;
plot(xl1,xl1./0.58,Color='k')

xlim([0 xmax]);
ylim([0 ymax]);
xlabel('x');
ylabel('y');
title('Streamlines of Moffatt vortices in a wedge');
hold off
[startx,starty] = meshgrid(0.05:0.1:0.85,0.01:0.08:1);

verts = stream2(X,Y,py,-px,startx,starty);

figure(2)
quiver(X,Y,py,-px);
lineobj = streamline(verts);


figure(3)
% Plot streamlines
hold on
contour(X,Y,Psi,100);
yline(0)
xl1 = 0:0.1:1;
plot(xl1,xl1./0.58,Color='k')
plot(xs,ys,'r.',MarkerSize=10)
xlim([0 xmax]);
ylim([0 ymax]);
xlabel('x');
ylabel('y');
title('Streamlines of Moffatt vortices in a wedge');
hold off
