%% Steady State Solution With Point Particle Swimmer
length_wall = 1; %in meters
n = 100;
alpha = pi/12; %half the angle of the wedge
lambda = [5]; % different eddies
theta = linspace(-alpha,alpha,n);
dtheta = theta(3)-theta(2);
r = linspace(0,1,n); r(1) = 1E-3; %buffer to avoid singularity
dr = r(3)-r(2); dt = linspace(0,1,n); %time steps
A = 1; B = 1; C = 1; D = 1;
a = 2E-6; %radius of microswimmer in meters

%streamline function
psi = ones(length(r),length(theta));

xx = ones(length(r),length(theta)); 
yy = ones(length(r),length(theta)); 

U = cell(n,n); Ur = cell(n,n); Urr = cell(n,n); Uthetatheta = cell(n,n);

for i3 = 1:n
    for j3 = 1:n
        U{i3,j3} = [0,0];Ur{i3,j3} = [0,0];Urr{i3,j3} = [0,0];Uthetatheta{i3,j3} = [0,0];
    end
end

for d = 1:length(lambda)
    for i2 = 1:length(r)
        for j = 1:length(theta)
            psi(i2,j) = (r(i2)^lambda(d))*(A*cos(lambda(d)*theta(j))+B*sin(lambda(d)*theta(j)) + ...
                C*cos((lambda(d)-2)*theta(j))+D*sin((lambda(d)-2)*theta(j)));
            xx(i2,j) = r(i2)*cos(theta(j));
            yy(i2,j) = r(i2)*sin(theta(j));
            if i2 < n && i2 > 1
                U{i2+1,j} = [u_rf(r(i2),lambda(d),theta(j),A,B,C,D),...
                    u_thetaf(r(i2),lambda(d),theta(j),A,B,C,D)]; %At (r,theta), returns <U(r),U(theta)>
                [drUr,drUt] = dru(r(i2),lambda(d),theta(j),A,B,C,D);
                Ur{i2+1,j} = [drUr,drUt]; %At (r,theta) returns <dU(r)/dr,dU(theta)/dr>
                [d2ru_r,d2ru_theta] = d2ru(r(i2),lambda(d),theta(j),A,B,C,D);
                Urr{i2+1,j} = [d2ru_r,d2ru_theta]; %At (r,theta) returns <d^2U(r)/dr^2,d^2U(theta)/dr^2>
                [d2tu_r,d2tu_theta] = d2tu(r(i2),lambda(d),theta(j),A,B,C,D);
                Uthetatheta{i2+1,j} = [d2tu_r,d2tu_theta]; %At (r,theta) returns <d^2U(r)/dtheta^2,d^2U(theta)/dtheta^2>
            else
                U{i2+1,j} = [0,0]; %BC
                %% ASSUMING Neumann conditions of 1st and 2nd order are 0 - NEED TO CHECK
                Ur{i2+1,j} = [0,0];
                Urr{i2+1,j} = [0,0];
                Uthetatheta{i2+1,j} = [0,0];
            end
        end
    end
    U(1,:) = []; Ur(1,:) = []; Urr(1,:) = []; Uthetatheta(1,:) = [];
    %% Influence of Faxen's Laws
    if i2 < length(r)
        U{i2+1,j} = U{i2+1,j}+(a^2/6)*Urr{i2+1,j}+(a^2/(6*r(i2)))*Ur{i2+1,j}+(a^2/(6*r(i2)^2))*Uthetatheta{i2+1,j};
    else
        U{n,j} = U{n,j}+(a^2/6)*Urr{n,j}+(a^2/(6*r(i2)))*Ur{n,j}+(a^2/(6*r(i2)^2))*Uthetatheta{n,j};
    end
    %% RK-4 Implementation 
    %initialize in middle of the array
    cr(1) = r(n/2); ctheta(1) = theta(n/2); 
    k1r = dt(1)*U{n/2,n/2}(1); newr = cr(1)+k1r; k1theta = dt(1)*U{n/2,n/2}(2); newtheta = ctheta(1)+k1theta;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k2r = dt(1)*U{indr,indtheta}(1); newr = cr(1)+k1r/2; k2theta = dt(1)*U{indr,indtheta}(2); newtheta = ctheta(1)+k1theta;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k3r = dt(1)*U{indr,indtheta}(1); newr = cr(1)+k2r/2; k3theta = dt(1)*U{indr,indtheta}(2); newtheta = ctheta(1)+k2theta/2;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k4r = dt(1)*U{indr,indtheta}(1); newr = cr(1)+k3r/2; k4theta = dt(1)*U{indr,indtheta}(2); newtheta = ctheta(1)+k3theta/2;
    for t = 2:length(dt)
        cr(t) = cr(t-1)+(1/6)*(k1r+2*k2r+2*k3r+2*k4r);
        ctheta(t) = ctheta(t-1)+(1/6)*(k1theta+2*k2theta+2*k3theta+2*k4theta);
        % use function to find the closest r index and theta index
        [indr,indtheta] = closestField(r,theta,cr(t),ctheta(t)); 
        % find the velocity field at that r and theta
        k1r = dt(t)*U{indr,indtheta}(1); newr = cr(t)+k1r; k1theta = dt(t)*U{indr,indtheta}(2); newtheta = ctheta(t)+k1theta;
        [indr,indtheta] = closestField(r,theta,newr,newtheta); 
        k2r = dt(t)*U{indr,indtheta}(1); newr = cr(t)+k1r/2; k2theta = dt(t)*U{indr,indtheta}(2); newtheta = ctheta(t)+k1theta;
        [indr,indtheta] = closestField(r,theta,newr,newtheta); 
        k3r = dt(t)*U{indr,indtheta}(1); newr = cr(t)+k2r/2; k3theta = dt(t)*U{indr,indtheta}(2); newtheta = ctheta(t)+k2theta/2;
        [indr,indtheta] = closestField(r,theta,newr,newtheta); 
        k4r = dt(t)*U{indr,indtheta}(1); newr = cr(t)+k3r/2; k4theta = dt(t)*U{indr,indtheta}(2); newtheta = ctheta(t)+k3theta/2;
    end
    for t = 2:length(dt)
    end
    psi = psi/norm(psi,2);
    figure(d)
    hold on
    contour(xx,yy,psi)
    scatter(cr,ctheta,'filled')
    f = gcf;
    title_temp = strcat('eigenvalue ',int2str(d), 'stream function','.jpg');
    saveas(f,title_temp)
    title{j} = title_temp;
    hold off
end

function index = closest(U_array,velocity)
    min = 1000;
    index = [0,0];
    for i = 1:length(U_array)
        for j = 1:length(U_array)
            min_temp = abs(velocity-U_array(i,j));
            if min_temp < min
                index(1,1) = i; index(1,2) = j;
            end
        end
    end
end

function [psi] = psif(r,lambda,theta,A,B,C,D)
    psi = (r^lambda).*(A*cos(lambda*theta)+B*sin(lambda*theta) + ...
              C*cos((lambda-2)*theta)+D*sin((lambda-2)*theta));
end

function [u_theta] = u_thetaf(r,lambda,theta,A,B,C,D)
    u_theta = -lambda*r^(lambda-1)*(A*cos(lambda*theta)+B*sin(lambda*theta) + ...
              C*(lambda-2)*cos((lambda-2)*theta)+D*(lambda-2)*sin((lambda-2)*theta));
end

function [u_r] = u_rf(r,lambda,theta,A,B,C,D)
    u_r = r^(lambda-1)*(-A*lambda*sin(lambda*theta)+B*lambda*cos(lambda*theta) - ...
              C*(lambda-2)*sin((lambda-2)*theta)+D*(lambda-2)*cos((lambda-2)*theta));
end

%[d_r u_r, d_r u_theta]
function [dru_r,dru_theta] = dru(r,lambda,theta,A,B,C,D)
    dru_r = (lambda-1)*r^(lambda-2)*(-A*lambda*sin(lambda*theta)+B*lambda*cos(lambda*theta)-...
        C*(lambda-2)*sin(theta*(lambda-2))+D*(lambda-2)*cos(theta*(lambda-2)));
    dru_theta = -lambda*(lambda-1)*r^(lambda-2)*(A*cos(lambda*theta)+B*sin(lambda*theta)+...
        C*cos(theta*(lambda-2))+D*sin(theta*(lambda-2)));
end

%[d^2_rr u_r, d^2_rr u_theta]
function [d2ru_r,d2ru_theta] = d2ru(r,lambda,theta,A,B,C,D)
    d2ru_r = (lambda-1)*(lambda-2)*r^(lambda-3)*(-A*lambda*sin(lambda*theta)+B*lambda*cos(lambda*theta)-...
        C*(lambda-2)*sin(theta*(lambda-2))+D*(lambda-2)*cos(theta*(lambda-2)));
    d2ru_theta = -lambda*(lambda-1)*(lambda-2)*r^(lambda-3)*(A*cos(lambda*theta)+B*sin(lambda*theta)+...
        C*cos(theta*(lambda-2))+D*sin(theta*(lambda-2)));
end

%[d^2_theta,theta u_r, d^2_theta,theta u_theta]
function [d2tu_r,d2tu_theta] = d2tu(r,lambda,theta,A,B,C,D)
    d2tu_r = r^(lambda-1)*(A*lambda^3*sin(lambda*theta)-B*lambda^3*(cos(lambda*theta))+...
        C*(lambda-2)^3*sin(theta*(lambda-2))-D*(lambda-2)^3*cos(theta*(lambda-2)));
    d2tu_theta = -lambda*r^(lambda-1)*(-A*lambda^2*cos(lambda*theta)-B*lambda^2*sin(lambda*theta)+...
        C*(lambda-2)^2*sin(theta*(lambda-2))-D*(lambda-2)^2*sin(theta*(lambda-2)));
end

%find the closest value based on r and theta, the corresponding velocity
%field
function [rindex,thetaindex] = closestField(r_array,theta_array,r_temp,theta_temp)
    [~, rindex] = min(abs(r_array-r_temp)); [~,thetaindex] = min(abs(theta_array-theta_temp));
end
