%% Steady State Solution With Point Particle Swimmer
length_wall = 1; %in meters
n = 100;
alpha = pi/12; %half the angle of the wedge
lambda = 1:50; % different eddies
theta = linspace(-alpha,alpha,n);
dtheta = theta(3)-theta(2);
r = linspace(0,1,n); r(1) = 1E-3; %buffer to avoid singularity
dr = r(3)-r(2);
A = 1; B = 1; C = 1; D = 1;

%streamline function
psi = ones(length(r),length(theta));

%velocity function
% U = ones(length(r),length(theta)); % U(r,theta)
% U_r = zeros(length(r),length(theta)); % U(r,theta)
% U_theta = zeros(length(r),length(theta)); % U(r,theta)

xx = ones(length(r),length(theta)); 
yy = ones(length(r),length(theta)); 

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
            else
                U{i2+1,j} = [0,0]; %BC
            end
        end
    end
    U(1,:) = [];
    %U = U/norm(U,2); 
    psi = psi/norm(psi,2);
    % randomize an initial position for the particle
    p_ar(1,1) = r(90); p_ar(1,2) = theta(90);
    indr = round((p_ar(1,1)/dr)+1); 
    if p_ar(1,2) < 0
        indtheta = round((abs(p_ar(1,2))/dtheta)+1);
    else
        indtheta = round((length(theta)-(abs(p_ar(1,2))/dtheta))+1);
    end
    vel = U{indr,indtheta};
    for counti = 2:1000
        if counti == 2
            p_ar(counti,1) = vel(1)+p_ar(1,1); 
            p_ar(counti,2) = vel(2)+p_ar(1,2);
        else
            % indr = round((p_ar(counti-1,1)/dr)+1); 
            % if p_ar(counti-1,2) < 0
            %     deg_ind = abs(p_ar(1,2))/dtheta;
            %     indtheta = round(deg_ind+1);
            % elseif p_ar(counti-1,2) > 0
            %     deg_ind = abs(p_ar(1,2))/dtheta;
            %     indtheta = round(deg_ind+1);
            % end
            % vel = U{indr,indtheta};
            % p_ar(counti,1) = vel(1)+p_ar(counti-1,1); 
            % p_ar(counti,2) = vel(2)+p_ar(counti-1,2);
        end
    end
    figure(d)
    contour(xx,yy,psi)
    f = gcf;
    title_temp = strcat('eigenvalue ',int2str(d), 'stream function','.jpg');
    saveas(f,title_temp)
    title{j} = title_temp;
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
    u_theta = lambda*r^(lambda-1)*(A*cos(lambda*theta)+B*sin(lambda*theta) + ...
              C*cos((lambda-2)*theta)+D*sin((lambda-2)*theta));
end

function [u_r] = u_rf(r,lambda,theta,A,B,C,D)
    u_r = r^lambda*(-A*lambda*sin(lambda*theta)+B*lambda*cos(lambda*theta) - ...
              C*(lambda-2)*sin((lambda-2)*theta)+D*(lambda-2)*cos((lambda-2)*theta));
end