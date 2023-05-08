close all; clear all;
N = 1; alpha = 15*pi/180;
[lambda,Ak] = data(N);
r = linspace(1E-3,1,100); theta = linspace(-alpha,alpha,100);
h = r(3)-r(2); k = theta(3)-theta(2);
timeiter_array = [100,10,1];

%% determine streamline functions
%[psi1,xx1,yy1] = stream(r,theta,lambda,Ak,alpha,N);
[psi2,xx2,yy2] = stream(r*0.9,theta,lambda,Ak,alpha,N);
[psi3,xx3,yy3] = stream(r*0.4,theta,lambda,Ak,alpha,N);
[psi4,xx4,yy4] = stream(r*0.2,theta,lambda,Ak,alpha,N);
%% determine velocity field
[velr2,velt2] = velocity_field(r*0.9,theta,lambda,Ak,alpha,N);
[velr3,velt3] = velocity_field(r*0.4,theta,lambda,Ak,alpha,N);
[velr4,velt4] = velocity_field(r*0.2,theta,lambda,Ak,alpha,N);
%% remove values that are too small to plot in the velocity field
[rvelr2,rvelt2] = reduce(psi2,r,theta,velr2,velt2); velr2 = rvelr2; velt2 = rvelt2; 
[rvelr3,rvelt3] = reduce(psi3,r,theta,velr3,velt3); velr3 = rvelr3; velt3 = rvelt3; 
[rvelr4,rvelt4] = reduce(psi4,r,theta,velr4,velt4); velr4 = rvelr4; velt4 = rvelt4; 

figure(2)
hold on
for w2 = 1:3
    timeiter = timeiter_array(w2);
    [cr,ctheta] = RK4(timeiter,alpha,r,theta,velr4,velt4,velr3,velt3,velr2,velt2);
    pos = [cr;ctheta]'; lp = length(pos);
    cr = pos(:,1); ctheta=pos(:,2);
    cx = cr.*cos(ctheta);
    cy = cr.*sin(ctheta);
    %scatter(cx,cy,'filled')
    plot(cx,cy)
    %quiver(xx1,yy1,velr1+-0.047*velr2+0.6E-4*velr3+-0.7E-7*velr4,velt1+-0.047*velt2+0.6E-4*velt3+-0.7E-7*velt4,1)
    % quiver(xx2,yy2,velr2,velt2,1)
    % quiver(xx2,yy2,velr2+velr3,velt2+velt3,1)
    quiver(xx2,yy2,velr2+velr3+velr4,velt2+velt3+velt4,1)
    CX{w2} = cy; CY{w2} = cy;
end
hold off
figure(1)
hold on
%contour(xx1,yy1,psi1)
contour(xx2,yy2,-0.047*psi2)
contour(xx3,yy3,0.6E-4*psi3)
contour(xx4,yy4,-0.7E-7*psi4)
hold off

figure(2)
hold on
scatter(cx,cy,'filled')
%quiver(xx1,yy1,velr1+-0.047*velr2+0.6E-4*velr3+-0.7E-7*velr4,velt1+-0.047*velt2+0.6E-4*velt3+-0.7E-7*velt4,1)
% quiver(xx2,yy2,velr2,velt2,1)
% quiver(xx2,yy2,velr2+velr3,velt2+velt3,1)
quiver(xx2,yy2,velr2+velr3+velr4,velt2+velt3+velt4,1)
hold off

%% L2 norm with 10K (w2 = 1) time step benchmark case
e1 = CX{2}(1:200)-CX{1}(1:200);
e2 = CX{3}(1:100)-CX{1}(1:100);
EL1 = sqrt(sum(abs(e1).^2));
EL2 = sqrt(sum(abs(e2).^2));

function [psi,xx,yy] = stream(r,theta,lambda,Ak,alpha,N)
    for count = 1:N
        for j = 1:length(r)
            for k = 1:length(theta)
                if k == 1 || k == length(theta) || j == 1 || j == length(r)
                    psi_temp(j,k) = 0;
                else
                    psi_temp(j,k) = (Ak(count,1)+1i*Ak(count,2))*r(j)^(lambda(count,1)+1i*lambda(count,2))*...
                        ((cos((lambda(count,1)+1i*lambda(count,2))*theta(k))/cos((lambda(count,1)+1i*lambda(count,2))*alpha))-...
                        (cos(((lambda(count,1)+1i*lambda(count,2))-2)*theta(k))/cos(((lambda(count,1)+1i*lambda(count,2))-2)*alpha)));
                end
                if count == 1
                    xx(j,k) = r(j)*cos(theta(k));
                    yy(j,k) = r(j)*sin(theta(k));
                end
            end
        end
        if count == 1
            psi = real(psi_temp);
            psi_temp = zeros(size(psi_temp));
        else
            psi = psi+real(psi_temp);
        end
    end
end
function [velr,velt] = num_vel_field(psi,r,theta,h,k)
    for m = 1:length(r)
        for n = 1:length(theta)
            if n == 1 || n == length(theta)
                velr(m,n) = 0; velt(m,n) = 0; 
                if n == 1
                    if m == 100
                        velr = 0; velt = 0;
                    else
                        velr(m,n+1) = 3/(4*k)*(psi(m,n+1)-psi(m,n)); velt(m,n+1) = -3/(4*h)*(psi(m+1,n)-psi(m,n));
                    end
                end
            elseif  m == 1 || m == length(r)
                velr(m,n) = 0; velt(m,n) = 0; 
                if m == 1 && n ~= 100
                    velr(m+1,n) = 0; velt(m+1,n) = 0;
                end
                if n == 100
                    velr = 0; velt = 0;
                end
            end
        end
    end
    for m2 = 1:length(r)
        for n2 = 1:length(theta) 
            if n2 ~= 1 && n2 ~= length(theta) && m2 ~= 1 && m2 ~= length(r)
                velr(m2,n2) = 3/(4*k)*(psi(m2,n2+1)-psi(m2,n2-1))-0.25*(velr(m2,n2+1)-velr(m2,n2-1));
                velt(m2,n2) = -3/(4*h)*(psi(m2+1,n2)-psi(m2-1,n2))-0.25*(velt(m2+1,n2)-velt(m2-1,n2));
            end
        end
    end
end
function [velr,velt] = velocity_field(r,theta,lambda,Ak,alpha,N)
    for count = 1:N
        strength = intensity(Ak,count); Ak(count,1) = strength*Ak(count,1); Ak(count,2) = strength*Ak(count,2);
        for j = 1:length(r)
            for k = 1:length(theta)
                if j == 1 || j == length(r) || k == 1 || k == length(theta)
                    velr_temp(j,k) = 0; velt_temp(j,k) = 0;
                else
                    velr_temp(j,k) = (Ak(count,1)+1i*Ak(count,2))*r(j)^(lambda(count,1)+1i*lambda(count,2)-1)*...
                        (-lambda(count,1)-1i*lambda(count,2))*sin((lambda(count,1)+1i*lambda(count,2)*theta(k)))/cos((lambda(count,1)+1i*lambda(count,2))*alpha)+...
                        (lambda(count,1)+1i*lambda(count,2)-2)*sin(((lambda(count,1)+1i*lambda(count,2)-2)*theta(k)))/cos((lambda(count,1)+1i*lambda(count,2)-2)*alpha);
                    velt_temp(j,k) = -(Ak(count,1)+1i*Ak(count,2))*(-2 + (lambda(count,1)+1i*lambda(count,2)))*sec(alpha*(-2+lambda(count,1)+1i*lambda(count,2))*sin((-2 + (lambda(count,1)+1i*lambda(count,2)))\theta(k)))-... 
                            (lambda(count,1)+1i*lambda(count,2))*r(j)^(lambda(count,1)+1i*lambda(count,2))*sec(alpha*(lambda(count,1)+1i*lambda(count,2)))*sin((lambda(count,1)+1i*lambda(count,2))\theta(k));
                end
                if count == 1
                    xx(j,k) = r(j)*cos(theta(k));
                    yy(j,k) = r(j)*sin(theta(k));
                end
            end
        end
        if count == 1
            velr = real(velr_temp); velt = real(velt_temp);
        else
            velr = velr+real(velr_temp); 
            velt = velt+real(velt_temp);
        end
    end
end

function [rvelr,rvelt] = reduce(psi,r,theta,velr,velt)
    sizer = length(velr); sizet = length(velt);
    for d = 1:sizer
        for e = 1:sizet
            if abs(velr(d,e)) < 0.001 && abs(velr(d,e)) > 0
                rvelr(d,e) = 0;
            else
                rvelr(d,e) = velr(d,e);
            end
            if abs(velt(d,e)) < 0.001 && abs(velt(d,e)) > 0
                rvelt(d,e) = 0;
            else
                rvelt(d,e) = velt(d,e);
            end
        end
    end
end

function [cr,ctheta] = RK4(n,alpha,r,theta,Ur1,Utheta1,Ur2,Utheta2,Ur3,Utheta3)
    dt = linspace(0.01,1,n*100);
    k = dt(3)-dt(2);
     %% RK-4 Implementation 
    %initialize in middle of the array
    %indr = 45; indtheta = 56;
    indr = 55; indtheta = 68;
    cr(1) = r(indr); ctheta(1) = theta(indtheta); 
    if cr(1) >= 0 && cr(1) <= 0.2
        ur = Ur1; utheta = Utheta1;
    elseif cr(1) > 0.2 && cr(1) <= 0.4
        ur = Ur2+Ur1; utheta = Utheta2+Utheta1;
    else
        ur = Ur3+Ur2+Ur1; utheta = Utheta3+Utheta2+Utheta1;
    end
    k1r = dt(1)*ur(indr,indtheta); newr = cr(1)+k1r; k1theta = dt(1)*utheta(indr,indtheta); newtheta = ctheta(1)+k1theta;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k2r = dt(1)*ur(indr,indtheta); newr = cr(1)+k1r/2; k2theta = dt(1)*utheta(indr,indtheta); newtheta = ctheta(1)+k1theta;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k3r = dt(1)*ur(indr,indtheta); newr = cr(1)+k2r/2; k3theta = dt(1)*utheta(indr,indtheta); newtheta = ctheta(1)+k2theta/2;
    [indr,indtheta] = closestField(r,theta,newr,newtheta); 
    k4r = dt(1)*ur(indr,indtheta); newr = cr(1)+k3r/2; k4theta = dt(1)*utheta(indr,indtheta); newtheta = ctheta(1)+k3theta/2;
    for t = 2:length(dt)
        cr(t) = cr(t-1)+(1/6)*(k1r+2*k2r+2*k3r+2*k4r);
        ctheta(t) = ctheta(t-1)+(1/6)*(k1theta+2*k2theta+2*k3theta+2*k4theta);
        %Stop particle if it hits a wall
        if ctheta(t) < -alpha || ctheta(t) > alpha
            cr(t) = cr(t-1)+(1/24)*(k1r+2*k2r+2*k3r+2*k4r);
            ctheta(t) = ctheta(t-1)+(1/24)*(k1theta+2*k2theta+2*k3theta+2*k4theta);
        else
            % use function to find the closest r index and theta index
            [indr,indtheta] = closestField(r,theta,cr(t),ctheta(t)); 
            % find the velocity field at that r and theta
            k1r = dt(t)*ur(indr,indtheta); newr = cr(t)+k1r/2; k1theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k1theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k2r = dt(t)*ur(indr,indtheta); newr = cr(t)+k1r/2; k2theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k1theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k3r = dt(t)*ur(indr,indtheta); newr = cr(t)+k2r/2; k3theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k2theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k4r = dt(t)*ur(indr,indtheta); newr = cr(t)+k3r/2; k4theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k3theta/2;
        end
        if cr(t) < 0 || cr(t) > 1
            cr(t) = cr(t-1)+(1/24)*(k1r+2*k2r+2*k3r+2*k4r);
            ctheta(t) = ctheta(t-1)+(1/24)*(k1theta+2*k2theta+2*k3theta+2*k4theta);
        else
            % use function to find the closest r index and theta index
            [indr,indtheta] = closestField(r,theta,cr(t),ctheta(t)); 
            % find the velocity field at that r and theta
            k1r = dt(t)*ur(indr,indtheta); newr = cr(t)+k1r/2; k1theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k1theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k2r = dt(t)*ur(indr,indtheta); newr = cr(t)+k1r/2; k2theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k1theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k3r = dt(t)*ur(indr,indtheta); newr = cr(t)+k2r/2; k3theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k2theta/2;
            [indr,indtheta] = closestField(r,theta,newr,newtheta); 
            k4r = dt(t)*ur(indr,indtheta); newr = cr(t)+k3r/2; k4theta = dt(t)*utheta(indr,indtheta); newtheta = ctheta(t)+k3theta/2;
        end
        
    end
end

%find the closest value based on r and theta, the corresponding velocity
%field
function [rindex,thetaindex] = closestField(r_array,theta_array,r_temp,theta_temp)
    [~, rindex] = min(abs(r_array-r_temp)); [~,thetaindex] = min(abs(theta_array-theta_temp));
end

function strg_ratio = intensity(Ak,k)
    pi = Ak(k,1); p1 = Ak(1,1); q1 = Ak(1,2);
    size_ratio = exp(pi/q1);      % Moffatt (3.11a)
    strg_ratio = exp(pi*p1/q1);   % Moffatt (3.12)
end

function [rpt,thetpt] = reflect(P,m)
    p = P(1);q = P(2); a = 1; b = m; 
    rpt = ((p*(a^2-b^2))-(2*b*(a*q)))/(a^2+b^2);
    thetpt = (q*(b^2-a^2)-(2*a*b*p))/(a^2+b^2);
end