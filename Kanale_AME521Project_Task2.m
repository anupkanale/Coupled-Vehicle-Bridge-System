%-----------------------------------------------------%
% Program for solving Task 2 is AME 521 final project  %
% Anup Kanale, November 2016                          %
%-----------------------------------------------------%

clear all;
close all; clc;

% Is there Impulse?
impulse = 0; % 0- no, 1- yes
I0 = 1200;

%% Define some parameters
%--------------------------
a=4.5;
M = 2.5e3;
I = 3.2e2;
rho = 800;
EI = 7e8;
L = 400;
v0 = 18;

nTimeSteps = 12000;
tTot = (L-a)/v0;
dt=(L-a)/(nTimeSteps*v0);
time = 0:dt:tTot;
nBeam = 10000;
xBeam = linspace(0,L,nBeam); % The number of points on the beam is same as the number of time steps

nPhiTerms = 50; % number of terms of admissible function

%Intial conditions
y10 = 0.02;
yDot10 = 0.6;
y20 = 0.05;
yDot20 = -0.3;
z(:,1) =[zeros(nPhiTerms,1); y10; y20; zeros(nPhiTerms,1); yDot10; yDot20] ;

%% Mass Matrix
%----------------
Mb = zeros(nPhiTerms);
Kb = zeros(nPhiTerms);
for ii=1:nPhiTerms
    for jj=1:nPhiTerms
        fun1 = (1-cos(2*ii*pi*xBeam/L)).*(1-cos(2*jj*pi*xBeam/L));
        Mb(ii,jj) = trapz(xBeam, rho*fun1);
        
        fun2 = cos(2*ii*pi*xBeam/L).*cos(2*jj*pi*xBeam/L);
        Kb(ii,jj) = trapz(xBeam, EI*(2*pi/L)^4*ii^2*jj^2*fun2);
    end;
end;

mMat = zeros(nPhiTerms+2);
mMat(1:nPhiTerms, 1:nPhiTerms) = Mb;
mMat(nPhiTerms+1,nPhiTerms+1) = M/4+I/a^2;
mMat(nPhiTerms+1,nPhiTerms+2) = M/4-I/a^2;
mMat(nPhiTerms+2,nPhiTerms+2) = M/4+I/a^2;
mMat(nPhiTerms+2,nPhiTerms+1) = M/4-I/a^2;

%% Solution using Runge-Kutta 4th Order
%-------------------------------------------
for kk=1:nTimeSteps
    [A1] = getA(time(kk), nPhiTerms, mMat, Kb);
    [A2] = getA(time(kk)+dt/2, nPhiTerms, mMat, Kb);
    [A3] = getA(time(kk)+dt, nPhiTerms, mMat,Kb);
    dMat = zeros((nPhiTerms+2)*2,1); 
    if kk==nTimeSteps/2
        if impulse==1
            b = [zeros(1,nPhiTerms+2) 0 I0];
            dMat = [zeros(nPhiTerms+2,1); -inv(mMat)*b'];
        end;
    end;
    
    % 4th Order Runge-Kutta
    f1 = A1* z(:,kk);
    f2 = A2* ( z(:,kk) + dt/2*f1 );
    f3 = A2* ( z(:,kk) + dt/2*f2 );
    f4 = A3* ( z(:,kk) + dt*f3 );
    z(:,kk+1) = z(:,kk) + dt/6*(f1+2*f2+2*f3+f4) + dMat;
end;

figure()
plot(time,z(nPhiTerms+1,:),'b');
title('Vibration response at the rear end of the car');
xlabel('time', 'FontSize', 12); ylabel('y_1', 'FontSize', 12);
xlim([0 21.97]);

figure()
plot(time,z(nPhiTerms+2,:),'r');
title('Vibration response at the front end of the car');
xlabel('time', 'FontSize', 12); ylabel('y_2', 'FontSize', 12);
xlim([0 21.97]);

figure()
plotTimeIndex = [1, nTimeSteps/6, 2*nTimeSteps/6, 3*nTimeSteps/6, 4*nTimeSteps/6, 5*nTimeSteps/6, nTimeSteps];
phiMat = zeros(1,nPhiTerms);
beamDisp = zeros(nTimeSteps,nBeam);
for kk=1:10:nTimeSteps
    for xi=1:length(xBeam)
        for ii=1:nPhiTerms
            phiMat(ii) = 1-cos(2*pi*ii*xBeam(xi)/L);
        end;
        beamDisp(kk,xi) = phiMat*z(1:nPhiTerms,kk);
    end;
    plot(xBeam,beamDisp(kk,:), 18*time(kk),0, 'ro', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    ylim([-9e-3 9e-3]);
    drawnow;
end;