%-----------------------------------------------------%
% Program for solving Task 1 is AME 521 final project  %
% Anup Kanale, November 2016                          %
%-----------------------------------------------------%

clear all;
close all; clc;

%% Define parameters and constants
%---------------------------------------
% Engine assembly
me = 200;
k2 = 1.3*10^5;
c2 = 1.02*10^3;
d = 1.65;

%Car Body
Mb = 1650;
Ib = 2330;
k1 = 2.5*10^5;
c1 = 2.73*10^3;

% Tire-rim assembly
m = 75;
k = 2.5*10^5;
c = 1.48*10^3;
l1 = 1.4;
l2 = 1.3;

%% Define Mass, Stiffness and damping matrices
%-----------------------------------------------------
mMat = [m m Mb Ib me];
mMat = diag(mMat);

kMat(1,:) = [k+k1 0 -k1 -k1*l2 0];
kMat(2,:) = [0 k+k1 -k1 k1*l1 0];
kMat(3,:) = [-k1 -k1 2*k1+k2 k2*d+k1*l2-k1*l1 -k2];
kMat(4,:) = [-k1*l2 k1*l1 -k1*l1+k1*l2+k2*d k1*l1^2+k1*l2^2+k2*d^2 -k2*d];
kMat(5,:) = [0 0 -k2 -k2*d k2];

cMat(1,:) = [c+c1 0 -c1 -c1*l2 0];
cMat(2,:) = [0 c+c1 -c1 c1*l1 0];
cMat(3,:) = [-c1 -c1 2*c1+c2 c2*d+c1*l2-c1*l1 -c2];
cMat(4,:) = [-c1*l2 c1*l1 -c1*l1+c1*l2+c2*d c1*l1^2+c1*l2^2+c2*d^2 -c2*d];
cMat(5,:) = [0 0 -c2 -c2*d c2];

%% 1(b) Eigenvalues and Eigenvectors of damped System
%------------------------------------------------------------
A = [zeros(5) eye(5); -inv(mMat)*kMat -inv(mMat)*cMat];
[eVecDamp, eValDamp] = eig(A);

%% 1(c) Eigenvectors with no damping
%------------------------------------------
A = inv(mMat)*kMat;
[eVec, eVal] = eig(A);
ProdMat = zeros(5);  % Check for orthogonality of eigen vectors of undamped system
for ii=1:5
    for jj=1:5
        if ii~=jj
            ProdMat(ii,jj) = eVec(:,ii)'*mMat*eVec(:,jj);
            if ProdMat(ii,jj)<1e-10
                ProdMat(ii,jj) = 0;
            end;
        end;
    end;
end;

%% 1(d) Frequency Response
%---------------------------
Lr = 10;
y0 = 0.08;  
vc = linspace(0,50,50);
omega = vc/Lr;
qs = zeros(1,5); qc = zeros(1,5);
yss = zeros(10,50);

for ii=1:length(vc)
    t0 = (l1+l2)/vc(ii);
    qs(1) = k*y0;
    qc(1) = c*y0*omega(ii);
    qs(2) = c*y0*omega(ii)*sin(omega(ii)*t0)+k*y0*cos(omega(ii)*t0);
    qc(2) =  c*y0*omega(ii)*cos(omega(ii)*t0)-k*y0*sin(omega(ii)*t0);
    mat = [-omega(ii)^2*mMat+kMat -omega(ii)*cMat; omega(ii)*cMat -omega(ii)^2*mMat+kMat];
    yss(:,ii) = inv(mat)*[qs, qc]';
    
    mag(ii) = sqrt(yss(5,ii)^2 + yss(10,ii)^2);
    phase(ii) = atan(yss(10,ii)/yss(5,ii));
end;

plot(vc,mag, 'LineWidth', 2);
xlabel('v_c', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);

figure()
plot(vc,phase, 'r','LineWidth', 2);
xlabel('v_c', 'FontSize', 12);
ylabel('Phase', 'FontSize', 12);