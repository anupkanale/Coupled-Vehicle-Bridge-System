function [A] = getA(tim, nPhiTerms, mMat, Kb)
%% DefinPhiTermse parameters anPhiTermsd conPhiTermsstanPhiTermsts
%---------------------------------------
% Car parameters
a = 4.5;
k1 = 5e5;
k2 = 5e5;
c1 = 3.6e3;
c2 = 3.6e3;
v0 = 18;

% Bridge parameters
L = 400;

%% Initialise Matrices
%-------------------------------
phiMatvt = zeros(1,nPhiTerms);
phiMatvtpa = zeros(1,nPhiTerms);
for ii=1:nPhiTerms
    phiMatvt(ii) = 1-cos(2*ii*pi*v0*tim/L);
    phiMatvtpa(ii) = 1-cos(2*ii*pi*(v0*tim+a)/L);
end;

kMat = zeros(nPhiTerms+2);
cMat = zeros(nPhiTerms+2);

% StiffnPhiTermsess Matrix
kMat(1:nPhiTerms,1:nPhiTerms) = Kb + k1*(phiMatvt'*phiMatvt) + k2*(phiMatvtpa'*phiMatvtpa);
kMat(1:nPhiTerms,nPhiTerms+1) = -k1*phiMatvt';
kMat(1:nPhiTerms,nPhiTerms+2) = -k2*phiMatvtpa';
kMat(nPhiTerms+1,1:nPhiTerms) = -k1*phiMatvt;
kMat(nPhiTerms+2,1:nPhiTerms) = -k2*phiMatvtpa;
kMat(nPhiTerms+1,nPhiTerms+1) = k1;
kMat(nPhiTerms+2,nPhiTerms+2) = k2;

% DampinPhiTermsg Matrix
cMat(1:nPhiTerms,1:nPhiTerms) = c1*(phiMatvt'*phiMatvt) + c2*(phiMatvtpa'*phiMatvtpa);
cMat(1:nPhiTerms,nPhiTerms+1) = -c1*phiMatvt';
cMat(1:nPhiTerms,nPhiTerms+2) = -c2*phiMatvtpa';
cMat(nPhiTerms+1,1:nPhiTerms) = -c1*phiMatvt;
cMat(nPhiTerms+2,1:nPhiTerms) = -c2*phiMatvtpa;
cMat(nPhiTerms+1,nPhiTerms+1) = c1;
cMat(nPhiTerms+2,nPhiTerms+2) = c2;

A = [zeros(nPhiTerms+2), eye(nPhiTerms+2);...
    -inv(mMat)*kMat, -inv(mMat)*cMat];