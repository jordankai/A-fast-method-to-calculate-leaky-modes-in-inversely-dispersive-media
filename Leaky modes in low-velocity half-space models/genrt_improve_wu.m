function [Td,Rd,e11,e12,e21,e22,du,mu2,nus,nup] = genrt_improve_wu(thk,dns,cvp,cvs,om,k,flag)

% This function calculates the E and Lambda matrices (up-going and
% down-going matrices) for the P-SV case. Note that a separate function,
% updown, is provided for calculating the Lambda matrices for use in
% determining the displacement-stress vectors.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

%% wu 2016

cvs2 = cvs.^2; cvp2 = cvp.^2;
mu = dns.*cvs2;
mu2 = mu/mean(mu);

nn=length(cvs);
e11 = zeros(2,2,nn);
e12 = zeros(2,2,nn);
e21 = zeros(2,2,nn);
e22 = zeros(2,2,nn);
du = zeros(2,2,nn-1);

k2 = k^2; om2 = om^2;

ks2 = om2./cvs2;
gamas=sqrt(k2-ks2);

if  strcmp(flag,'normal')

    if real(gamas(end))<0   %% normal modes
        gamas(end)=-gamas(end);
    end
elseif strcmp(flag,'leaky')

    if real(gamas(end))>0   %% leaky modes
        gamas(end)=-gamas(end);
    end
else
    error('Wrong input parameters!');
end

nus=gamas/k;

kp2 = om2./cvp2;
gamap=sqrt(k2-kp2);
if real(gamap(end))<0
    gamap(end)=-gamap(end);
end
nup=gamap/k;

x_wu=1-0.5*ks2/k2;

e11(1,1,:) = ones(nn,1);
e11(1,2,:) = nus;
e12(1,1,:) =  ones(nn,1);
e12(1,2,:) = nus;

e11(2,1,:) = nup;
e11(2,2,:) =  ones(nn,1);
e12(2,1,:) = -nup;
e12(2,2,:) = -1*ones(1,nn);

e21(1,1,:) = -mu2.*nup;
e21(1,2,:) = -mu2.*x_wu;
e22(1,1,:) = -e21(1,1,:);
e22(1,2,:) = -e21(1,2,:);

e21(2,1,:) = e21(1,2,:);
e21(2,2,:) = -mu2.*nus;
e22(2,1,:) = e21(2,1,:);
e22(2,2,:) = e21(2,2,:);


% 矩阵逆
mofa = 1./(2*mu2.*(1-x_wu));

e11_inv(1,1,:) = mu2.*mofa;
e11_inv(1,2,:) = -mu2.*x_wu./nup.*mofa;
e12_inv(1,1,:) = -1./nup.*mofa;
e12_inv(1,2,:) = ones(nn,1).*mofa;

e11_inv(2,1,:) =-mu2.*x_wu./nus.*mofa;
e11_inv(2,2,:) = mu2.*mofa;
e12_inv(2,1,:) =  ones(nn,1).*mofa;
e12_inv(2,2,:) = -1./nus.*mofa;

e21_inv(1,1,:) = mu2.*mofa;
e21_inv(1,2,:) = - e11_inv(1,2,:);
e22_inv(1,1,:) = -e12_inv(1,1,:);
e22_inv(1,2,:) = ones(nn,1).*mofa;

e21_inv(2,1,:) = e11_inv(2,1,:);
e21_inv(2,2,:) = -mu2.*mofa;
e22_inv(2,1,:) = -1*ones(nn,1).*mofa;
e22_inv(2,2,:) = e12_inv(2,2,:);


du(1,1,:) = exp(-gamap(1:length(thk)).*thk);
du(2,2,:) = exp(-gamas(1:length(thk)).*thk);

N = length(thk);

Td = zeros(2,2,N);
Rd = zeros(2,2,N+1);

% Calculate the Td and Rd matrices for the Nth layer
du(1:2,1:2,N+1) = 1; % 任意值都行
Rd(1:2,1:2,N+1) = 0;

% Loop through the first N-1 layers in reverse order
for j = N:-1:1
    A = [e11_inv(:,:,j) e12_inv(:,:,j); e21_inv(:,:,j) e22_inv(:,:,j)];
    B = [e11(:,:,j+1) e12(:,:,j+1); e21(:,:,j+1) e22(:,:,j+1)];
    L = A*B;
    L11 = L(1:2,1:2);
    L12 = L(1:2,3:4);
    L21 = L(3:4,1:2);  L22 = L(3:4,3:4);

    Td(:,:,j) = (L11+L12*du(:,:,j+1)*Rd(:,:,j+1))\du(:,:,j);
    Rd(:,:,j) = L21*Td(:,:,j) + L22*du(:,:,j+1)*Rd(:,:,j+1)*Td(:,:,j);
end



