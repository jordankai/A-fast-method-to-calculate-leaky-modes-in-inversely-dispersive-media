function y=Re_Haskell_Rayleigh(k,om,thk,dns,vp,vs,flag)


%% Reference: 1 Tianyun Liu, The computations of reflection coefficients of multilayer structure based
%%              on the reformulation of Thomson-Haskell method
%%            2 Jian Chen et.al, Efficient Reformulation of the Thomson–Haskell Method for
%%              Computation of Surface Waves in Layered Half-Space
%%            3 P. W. Buchen and R. Ben-Hador, Free-mode surface-wave computations

%**************************************************************************
%    Written and Copyleft by Zhang Kai,IGGE, zhangkai@igge.cn
%    Version 1, 2020/4/11
%    Version 2, 2022/11/25
%    You may use and modify this code, provided you acknowledge the source
%**************************************************************************

n=length(vs);

% epsilon = 0.0001;
% while any(abs(om/k-vs)<epsilon) || any(abs(om/k-vp)<epsilon)
%     k = k * (1+epsilon);
% end

rms = sqrt((om./vs).^2/(k^2)-1);
rmp = sqrt((om./vp).^2/(k^2)-1);

if nargin == 6
    if imag(rms(end)) < 0  %% normal modes
        rms(end)=-rms(end);
    end
    if imag(rmp(end)) < 0  %% normal modes
        rmp(end)=-rmp(end);
    end

elseif nargin == 7
    if strcmp(flag,'leaky')
        if imag(rms(end)) > 0   %% leaky modes
            rms(end)=-rms(end);
        end
        if imag(rmp(end)) < 0
            rmp(end)=-rmp(end);
        end

    elseif strcmp(flag,'leakyP')
        if imag(rms(end)) > 0   %% leakyP modes
            rms(end)=-rms(end);
        end
        if imag(rmp(end)) > 0
            rmp(end)=-rmp(end);
        end

    else
        error('Wrong input parameter names!');
    end
else
    error('Wrong input parameter number!');
end


rm=2*(k*vs/om).^2;
T31=dns.*rm.*rmp;
T32=dns.*(1-rm);
T34=dns.*rm.*rms;

Qp=(rm-1)./rmp;
Qs=(rm-1)./rms;
p1=1./dns;
%% half-space matrix
T1 =[1 -rms(n); rmp(n) 1];
T3 =  [-T31(n) T32(n);-T32(n) -T34(n)];
S = T3/(T1);
Edu=zeros(2,2);

for i=(n-1):-1:1
    
    Edu(1,1)=exp(1i*rmp(i)*k*thk(i));
    Edu(2,2)=exp(1i*rms(i)*k*thk(i));
    
    T1 = [1 -rms(i); rmp(i) 1];
    T2=  [1  rms(i);-rmp(i)  1];
    T3 = [-T31(i) T32(i);-T32(i) -T34(i)];
    T4=  [T31(i)  T32(i);-T32(i) T34(i)];
 
    
    Q_inv=0.5*[ rm(i) -Qp(i) -p1(i)/rmp(i) -p1(i);
                Qs(i)  rm(i)  p1(i)        -p1(i)/rms(i);
                rm(i)  Qp(i)  p1(i)/rmp(i) -p1(i);
               -Qs(i)  rm(i)  p1(i)         p1(i)/rms(i)];
    
    C = Q_inv*[eye(2);S];
    Sg = Edu*C(3:4,:)/C(1:2,:)*Edu;    
    S = (T3+T4*Sg)/(T1+T2*Sg);
    
end
y=(det(S));