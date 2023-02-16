function   [cr,cr_real,cr_imag] = leaky_fast(freq,thk,dns,cvs,cvp,Qs,Qp)

% This function calculates the modal phase velocities in ealstic and viscoelastic
% vertically heterogeneous medium using Muller search techniques.
% Coding by Kai 2022/12/29  naturekai@126.com
% cr=wavenumber for surface waves
% cr_real= surface wave velocities
% cr_imag= attenuation coefficient


NUMINC = 200; % number of subsections divided in the wavenumber domain, so we can
% search NUMINC times which cover the whole area as much as possible.

TOL=1e-10;   %  Tolerance

%% If some solutions are missing, increase the NUMINC (more dense grid)
%% If some solutions are possibly not the true surface wave eignvalues, try increase the TOL.

muller_iter_num=30;  % max iterations for muller method
max_mode=5;
v_int=0.1;%1;  %% 根之间的最小速度间隔

if nargin == 5
    nf=length(freq);
    % Convert all input parameters to column vectors
    thk = reshape(thk,length(thk),1);
    dns = reshape(dns,length(dns),1);
    cvp = reshape(cvp,length(cvp),1);
    cvs = reshape(cvs,length(cvs),1);
    freq = reshape(freq,nf,1);
    cvs= repmat(cvs,1,nf);
    cvp= repmat(cvp,1,nf);

elseif nargin == 7
    nf=length(freq);
    % Convert all input parameters to column vectors
    thk = reshape(thk,length(thk),1);
    dns = reshape(dns,length(dns),1);
    cvp = reshape(cvp,length(cvp),1);
    cvs = reshape(cvs,length(cvs),1);
    freq = reshape(freq,nf,1);
    %% various viscoelastic models
    [cvs, cvp] = model_simple(cvs,cvp,Qs,Qp,freq);  %% constant Q viscoealstic model

    % [cvs, cvp] = visco_model(cvs,cvp,freq); % Carcione model
    % [cvs, cvp] = model_KK(cvs,cvp,Qs,Qp,freq);  % Lai model
    % [cvs, cvp] = model_KD(cvs,cvp,Qs,Qp,freq);
    % [cvs, cvp] = model_KD_constant(cvs,cvp,Qs,Qp,freq);
    % [cvs, cvp] = model_KK_constant(cvs,cvp,Qs,Qp,freq);  % Lai model
    % [cvs, cvp, Rs, Rp] = constant_Q(cvs,cvp,Qs,Qp,freq,0.00001);

else
    error('Wrong input parameter number!');
end

% Loop through the frequencies
% Muller method

num_imag=-0.001; %% num_imag=0 is usually enough for both the elastic and viscoelastic media.

cr = zeros(nf,max_mode);
cr_real = zeros(nf,max_mode);
cr_imag=cr_real;

fbar = waitbar(0,'Please wait...');
for j = 1:nf
    str=strcat('Calculating over frequency:  ',num2str(j/nf*100),'%');
    waitbar(j/nf,fbar,str);

    om = 2*pi*freq(j);
    kvs=om./cvs(:,j);
    %kvp=om./cvp(:,j);

    kk_min=min(real(kvs));
    %kk_min=min(real(kvp));
    kk_max=max(real(kvs));
    kk=linspace(0.8*kk_min,1.2*kk_max,NUMINC);
    dk=kk(2)-kk(1);

    cvalue=zeros(max_mode,1);
    root_num=1;
    for ii=1:(NUMINC)

        Z0(1)=kk(ii)-dk/4+1i*num_imag;
        Z0(2)=kk(ii)+1i*num_imag;
        Z0(3)=kk(ii)+dk/4+1i*num_imag;

        %RES = muller(@(x) Re_Haskell_Rayleigh(x,om,thk,dns,cvp(:,j),cvs(:,j),'leaky'),Z0,muller_iter_num,TOL,TOL); %% the modified Thomson-Haskell method
        %   RES = muller(@(x) secular_improve(x,om,thk,dns,cvp(:,j),cvs(:,j),'leaky'),Z0,muller_iter_num,TOL,TOL);%% the generalized reflection/transmission method
         RES = muller(@(x) Fast_Delta(x,om,thk,dns,cvp(:,j),cvs(:,j),'leaky'),Z0,muller_iter_num,TOL,TOL); %% the fast delta method

        %% viscoelastic case and leaky case
        if   (RES~=0) && real(RES) <= 1.2*kk_max && all(abs(om/real(RES)-om./real(cvalue))>v_int) && real(RES) >= kk_min && (imag(RES))<=0 && abs(real(RES))>abs(imag(RES))
            cvalue(root_num,1) = RES;
            root_num = root_num+1;

        end

        %% Break out of loop of maxroots is reached
        if root_num >= max_mode+1
            break;
        end

    end

    [~,index] = sort(real(cvalue),'descend');  %% sort the solutions based on phase velocity
    cr(j,:)=cvalue(index);


end

close(fbar);
delete(fbar);

for ii=1:length(freq)
    cr_real(ii,:) = 2*pi*freq(ii)./real(cr(ii,:));
    cr_imag(ii,:) = imag(cr(ii,:));
end
cr_real(isinf(cr_real))=NaN;
cr_imag(cr_imag==0)=NaN;




