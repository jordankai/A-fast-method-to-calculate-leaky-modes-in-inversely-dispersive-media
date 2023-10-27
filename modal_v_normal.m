function cr = modal_v_normal(freq,thk,dns,cvs,cvp,Qs,Qp)

% This function calculates the modal phase velocities in an elastic,
% vertically heterogeneous medium using search techniques.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Establish global parameters
%TOL=1e-3;
MAXROOT=5;

NUMINC=200;

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
    %[cvs, cvp] = model_simple(cvs,cvp,Qs,Qp,freq);  %% constant Q viscoealstic model

    % [cvs, cvp] = visco_model(cvs,cvp,freq); % Carcione model
    % [cvs, cvp] = model_KK(cvs,cvp,Qs,Qp,freq);  % Lai model
    % [cvs, cvp] = model_KD(cvs,cvp,Qs,Qp,freq);
    % [cvs, cvp] = model_KD_constant(cvs,cvp,Qs,Qp,freq);
    % [cvs, cvp] = model_KK_constant(cvs,cvp,Qs,Qp,freq);  % Lai model
    [cvs, cvp, Rs, Rp] = constant_Q(cvs,cvp,Qs,Qp,freq,0.00001);

else
    error('Wrong input parameter number!');
end



% Initialize a matrix to store modal phase velocities

cr = zeros(nf,MAXROOT);
fbar = waitbar(0,'Please wait...');

% Loop through the frequencies
for j = 1:nf


    str=strcat('Calculating over frequency:  ',num2str(j/nf*100),'%');
    waitbar(j/nf,fbar,str);

    if length(cvs(:,j)) == 1

        cr(j,1) = homogeneous_visco(cvp(:,j),cvs(:,j));
    else
        numroot = 0;
        om = 2*pi*freq(j);

        % Establish the search parameters

        vmin=0.8*min(real(cvs(:,j)));  vmax=max(real(cvs(:,j)));

        dv = (vmax - vmin)/NUMINC;

        % Establish the dispersion equation
        %fun=@(x) abs((Re_Haskell_Rayleigh(x,om,thk,dns,cvp(:,j),cvs(:,j))));  %% the modified Thomson-Haskell method

          fun=@(x) abs((Fast_Delta(x,om,thk,dns,cvp(:,j),cvs(:,j)))); %% the fast delta method


        %         fun=@(x) abs((secular_improve(x,om,thk,dns,cvp(:,j),cvs(:,j),'normal')));   %% the generalized reflection/transmission method


        k1 = vmin;
        f1 = fun(om/k1);
        k2 = vmin + dv;
        f2 = fun(om/k2);

        % Loop through the remaining points
        for m = 2:NUMINC-1
            k3 = vmin + m*dv;
            f3 = fun(om/k3);

            % Determine if a minimum is bracketed
            if (f2 < f1) && (f2 < f3)

                % Use golden search/parabolic interpolation to refine minimun
                [ktrial,ftrial] = fminbnd(fun,om/k3,om/k1,optimset('TolX',1e-8,'Display','off'));

                % Check to see if ktrial is a zero and different from the previous zero
                mmt=om/ktrial;
                if  (mmt<vmax) && (mmt>=vmin)
                    numroot = numroot + 1;
                    cr(j,numroot) = mmt;

                end
            end

            % Break out of loop of maxroots is reached
            if numroot == MAXROOT
                break;
            end

            k1 = k2; f1 = f2;
            k2 = k3; f2 = f3;

        end
    end
end

close(fbar);
% delete(fbar);
