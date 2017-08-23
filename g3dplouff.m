function gp = g3dplouff(Pobs,Pp,zp,rho)

% gp = plouff(Pobs,Pp,zp,rho,dx,dy)
% This is a forward model of gravity due to rectangluar prisms with
% specified number of layers (L). For details see Plouff (1976) and Seber 
% et al. (2001).
% 
% Input
%   Pobs = N by 2 matrix of data points, x in 1st column, y in 2nd column
%   Pp = M by 2 matrix of lower left corner position of the prism, 
%        x in 1st column, y in 2nd column
%   zp = L+1 by 1 vector of depth of layer interfaces
%   rho = L by 1 vector of density of layers
% 
% Output
%   gp = N by 1 vector of gravity at each data point

% Written by Atsuhiro Muto
% Dept. of Earth & Env. Sci., Temple University
% amuto@temple.edu
% Last updated Jan. 24, 2017

% gravitational constant, (m^3/kg/sec)
G = 6.673e-11;

S = size(Pobs);
N = S(1);               % number of data points
S = size(Pp);
M = S(1);               % number of prisms
L = length(rho);        % number of layers

% xp = nan(M,1);
% yp = nan(M,1);
R = nan(M,4,L+1);
U = nan(M,4,L+1);
gl = nan(M,1);
gp = zeros(N,1);
                
% -------------------------------------------------------------------------
% loop over data points
for n = 1:N
    
    % calculate x-positions
    xp = [Pp(:,1)-Pobs(n,1) Pp(:,2)-Pobs(n,1)];
    % flip if obs. point is further from origin
    id = find(xp(:,1)<0 & xp(:,1)<xp(:,2) & xp(:,1).*xp(:,2)>=0);
    xp(id,:) = abs(fliplr(xp(id,:)));
    
    % calculate y-positions
    yp = [Pp(:,3)-Pobs(n,2) Pp(:,4)-Pobs(n,2)];
    % flip if obs. point is further from origin
    id = find(yp(:,1)<0 & yp(:,1)<yp(:,2) & yp(:,1).*yp(:,2)>=0);
    yp(id,:) = abs(fliplr(yp(id,:)));
    
    si = 1; sj = 1;     % initialize s for i and j

    % -----------------------------------------------------------------
    % loop over corners of the prism, s for only x and y are
    % incorporated here
    id = 1;
    for i=1:2   % x-points
        si = si*-1;
        for j=1:2   % y-points
            sj = sj*-1;
            s = si*sj;
            for k = 1:L+1   % z-points, each interface
                s = s*-1;
                R(:,id,k) = sqrt(xp(:,i).^2+yp(:,j).^2+zp(:,k).^2);
                neg = find(zp(:,k)<=0);
                A = atan(xp(:,i).*yp(:,j)./zp(:,k)./R(:,id,k));
                A(neg) = atan(-xp(neg,i).*yp(neg,j)./zp(neg,k)./R(neg,id,k));
                U(:,id,k) = s * (zp(:,k) .* A ...
                    - xp(:,i).*log(R(:,id,k)+yp(:,j))...
                    - yp(:,j).*log(R(:,id,k)+xp(:,i)));
            end
            id = id+1;
        end
    end
    % -----------------------------------------------------------------
    
    % gravity effect of each layer in the prism, s for z is
    % incorporated in this part
%     gl = nan(M,2);
    for l=1:L
        gl(:,l) = rho(l) * (sum(U(:,:,l),2) + sum(U(:,:,l+1),2));
    end
    
    % gravity effect of all layers in the prism
    gp(n) =  G * sum(sum(gl,1));
    
    % ---------------------------------------------------------------------

end
% -------------------------------------------------------------------------
