function m = leastsq_plane(x,y,z)

% m = leastsq_plane(x,y,z)
% This m-file performs the least squares linear regressio to find the plane
% that fits the data z in the least-square sense.
% 
% Input
%   x = N by 1 vector of x-coordinate
%   y = N by 1 vector of y-coordinate
%   z = N by 1 vector of z values
% 
% Output
%   m = 2 by 1 vecctor of coefficients a and b

% Written by Atsuhiro Muto
% Dept. of Earth & Env. Sci., Temple Univ.
% amuto@temple.edu
% Last updated March 10, 2017

G = [x y ones(size(x))];

m = G\z;
