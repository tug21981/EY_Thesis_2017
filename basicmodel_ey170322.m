%% RESIDUAL GRAVITY AND SNOW EFFECT MODEL
% Creation of the residual gravity using "SB_2016"
% "km10PRISM" is a grided DEM up to 10 km away from Wolverine Glacier, created in ArcGIS, pulled DEM values at 100 x 100 m spacing

% "km10PRISM" contains:
% column 1: easting
% column 2: northing
% column 3: elevation

% because of the way I wrote this code, it can be modified easily for the snow effect model, you must uncomment and replace the sections
% where it says "%this changes when using SWED values for the snow effect model" the "newSWED10km" is the values created with the MSWED
% models discussed in the thesis. It uses the same 100 x 100 m grid over the 10 km area surrounding Wolverine Glacier

% "newSWED10km" contains:
% column 1: easting
% column 2: northing
% column 3: modeled snow water equivalent depth for seasonal snow

% this section produces "gravityRESIDUAL" or "gravitySWE" depending on the use and contains:
% column 1: gravity signal, in mGal
% the easting and northing match that of "ObsMay2" for plotting purposes


tic
%% have to run USGS_calibrate first to get SB_2016

% removal of extraneous values (this was a product of the DEM size in ArcGIS)
load km10PRISM
km10PRISM(1,:) = [];
remove = find(km10PRISM(:,3) == -9999) ;
km10PRISM(remove,:) = [];

% grid of depth of the Wolverine prisms 

load ObsMay2
ObsMay2(1,:) = [];
remove = find(ObsMay2(:,3) == -9999) ;
ObsMay2(remove,:) = [];

L = length(ObsMay2(:,1));


%% using calibrated LiDAR data to create depth of prisms
load SB_2016


%% creating prisms around the center point

X = km10PRISM(:,1);
Y = km10PRISM(:,2);
Z = km10PRISM(:,3); 

% [X Y Z] = xyz2grid(X,Y,Z);
% mesh(X, Y, Z) 

dx = 100;
dy = dx;

% corner of prisms
Xbl = X-dx/2;         % x-coordinates of bottom left corner
Ybl = Y-dy/2;         % y-coordinates of bottom left corner
Xtr = X+dx/2;         % x-coordinates of top right corner
Ytr = Y+dy/2;         % y-coordinates of top right corner


S = size(X);
ny = S(1);
nx = S(2);

% positions of prisms
Pp = [reshape(Xbl',ny*nx,1) reshape(Xtr',ny*nx,1)...
        reshape(Ybl',ny*nx,1) reshape(Ytr',ny*nx,1)];

%% you have to resize the prisms so that they are in the right matrices for the code.

gravityRESIDUAL = nan(L,1);                                    %this changes to "gravitySWE" when using SWED values for the snow effect model

m = nan(ny*nx,2);

% uncomment if using outside SWEDs:                           %this changes when using SWED values for the snow effect model
% load newSWED10km
% Zb = km10PRISM(:,3) - newSWED10km(:,3);

% Depths of interfaces of density prisms for the entire 10 km area
prism = [reshape(Z',ny*nx,1) (reshape(Z',ny*nx,1)-(1e-10))];  %reshape(Zb',ny*nx,1)]; %this changes when using SWED values for the snow effect model
                                                              %reshape(Z',ny*nx,1)-(1e-10))]; %use this when creating residual gravity (creates negligible prisms depth for all of the prisms OFF of the glacier)

% changing the heights of the density prisms on the surface of the glacier

for j = 1:L
    replace = find(X == ObsMay2(j,1) & Y == ObsMay2(j,2));
    prism(replace,1) = ObsMay2(j,3);
    prism(replace,2) = ObsMay2(j,3) + SB_2016(j,3); %1e-10;    %this changes when using SWED values for the snow effect model (creates a negligible prism depth for the prisms directly under the surface of the glacier)
                                                    %SB_2016(j,3); %use this when finding residual gravity
    Z(replace) = ObsMay2(j,3);                                 %for a better resolution DEM starting height of the calculation point

end
 
 rho = -1000; %water density

for i = 1:L
    % XY corrdinates for the calculation points
    Pobs = [ObsMay2(i,1) ObsMay2(i,2)];

    % elevation at Pobs
    zPobs = ObsMay2(i,3);
   
    % elevations at Pp relative to zPobs
    % Elevation above Pobs is negative so the gravity effect from terrain above
    % Pobs will be negative, with z towards the center of Earth.
        
    m(:,1) = zPobs - prism(:,1);
    m(:,2) = zPobs - prism(:,2);
    

    gravityRESIDUAL(i,1) = g3dplouff(Pobs,Pp,m,rho);          %this changes to "gravitySWE" when using SWED values for the snow effect model

end

gravityRESIDUAL = gravityRESIDUAL*1e5;

save gravityRESIDUAL gravityRESIDUAL

%% plotting gravity

figure
[x,y,z] = xyz2grid(ObsMay2(:,1),ObsMay2(:,2),gravityRESIDUAL);
contourf(x,y,z,100,'LineStyle','none')
axis equal

%%
toc
