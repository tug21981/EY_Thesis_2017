%% VERTICAL GRAVITY GRADIENT
% The next section is in parts to create the VGG and requires multiple DEM grids. They are 100 x 100 m grid up to 10 km away from the
% center of the glacier. 200 x 200 m grid up to 50 km away, 1000 x 1000 m grid up to 100 km away, and 2000 x 2000 m grid up to 200 km away.
% The larger grids have the center cut out to not reproduce the gravity signal each time (instead creating a single signal)

% Each section uses the same code to produce the gravity anomaly from each section and each grided DEM value numeric matrix contains:
% column 1: easting
% column 2: northing
% column 3: elevation

% Each section produces a file "g##km2695" that shows the gravity signal in mGal at each calculation point at various heights above the ground:
% column 1: at 0.0 m
% column 2: at 0.1 m
% column 3: at 0.2 m
% column 4: at 0.3 m
% column 5: at 0.4 m
% column 6: at 0.5 m
% column 7: at 0.6 m
% column 8: at 0.7 m
% column 9: at 0.8 m
% column 10: at 0.9 m
% column 11: at 1.0 m

% the easting and northing match that of "ObsMay2" for plotting purposes

% The final VGG produced "vgg2695" shows the VGG at each height above the ground up to 1m:
% column 1: at 0.05 m
% column 2: at 0.15 m
% column 3: at 0.25 m
% column 4: at 0.35 m
% column 5: at 0.45 m
% column 6: at 0.55 m
% column 7: at 0.65 m
% column 8: at 0.75 m
% column 9: at 0.85 m
% column 10: at 0.95 m

% the easting and northing match that of "ObsMay2" for plotting purposes


% run all prism codes separately

g200km_170317
clear

g100km_170317
clear

g50km_170317
clear

g10km_170317
clear

%%

load g200km2695
load g100km2695
load g50km2695
load g10km2695

gtotal2695 = g200km2695+g100km2695+g50km2695+g10km2695;

dh = 0.1;

vgg2695 = -0.308+diff(gtotal2695,1,2)/dh;
dhp = 0.05:dh:0.95;

save vgg2695.mat vgg2695

%% plotting VGG

load ObsMay.mat
ObsMay(1,:) = [];

[X,Y,Z] = xyz2grid(ObsMay(:,1),ObsMay(:,2),vgg2695(:,1));

figure
contourf(X,Y,Z,100,'linestyle','none')

%% density sensitivity
% the previous code for VGG was modified by changing the rho in each
% prisms' code. The value saved then had a version of 

load vgg2750
vgg2750new = vgg2750 - vgg2695;

[X2,Y2,Z2] = xyz2grid(ObsMay(:,1),ObsMay(:,2),vgg2750new(:,1));

figure
contourf((X2),Y2,Z2,100,'linestyle','none')

load vgg2660

load vgg2695

vggdiff1 = vgg2695(:,1) - vgg2660(:,1);
vggdiff2 = vgg2695(:,1) - vgg2695(:,1);

vggdiff1(1760,:) = [];
vggdiff2(1760,:) = [];

a = mean(vggdiff1(:,1))
b = mean(vggdiff2(:,1))


%% FREE-AIR EFFECT
% to produce the FAE the VGG is multiplied by the LiDAR uncorrected data
% (normal height change without density correction)

% both "FAEvgg2695" and "FAEnormal" produce a single column of FAE values in mGal
% easting and northings for plotting were taken from "SB_2016"

load SB_2016
load vgg2695

SB_2016(1760,:) = [];
vgg2695(1760,:) = [];


% FAEnormal created for comparison to the modeled results of VGG using the normal VGG
FAEnormal = -0.3085 * SB_2016(:,3);

L = length(vgg2695);
FAEvgg2695 = nan(L,1);

for i = 1:L
    FAEvgg2695(i,1) = vgg2695(i,1) * SB_2016(i,3);
end

% plotting each FAE
figure
[x,y,z] = xyz2grid(SB_2016(:,1),SB_2016(:,2),FAEnormal);
contourf(x,y,z,100,'LineStyle','none')
axis equal

figure
[x,y,z] = xyz2grid(SB_2016(:,1),SB_2016(:,2),FAEvgg2695);
contourf(x,y,z,100,'LineStyle','none')
axis equal

