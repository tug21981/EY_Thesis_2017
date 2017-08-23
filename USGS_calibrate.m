%% CORRECTION OF LiDAR DATA TO USGS SUMMER MASS BALANCE

% This section runs the LiDAR data "ObsMay2" and "ObsSept" alongside the USGS stake summer balance ("USGSsumbal") to create a 
% planar least squares, so the LiDAR is corrected for density.

% The data is all in reference to Wolverine Glacier, AK.

% The LiDAR data was created using ARCGIS on a 100 x 100 m grid spacing over the glacier and pulled the DEM points onto the grid.
% The USGS summer balance data was collected at 23 stakes over the glacier at various locations.

% "ObsMay2" and "ObsSept" contain:
% column 1: easting
% column 2: northing
% column 3: elevation

% "USGSsumbal" contains:
% column 1: easting
% column 2: northing
% column 3: elevation
% column 4: summer mass balance


% Produces the numeric matrix of "SB_2016"

% "SB_2016" contains:
% column 1: easting
% column 2: northing
% column 3: corrected LiDAR summer mass balance


%% CORRECTION OF LiDAR DATA
%% Working with USGS data

load USGS_sumbal
remove = [5,17,21,22];
USGSsumbal(remove,:) =[];


% finding closest grid points to the USGS stake locations

N = length(USGSsumbal(:,1));

load ObsMay2
load ObsSept
ObsSept(1,:) = [];
rm = find(ObsSept(:,3) == -9999) ;
ObsSept(rm,:) = [];

ObsDiff(:,1) = ObsMay2(:,1);
ObsDiff(:,2) = ObsMay2(:,2);
ObsDiff(:,3) = -(ObsMay2(:,3) - ObsSept(:,3));

M = length(ObsMay2(:,1));

e = nan(N,1);
n = nan(N,1);
df = nan(N,1);


for i = 1:N
    
    for j = 1:M
        
        distance(j,1) = sqrt(( (USGSsumbal(i,1) - ObsMay2(j,1)) ^2) + ((USGSsumbal(i,2) - ObsMay2(j,2))^2));
    end
    
        [Q,I] = min(distance);
        e(i,1) = ObsDiff(I,1);
        n(i,1) = ObsDiff(I,2);
        df(i,1) = ObsDiff(I,3);
  
end


df(:,1) =  df(:,1)- USGSsumbal(:,4);

m = leastsq_plane(e,n,df);
sb = ObsDiff(:,3) - (m(1)*ObsDiff(:,1)+m(2)*ObsDiff(:,2)+m(3));
SB_2016 = [ObsDiff(:,1:2) sb];
save SB_2016 SB_2016

%% creation of smoothed LiDAR summer balance map

[X,Y,Z] = xyz2grid(SB_2016(:,1),SB_2016(:,2),SB_2016(:,3));

figure
contourf(X,Y,Z,100,'linestyle','none')
axis equal

%% plotting elevation vs summer balance

N = length(USGSsumbal(:,1));
M = length(SB_2016(:,1));

forplot = nan(N,2);


for i = 1:N
    
    for j = 1:M
        
        distance(j,1) = sqrt(( (USGSsumbal(i,1) - SB_2016(j,1)) ^2) + ((USGSsumbal(i,2) - SB_2016(j,2))^2));
    end    
        [Q,I] = min(distance);

        forplot(i,1) = SB_2016(I,3);
        forplot(i,2) = USGSsumbal(i,3);   
end

% plot of the USGS summer balance vs the LiDAR summer balance against elevation
% USGS data in red, LiDAR data in blue
figure
plot(USGSsumbal(:,3),USGSsumbal(:,4),'ro',forplot(:,2),forplot(:,1),'bo')