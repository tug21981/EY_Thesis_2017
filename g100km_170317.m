tic

load newBATH100
% newBATH100(1,:) = [];
remove = find(newBATH100(:,3) == -9999);
newBATH100(remove,:) =[];

% create a find and replace for the elevations above sea level
load km100RIM
km100RIM(1,:) = [];
remove = find(km100RIM(:,3) == -9999);
km100RIM(remove,:) = [];
replace = find(km100RIM(:,3) <= 0);
km100RIM(replace,:) = newBATH100(replace,:);

dx = 1000;
dy = dx;

%% Above sea level

idasl = find(km100RIM(:,3)>=0);
ASL = km100RIM(idasl,:);

% corner of prisms
Xbl1 = ASL(:,1)-dx/2;         % x-coordinates of bottom left corner
Ybl1 = ASL(:,2)-dy/2;         % y-coordinates of bottom left corner
Xtr1 = ASL(:,1)+dx/2;         % x-coordinates of top right corner
Ytr1 = ASL(:,2)+dy/2;         % y-coordinates of top right corner

S = size(ASL(:,1));
ny1 = S(1);
nx1 = S(2);

% positions of prisms relative to Pobs
Pp1 = [reshape(Xbl1',ny1*nx1,1) reshape(Xtr1',ny1*nx1,1)...
      reshape(Ybl1',ny1*nx1,1) reshape(Ytr1',ny1*nx1,1)];

%% Below sea level

idbsl = find(km100RIM(:,3)<0);
BSL = km100RIM(idbsl,:);

% corner of prisms
Xbl2 = BSL(:,1)-dx/2;         % x-coordinates of bottom left corner
Ybl2 = BSL(:,2)-dy/2;         % y-coordinates of bottom left corner
Xtr2 = BSL(:,1)+dx/2;         % x-coordinates of top right corner
Ytr2 = BSL(:,2)+dy/2;         % y-coordinates of top right corner

S = size(BSL(:,1));
ny2 = S(1);
nx2 = S(2);

% positions of prisms relative to Pobs
Pp2 = [reshape(Xbl2',ny2*nx2,1) reshape(Xtr2',ny2*nx2,1)...
      reshape(Ybl2',ny2*nx2,1) reshape(Ytr2',ny2*nx2,1)];

%%
% Load observation points
load ObsMay2.mat
ObsMay2(1,:) = [];
N = length(ObsMay2);

rhoc = 2695;        % crustal density
rhoo = 1028-2695;   % ocean-water density, contrast from the crust

dh = 0.1;
h = 0:dh:1;
L = length(h);
g100km2695 = nan(1760,11);

%%
% Loop through observation points
for j = 1:N
    
    Pobs = [ObsMay2(j,1) ObsMay2(j,2)];
    Pobs = Pobs+0.001;
    
    % elevation at Pobs
    zPobs = ObsMay2(j,3);
    
    % elevations at Pp relative to zPobs
    % Elevation above Pobs is negative so the gravity effect from terrain above
    % Pobs will be negative, with z towards the center of Earth.
    zPp1 = zPobs-ASL(:,3);
    m01 = [ones(ny1*nx1,1)*1e-10 zPp1];
    m02 = [ones(ny2*nx2,1)*1e-10 BSL(:,3)];
    m02 = zPobs-m02;

    g = nan(1,L);
    
    % Loop through heights above the observation points
    for i=1:L
        m1 = m01+h(i);
        m1(:,1) = 1e-10;
        m2 = m02+h(i);
        
        % above sea level
        g1 = g3dplouff(Pobs,Pp1,m1,rhoc);
        
        % below sea level
        g2 = g3dplouff(Pobs,Pp2,m2,rhoo);
        
        g(i) = g1+g2;
    end
    
    g100km2695(j,:) = g*1e5;
    
end
toc

save g100km2695 g100km2695