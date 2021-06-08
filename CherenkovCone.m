
%{
This script was used to analyze raw output from a GAMOS simulation of Cherenkov light 
emitted from a thin slab of dielectric and determine the dimensions of the cone of light 
incident on a plane some distance from the slab.
%}

%% Path

addpath(genpath('/Users/danielalexander/Documents/Analysis/Cherenkov Cone'));
addpath(genpath('/Users/danielalexander/Documents/Matlab Functions'));
cd '/Users/danielalexander/Documents/Analysis/Cherenkov Cone'

%% Rename CSVs

% Adds .csv to file names

datapath = '/Users/danielalexander/Documents/Analysis/Cherenkov Cone/Sim Data'; 

files = dir(datapath); % get all files

for i = 1:numel(files)
    if contains(files(i).name, 'inclusion')
        file = fullfile(datapath, files(i).name);
        movefile(file, strcat(file, '.csv')) % ad .csv
    end
    disp(num2str(100*i/numel(files)));
end

%% Scored Values

fields = {'InitialPosX', 'InitialPosY', 'InitialPosZ', 'InitialDirX', 'InitialDirY', 'InitialDirZ' ...
    'FinalPosX', 'FinalPosY', 'FinalPosZ', 'FinalDirX', 'FinalDirY', 'FinalDirZ', 'FinalTotalEnergy'};

%% Read Data

datapath = '/Users/danielalexander/Documents/Analysis/Cherenkov Cone/Sim Data'; 
ttds = tabularTextDatastore(datapath,'FileExtensions',{'.csv'});
ttds.VariableNames = fields;

ttds.SelectedVariableNames = {'InitialPosX', 'InitialPosY', 'InitialPosZ', 'InitialDirX', 'InitialDirY', 'InitialDirZ'};

reset(ttds)

data = read(ttds);

counter = 1;
while hasdata(ttds)
    if counter > 100
        break
    end
    % Read in Chunk
    dataChunk = read(ttds);
    data = [data; dataChunk];
    counter = counter + 1;
    disp(num2str(counter));
end
 

%% rearrange data

P0 = data{:,1:3}'; % initial positions of photons
u0 = data{:,4:6}'; % directions of photons

%% adaption of plane_line_intersect function for matrix multiplication

z_interface = 5; % water/air interface z-distance from (0,0,0) in mm

n = [0;0;-1]; % normal vector of water/air interface
n = repmat(n, 1, size(P0, 2));

V0 = [0; 0; z_interface]; % point in plane
V0 = repmat(V0, 1, size(P0, 2));

w0 = P0 - V0; % vector from photon origin to plane
D0 = dot(n,u0);
N0 = -dot(n,w0);

sI0 = N0./D0;
P1 = P0 + sI0.*u0;

%% snells law

n_i = 1.33; % watter index of refraction
n_f = 1; % air index of refraction

mu = n_i/n_f;

% refract rays (only works on unit vectors)
u1 = (sqrt(1 - ((mu^2) * ( 1 - (dot(n,u0).^2) ) )) .* n) + mu.*(u0 - dot(n,u0) .* n);
u1(imag(u1)~=0) = nan; % incident angles > theta_crit come out as complex, remove
nan_ind = isnan(u1);
P1(nan_ind) = nan;
clear nan_ind
% Now, P1 are positions of the photons at the water/air interface, with
% vectors u1

%% Now, plane - line intersection with refracted rays

z_screen = 10; % water/air interface z-distance from (0,0,0) in mm

n = [0;0;-1]; % normal vector of screen plane (x-y)
n = repmat(n, 1, size(P1, 2));

V1 = [0; 0; z_screen]; % point in plane
V1 = repmat(V1, 1, size(P1, 2));

w1 = P1 - V1; % vector from photon origin to plane
D1 = dot(n,u1);
N1 = -dot(n,w1);

sI1 = N1./D1;
P2 = P1 + sI1.*u1;

%% Data

scatter_data_ref = P2(1:2,:)';
s = num2str(z_screen/10);

%% Create image

edges = -300:1:300;

intMap_ref = hist3(scatter_data_ref, 'Edges',{edges, edges});

figure;
imagesc(intMap_ref);
axis image
colormap(jet); colorbar;
title([s, ' cm - with refraction'])
xlabel('x (mm)');
ylabel('y (mm)');
set(gca, 'Fontsize', 16)
saveas(gcf,['figs/',s,'cm_refraction.tiff'])


%% Without refraction

n = [0;0;-1]; % normal vector of water/air interface
n = repmat(n, 1, size(P0, 2));

V1 = [0; 0; z_screen]; % point in plane
V1 = repmat(V1, 1, size(P0, 2));

w_alt = P0 - V1; % vector from photon origin to plane
D0 = dot(n,u0);
N_alt = -dot(n,w_alt);

sI_alt = N_alt./D0;
P_alt = P0 + sI_alt.*u0;


%% Data - alt

scatter_data_alt = P_alt(1:2,:)';

%% Create image - alt

edges = -300:1:300;

intMap_ref = hist3(scatter_data_alt, 'Edges',{edges, edges});

figure;
imagesc(intMap_ref);
axis image
colormap(jet); colorbar;
title([s, ' cm - without refraction'])
xlabel('x (mm)');
ylabel('y (mm)');
set(gca, 'Fontsize', 16)
saveas(gcf,['figs/',s,'cm_norefraction.tiff'])






































