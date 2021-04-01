close all
clear all
clc
%% Philip Mocz (2020), Harvard University
% check output of arepo runs

% parameters
N = 128;          % resolution
%snap = 20;
snaps = [0 95];

outpath = '/mnt/tigress/MHD/subcriticalMHD_2020/output/mhd128GB4M10/';


%%

Mach = 10;
boxSize =  1;
dx = boxSize/N;

clim = [-2 2];

Nsim = numel(snaps);


%% read info

rhoProj = cell(Nsim,1);


for i = 1:Nsim
    i
    
    snap = snaps(i);
    
    filename = [outpath  'snap_' sprintf('%03d',snap) '.hdf5'];
    
    pos = h5read(filename,'/PartType0/Coordinates');
    x = pos(1,:)';
    y = pos(2,:)'; y(y==0) = dx/100;
    z = pos(3,:)';
    clear pos;
    m = h5read(filename,'/PartType0/Masses');
    rhoProj{i} = accumarray([ceil(y/dx) ceil(x/dx)], m/dx^3, [N N]) / N;
    clear x;
    clear y;
    clear z;
    clear m;
    rho{i} = h5read(filename,'/PartType0/Density');
    mass{i} = h5read(filename,'/PartType0/Masses');
    
end

%% Plot Projected Densities
for i = 1:Nsim
    
    snap = snaps(i);
    
    figure;
    imagesc(log(rhoProj{i})')
    caxis(clim)
    title(['snap=' num2str(snap)],'interpreter','latex')
    set(gca,'YDir','normal')
    axis square
    axis off
    
    cmap = parula(256);
    savname = ['../writeup/snap' num2str(num2str(snap)) '.png'];
    my_imwrite(log(rhoProj{i}),cmap,clim,savname)
end


%%  
fh = figure;

bins = -8:0.1:8;

for i = 1:Nsim
    
    snap = snaps(i);
    
    H = myBin1D(log(rho{i}),mass{i},bins);
    
    ph{i} = semilogy(bins, H);
    hold on
    
end

axis([-8 8 1e-8 1e-1])
lh = legend(['snap=' num2str(snaps(1))],['snap=' num2str(snaps(2))]);
set(lh,'location','northwest')

xlabel('$\ln(\rho)$','interpreter','latex')
ylabel('pdf','interpreter','latex')

saveas(fh,'../writeup/densitypdf.eps','epsc2')