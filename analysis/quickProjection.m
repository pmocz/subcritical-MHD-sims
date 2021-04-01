close all
clear all
clc
%% Philip Mocz (2020), Harvard University
% check output of arepo runs

% parameters
N = 128;          % resolution
Bs = [0.25 0.5 1 2 4] ;
Mach = 10;
snap = 20;

outpath = '/mnt/tigress/MHD/subcriticalMHD_2020/output/';


%%

Nsim = numel(Bs);

boxSize =  1;
dx = boxSize/N;



%% read info

rhoProj = cell(Nsim,1);

for i = 1:Nsim
    
    B = Bs(i)
    
    filename = [outpath 'mhd' num2str(N) 'B' num2str(B) 'M' num2str(Mach) '/' 'snap_' sprintf('%03d',snap) '.hdf5'];
    
    pos = h5read(filename,'/PartType0/Coordinates');
    x = pos(1,:)';
    y = pos(2,:)';
    z = pos(3,:)';
    clear pos;
    m = h5read(filename,'/PartType0/Masses');
    rhoProj{i} = accumarray([round(y/dx-.5)+1 round(x/dx-.5)+1], m/dx^3, [N N]) / N;
    clear x;
    clear y;
    clear z;
    clear m;
    
end

%%
figure;
for i = 1:Nsim
    subplot(1,Nsim,i)
    imagesc(log(rhoProj{i})')
    caxis([-2 2])
    title(['$B_0=$' num2str(Bs(i))],'interpreter','latex')
    set(gca,'YDir','normal')
    axis square
    axis off
end

