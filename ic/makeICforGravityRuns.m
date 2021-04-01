clear all;
close all;
clc;

startsnap = 2;
N_linear = 256;  64;
Bstrength =  100; 1; 10; 30; 100; % P_B = B^2/ 8 / pi
Mach = 10;
isCT = 0; 1;
switch Bstrength
    case 0.1
        Bname = 'p1';
    otherwise
        Bname = num2str(Bstrength);
end

switch isCT
    case 0
        filein   = ['/mnt/hernquistfs3/Output/Arepo/ISMTurbulence/mhd' num2str(N_linear) 'B' Bname 'M' num2str(Mach) '/snap_' sprintf('%.3d',startsnap) '.hdf5'];
    case 1
        filein   = ['/mnt/hernquistfs3/Output/Arepo/ISMTurbulence/mhdCT' num2str(N_linear) 'B' Bname 'M' num2str(Mach) '/snap_' sprintf('%.3d',startsnap) '.hdf5'];
end
icname   = ['ics_MHDTurbulence_G_B' Bname 'M' num2str(Mach) '_' num2str(N_linear) '_' sprintf('%.3d',startsnap) '.dat.hdf5'];

boxsize = 1.0;
N_total  = N_linear^3;

%% Load particles from simulation
datasetname = '/PartType0/Coordinates';
pos = double(h5read(filein,datasetname));
x = pos(1,:)';
y = pos(2,:)';
z = pos(3,:)';
datasetname = '/PartType0/Velocities';
vel = double(h5read(filein,datasetname));
vx = vel(1,:)';
vy = vel(2,:)';
vz = vel(3,:)';
datasetname = '/PartType0/MagneticField';
Bfld = double(h5read(filein,datasetname));
Bx = Bfld(1,:)';
By = Bfld(2,:)';
Bz = Bfld(3,:)';
datasetname = '/PartType0/InternalEnergy';
utherm = double(h5read(filein,datasetname));
datasetname = '/PartType0/Masses';
m = double(h5read(filein,datasetname));
datasetname = '/PartType0/ParticleIDs';
ids = int64(h5read(filein,datasetname));

if isCT
    datasetname = '/PartType0/MagneticVectorPotential';
    A = double(h5read(filein,datasetname));
    Ax = A(1,:)';
    Ay = A(2,:)';
    Az = A(3,:)';
    datasetname = '/PartType0/MagneticAShiftX';
    ASx = double(h5read(filein,datasetname));
    ASxx = ASx(1,:)';
    ASxy = ASx(2,:)';
    ASxz = ASx(3,:)';
    datasetname = '/PartType0/MagneticAShiftY';
    ASy = double(h5read(filein,datasetname));
    ASyx = ASy(1,:)';
    ASyy = ASy(2,:)';
    ASyz = ASy(3,:)';
    datasetname = '/PartType0/MagneticAShiftZ';
    ASz = double(h5read(filein,datasetname));
    ASzx = ASz(1,:)';
    ASzy = ASz(2,:)';
    ASzz = ASz(3,:)';
end


%% Write snapshot


filename = icname;
delete(filename);

N = N_total;
Lx = boxsize;
pos = [x'; y'; z'];
vel = [vx'; vy'; vz'];
B = [Bx'; By'; Bz'];
if isCT
    A = [Ax'; Ay'; Az'];
    ASx = [ASxx'; ASxy'; ASxz'];
    ASy = [ASyx'; ASyy'; ASyz'];
    ASz = [ASzx'; ASzy'; ASzz'];
end

NumPart_ThisFile       = [N 0 0 0 0 0];
NumPart_Total          = [N 0 0 0 0 0];
NumPart_Total_HighWord = [0 0 0 0 0 0];
MassTable              = [0 0 0 0 0 0];
Time                   = 0;
Redshift               = 0;
BoxSize                = Lx;
NumFilesPerSnapshot    = 1;
Omega0                 = 0;
OmegaLambda            = 0;
HubbleParam            = 0;
Flag_Sfr               = 0;
Flag_Cooling           = 0;
Flag_StellarAge        = 0;
Flag_Metals            = 0;
Flag_Feedback          = 0;


h5create(filename,'/Header/data',[1],'datatype','single');
h5create(filename,'/PartType0/Coordinates',[3 N],'datatype','single');
h5create(filename,'/PartType0/ParticleIDs',[N],'datatype','uint32');
h5create(filename,'/PartType0/Masses',[N],'datatype','single');
h5create(filename,'/PartType0/Velocities',[3 N],'datatype','single');
h5create(filename,'/PartType0/MagneticField',[3 N],'datatype','single');
h5create(filename,'/PartType0/InternalEnergy',[N],'datatype','single');
if isCT
    h5create(filename,'/PartType0/MagneticVectorPotential',[3 N],'datatype','single');
    h5create(filename,'/PartType0/MagneticAShiftX',[3 N],'datatype','single');
    h5create(filename,'/PartType0/MagneticAShiftY',[3 N],'datatype','single');
    h5create(filename,'/PartType0/MagneticAShiftZ',[3 N],'datatype','single');
end

h5writeatt(filename,'/Header','NumPart_ThisFile',uint32(NumPart_ThisFile));
h5writeatt(filename,'/Header','NumPart_Total',uint32(NumPart_Total));
h5writeatt(filename,'/Header','NumPart_Total_HighWord',uint32(NumPart_Total_HighWord));
h5writeatt(filename,'/Header','MassTable',single(MassTable));
h5writeatt(filename,'/Header','Time',single(Time));
h5writeatt(filename,'/Header','Redshift',single(Redshift));
h5writeatt(filename,'/Header','BoxSize',single(BoxSize));
h5writeatt(filename,'/Header','NumFilesPerSnapshot',uint32(NumFilesPerSnapshot));
h5writeatt(filename,'/Header','Omega0',single(Omega0));
h5writeatt(filename,'/Header','OmegaLambda',single(OmegaLambda));
h5writeatt(filename,'/Header','HubbleParam',single(HubbleParam));
h5writeatt(filename,'/Header','Flag_Sfr',uint32(Flag_Sfr));
h5writeatt(filename,'/Header','Flag_Cooling',uint32(Flag_Cooling));
h5writeatt(filename,'/Header','Flag_StellarAge',uint32(Flag_StellarAge));
h5writeatt(filename,'/Header','Flag_Metals',uint32(Flag_Metals));
h5writeatt(filename,'/Header','Flag_Feedback',uint32(Flag_Feedback));
h5writeatt(filename,'/Header','Flag_DoublePrecision',uint32(0));

h5write(filename,'/PartType0/Coordinates',single(pos));
h5write(filename,'/PartType0/ParticleIDs',uint32(ids));
h5write(filename,'/PartType0/Masses',single(m));
h5write(filename,'/PartType0/Velocities',single(vel));
h5write(filename,'/PartType0/MagneticField',single(B));
h5write(filename,'/PartType0/InternalEnergy',single(utherm));

if isCT
    h5write(filename,'/PartType0/MagneticVectorPotential',single(A));
    h5write(filename,'/PartType0/MagneticAShiftX',single(ASx));
    h5write(filename,'/PartType0/MagneticAShiftY',single(ASy));
    h5write(filename,'/PartType0/MagneticAShiftZ',single(ASz));
end


h5disp(filename);

