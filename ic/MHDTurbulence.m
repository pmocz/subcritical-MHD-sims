clear all;
close all;
clc;

N_linear = 256;128; % 64
N_total  = N_linear^3;
Bstrength = 4;  1;2;4;0.5;0.25; % P_B = B^2/2
icname   = ['ics_MHDTurbulence_B' num2str(Bstrength) '_' num2str(N_linear) '.dat.hdf5'];

boxsize = 1.0;


x = zeros(N_total,1);
y = zeros(N_total,1);
z = zeros(N_total,1);
vx = zeros(N_total,1);
vy = zeros(N_total,1);
vz = zeros(N_total,1);
Bx = zeros(N_total,1);
By = zeros(N_total,1);
Bz = zeros(N_total,1);
utherm = zeros(N_total,1);
m = ones(N_total,1)/N_total;
ids = (1:N_total)';


%% Setting Particles
index = 1;
for i=1:N_linear
    i
    for j=1:N_linear
        for k=1:N_linear
            x(index) = (1.0*(i-1)+0.5)/(1.0*N_linear) * boxsize;
            y(index) = (1.0*(j-1)+0.5 + 0.48*mod((i-1),2))/(1.0*N_linear) * boxsize;
            z(index) = (1.0*(k-1)+0.5 + 0.46*mod((i-1),2))/(1.0*N_linear) * boxsize;
            vx(index) = 0.0;
            vy(index) = 0.0;
            vz(index) = 0.0;
            Bx(index) = 1.0*Bstrength;
            By(index) = 0.0;
            Bz(index) = 0.0;
            utherm(index) = 1000.0;
            index = index + 1;
        end
    end
end


%% Write snapshot


filename = icname;
delete(filename);

N = N_total;
Lx = boxsize;
pos = [x'; y'; z'];
vel = [vx'; vy'; vz'];
B = [Bx'; By'; Bz'];

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

h5disp(filename);

