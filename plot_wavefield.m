dt=0.001

fn = 'OUTPUT/wavefield_xvec1_for_dumps.bin';
fid = fopen(fn);
xvec1 = fread(fid, inf, 'float');
fclose(fid);
nx = length(xvec1)

fn = 'OUTPUT/wavefield_zvec1_for_dumps.bin';
fid = fopen(fn);
zvec1 = fread(fid, inf, 'float');
fclose(fid);
nz = length(zvec1)

istep=input('Please input istep:');

fn = sprintf('OUTPUT/wavefield%07d_Vz1.bin',istep)
title_1 = (['Vz1: ', num2str(istep*dt),'s']) 
fid=fopen(fn);
Vx1 = fread(fid, [nz nx], 'float');
fclose(fid);

figure;
subplot(1,2,1)
imagesc(xvec1, zvec1, Vx1);
axis ij; axis image;
colorbar;
title(title_1);


fn = sprintf('OUTPUT/wavefield%07d_Vz2.bin',istep)

title_2 = (['Vz2: ', num2str(istep*dt),'s']) 
fid=fopen(fn);
Vz2 = fread(fid, [nz nx], 'float');
fclose(fid);

subplot(1,2,2);
imagesc(xvec1, zvec1, Vz2);
axis ij; axis image;
title(title_2);
colorbar;

