clear all;
close all;

cmap=load('jet_modn');         % colorscale


fid = fopen('configure','r');
dt = fscanf(fid,'%f',1);fgets(fid);
nt = fscanf(fid,'%d',1);fgets(fid);
xmin = fscanf(fid,'%f',1); dx = fscanf(fid,'%f',1); nx = fscanf(fid,'%f',1); fgets(fid);
zmin = fscanf(fid,'%f',1); dz = fscanf(fid,'%f',1); nz = fscanf(fid,'%f',1); fgets(fid);
fclose(fid);

x_range = xmin+dx*(nx-1);
z_range = zmin+dz*(nz-1);

fn = '../OUTPUT/wavefield_xvec1_for_dumps.bin';
fid = fopen(fn);
xvec1 = fread(fid, inf, 'float');
fclose(fid);
nx = length(xvec1)

fn = '../OUTPUT/wavefield_zvec1_for_dumps.bin';
fid = fopen(fn);
zvec1 = fread(fid, inf, 'float');
fclose(fid);
nz = length(zvec1)

fn = '../OUTPUT/wavefield_xvec2_for_dumps.bin';
fid = fopen(fn);
xvec2 = fread(fid, inf, 'float');
fclose(fid);

fn = '../OUTPUT/wavefield_zvec2_for_dumps.bin';
fid = fopen(fn);
zvec2 = fread(fid, inf, 'float');
fclose(fid);

istep=input('Please input istep:');

fn = sprintf('../OUTPUT/wavefield%07d_Vz1.bin',istep);
fid=fopen(fn);
Vz1 = fread(fid, [nz nx], 'float');
fclose(fid);

fn = sprintf('../OUTPUT/wavefield%07d_Vz2.bin',istep);
fid=fopen(fn);
Vz2 = fread(fid, [nz nx], 'float');
fclose(fid);

fn = sprintf('../OUTPUT/wavefield%07d_Vx1.bin',istep);
fid=fopen(fn);
Vx1 = fread(fid, [nz nx], 'float');
fclose(fid);

fn = sprintf('../OUTPUT/wavefield%07d_Vx2.bin',istep);
fid=fopen(fn);
Vx2 = fread(fid, [nz nx], 'float');
fclose(fid);

xvec_all=xvec1(1):dx/2:xvec2(end);
zvec_all=zvec1(1):dz/2:zvec2(end);
nx_all = length(xvec_all);
nz_all = length(zvec_all);
n_plot = min(length(xvec_all),length(zvec_all));

% can not handle dx != dz
x_plot_45 = xvec_all(2:end);
z_plot_45 = x_plot_45-dx/2;
x_plot_135 = xvec_all(2:end);
z_plot_135 = zvec_all(end)+xvec_all(1) - x_plot_45;
for ix = 1:nx_all
    for iz = 1:nz_all
        if ~isempty(find(zvec1==zvec_all(iz))) && ~isempty(find(xvec2==xvec_all(ix)))
            Vx_all(iz,ix) = Vx1(find(zvec1==zvec_all(iz)),find(xvec2==xvec_all(ix)));
            Vz_all(iz,ix) = Vz2(find(zvec1==zvec_all(iz)),find(xvec2==xvec_all(ix)));
        end
        if ~isempty(find(zvec2==zvec_all(iz))) && ~isempty(find(xvec1==xvec_all(ix)))
            Vx_all(iz,ix) = Vx2(find(zvec2==zvec_all(iz)),find(xvec1==xvec_all(ix)));
            Vz_all(iz,ix) = Vz1(find(zvec2==zvec_all(iz)),find(xvec1==xvec_all(ix)));
        end
    end
end

for i = 1:n_plot-1
    Vx_plot_45(i) = Vx_all(find(zvec_all==z_plot_45(i)),find(xvec_all==x_plot_45(i)));
    Vz_plot_45(i) = Vz_all(find(zvec_all==z_plot_45(i)),find(xvec_all==x_plot_45(i)));
    Vx_plot_135(i) = Vx_all(find(zvec_all==z_plot_135(i)),find(xvec_all==x_plot_135(i)));
    Vz_plot_135(i) = Vz_all(find(zvec_all==z_plot_135(i)),find(xvec_all==x_plot_135(i)));
end

%============== plot
figure;

nSeries = 8;
% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
rowH = 0.73 / nRow ;  colW = 0.95 / nCol ;
colX = 0.02 + linspace( 0, 0.9 , nCol+1 ) ;  colX = colX(1:end-1) ;
rowY = 0.06 + linspace( 1, 0.05, nRow+1 ) ;  rowY = rowY(2:end) ;

%subplot(2,4,1)
dId = 1;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec2, zvec1, Vx1);
axis([xmin x_range zmin z_range])
axis ij; axis image;
colormap(cmap);
title(['Vx1: ', num2str(istep*dt),'s']);
lim=caxis

%subplot(2,4,2);
dId = 2;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec1, zvec2, Vx2);
axis([xmin x_range zmin z_range]);
axis ij; axis image;
title(['Vx2: ', num2str(istep*dt),'s']);
caxis(lim);
colormap(cmap);

%subplot(2,4,3);
dId = 3;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec_all, zvec_all, Vx_all);
axis ij; axis image;
axis([xmin x_range zmin z_range])
caxis(lim);
colormap(cmap);
hold on;
plot(x_plot_45,z_plot_45,'r--','linewidth',1.2);
hold on;
plot(x_plot_135,z_plot_135,'b--','linewidth',1.2);
title(['Vx\_all']);

%--- waveform
axes('position',[colX(colId)+colW+0.02,rowY(rowId)-0.01,0.25,0.17]);
plot(x_plot_45,Vx_plot_45,'r-','linewidth',1.2);
%hold on;
%plot([0 x_range], [0 0],'k--','linewidth',0.8);
xlim([xmin x_range]);
xlabel('Distance (m)','fontsize',12);
%xlim([X0-(it*dt-0.05)*Vp_max X0+(it*dt-0.05)*Vp_max]);
axes('position',[colX(colId)+colW+0.02,rowY(rowId)+0.2,0.25,0.17]);
%plot([0 x_range], [0 0],'k--','linewidth',0.8);
%hold on;
plot(x_plot_135,Vx_plot_135,'b-','linewidth',1.2);
xlim([xmin x_range]);
title('Vx')
%xlim([X0-(it*dt-0.05)*Vp_max X0+(it*dt-0.05)*Vp_max]);



%subplot(2,4,5)
dId = 5;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec1, zvec2, Vz1);
axis([xmin x_range zmin z_range])
axis ij; axis image;
colormap(cmap);
title(['Vz1: ', num2str(istep*dt),'s']);
lim=caxis

%subplot(2,4,6);
dId = 6;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec2, zvec1, Vz2);
axis([xmin x_range zmin z_range])
axis ij; axis image;
title(['Vz2: ', num2str(istep*dt),'s']);
caxis(lim);
colormap(cmap);

%subplot(2,4,7);
dId = 7;
rowId = ceil( dId / nCol ) ;
colId = dId - (rowId - 1) * nCol ;
axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
imagesc(xvec_all, zvec_all, Vz_all);
axis([xmin x_range zmin z_range])
axis ij; axis image;
caxis(lim);
colormap(cmap);
hold on;
plot(x_plot_45,z_plot_45,'r--','linewidth',1.2);
hold on;
plot(x_plot_135,z_plot_135,'b--','linewidth',1.2);
title(['Vz\_all']);
hbar = colorbar('horiz');
get(hbar,'position');
set(hbar,'position',[0.15 0.05 0.4 0.02],'fontsize',10);

%---- waveforms 45 deg and 135 deg
axes('position',[colX(colId)+colW+0.02,rowY(rowId)-0.02,0.25,0.17]);
plot(x_plot_45,Vz_plot_45,'r-','linewidth',1.2);
%hold on;
%plot([0 x_range], [0 0],'k--','linewidth',0.8);
xlim([xmin x_range]);
xlabel('Distance (m)','fontsize',12);
%xlim([X0-(it*dt-0.05)*Vp_max X0+(it*dt-0.05)*Vp_max]);
axes('position',[colX(colId)+colW+0.02,rowY(rowId)+0.19,0.25,0.17]);
%plot([xmin x_range], [0 0],'k--','linewidth',0.8);
%hold on;
plot(x_plot_135,Vz_plot_135,'b-','linewidth',1.2);
xlim([xmin x_range]);
title('Vz')
%xlim([X0-(it*dt-0.05)*Vp_max X0+(it*dt-0.05)*Vp_max]);


set(gcf,'position',[1024 1024 2000 1024]);
set( gcf, 'Color', 'White') ;


fn=['lebedev_it',num2str(istep),'.eps'];
print(gcf,'-depsc',fn);
