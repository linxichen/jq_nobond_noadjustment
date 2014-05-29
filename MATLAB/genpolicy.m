clear all
clc 
close all


%% Parameters that are used
% Bgrid = csvread('Bgridcuda.csv');
Kgrid = csvread('../Kgridcuda.csv');
Zgrid = csvread('../Zgridcuda.csv');
XXIgrid = csvread('../XXIgridcuda.csv');
P = csvread('../Pcuda.csv');

nk = length(Kgrid); % nb = length(Bgrid);
nz = length(Zgrid); nxi = length(XXIgrid);
parameters;
%% Read in shadow value and policy functions
MK_low = csvread('../V1_low.csv');
% MC_low = csvread('V3_low.csv');
% 
MK_high = csvread('../V1_high.csv');
% MC_high = csvread('V3_high.csv');

koptcuda = csvread('../koptcuda.csv');
coptcuda = csvread('../coptcuda.csv');
Rcuda = csvread('../Rcuda.csv');
dcuda = csvread('../dcuda.csv');
ncuda = csvread('../ncuda.csv');
mmucuda = csvread('../mmucuda.csv');
wagecuda = csvread('../wagecuda.csv');

flagcuda = csvread('../flagcuda.csv');
% EEerrorcuda = csvread('EEerror.csv');

%% Reshape
MK_low = reshape(MK_low,[nk,nz,nxi]);
MK_high = reshape(MK_high,[nk,nz,nxi]);

koptcuda = reshape(koptcuda,[nk,nz,nxi]);
coptcuda = reshape(coptcuda,[nk,nz,nxi]);
%Rcuda = reshape(Rcuda,[nk,nz,nxi]);
dcuda = reshape(dcuda,[nk,nz,nxi]);
ncuda = reshape(ncuda,[nk,nz,nxi]);
wagecuda = reshape(wagecuda,[nk,nz,nxi]);
mmucuda = reshape(mmucuda,[nk,nz,nxi]);

P = reshape(P,[nz*nxi nz*nxi]);

flagcuda = reshape(flagcuda,[nk,nz,nxi]);

%% Find the gap between high and low
dist_K = MK_high - MK_low;
figure
surf(Zgrid,Kgrid(end-20:end),squeeze(dist_K(end-20:end,:,ceil(nxi/2))),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('Dist of Shadow Val Investment')

dist_K = MK_high - MK_low;
figure
surf(Zgrid,Kgrid(end-20:end),squeeze(flagcuda(end-20:end,:,ceil(nxi/2))),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('Flags')


%%
% Linxi: Paint Kopt for i_z=5, i_xxi=5

% One slice of captial policy
figure
tightness = 1;
tfp = 7;
subplot(3,3,1)
plot(Kgrid,koptcuda(:,tfp,tightness)-Kgrid);
xlabel('Capital'); ylabel('Change in Capital')
subplot(3,3,2)
plot(Kgrid,coptcuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Consumption')
subplot(3,3,3)
plot(Kgrid,wagecuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Wage');
subplot(3,3,4)
plot(Kgrid,dcuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Dividend');
subplot(3,3,5)
plot(Kgrid,koptcuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Capital Tmr');
subplot(3,3,6)
plot(Kgrid,ncuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Hours');
subplot(3,3,7)
plot(Kgrid,mmucuda(:,tfp,tightness));
xlabel('Capital'); ylabel('LM');
subplot(3,3,8)
plot(Kgrid,koptcuda(:,tfp,tightness)-(1-ddelta)*Kgrid);
xlabel('Capital'); ylabel('Investment')
subplot(3,3,9)
plot(Kgrid,MK_low(:,tfp,tightness));
xlabel('Capital'); ylabel('Shadow Value')

leg=legend('Lower Shadow Value','High Shadow Value');
set(leg,'Orientation','horizontal','Position',[0.3144 0.0123 0.4014 0.02725],'FontSize',10);
legend boxoff;

Kgrid_occa = Kgrid;
Zgrid_occa = Zgrid;
XIgrid_occa = XIgrid;
koptcuda_occa = koptcuda;
coptcuda_occa = coptcuda;
wagecuda_occa = wagecuda;
dcuda_occa = dcuda;
ncuda_occa = ncuda;
mmucuda_occa = mmucuda;

figure
surf(Zgrid,Kgrid,squeeze(dist_K(:,:,ceil(nxi/2))),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('Dist of Shadow Val Investment')


%%
% Shrinkage
dist_K = MK_high-MK_low;
% dist_B = MB_high-MB_low;

% 3D distance check
figure
plot(Kgrid,squeeze(dist_K(:,nz,1)))
xlabel('Capital'); ylabel('Distance K');
print -f -depsc2 '2ddistance1.eps'
title('i_z=5, i_b=1')

% 3D distance check
figure
plot(Kgrid,squeeze(dist_K(:,1,nxi)))
xlabel('Capital'); ylabel('Distance K');
print -f -depsc2 '2ddistance2.eps'
title('i_z=5, i_b=1')

% 3D distance check
figure
plot(Kgrid,squeeze(dist_K(:,tfp,ceil(nxi/2))))
xlabel('Capital'); ylabel('Distance K');
print -f -depsc2 '2ddistance2.eps'
title('i_z=5, i_b=1')

figure
surf(Zgrid,Kgrid,squeeze(dist_K(:,:,ceil(nxi/2))),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('Dist of Shadow Val Investment')
print -f -depsc2 '3ddistance_new.eps'

figure
surf(Zgrid,Kgrid,squeeze(dist_K(:,:,end)),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('Dist of Shadow Val Investment')
print -f -depsc2 '3ddistance_new.eps'

figure
surf(Zgrid,Kgrid,squeeze(flagcuda(:,:,ceil(nxi/2))))
xlabel('TFP Shock'); ylabel('Capital'); zlabel('flag')
print -f -depsc2 '3dflag_new.eps'

figure
surf(Zgrid,Kgrid,squeeze(mmucuda(:,:,ceil(nxi/2))),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('TFP Shock'); ylabel('Capital'); zlabel('\mu_t')
print -f -depsc2 '3dmmu_new.eps'
