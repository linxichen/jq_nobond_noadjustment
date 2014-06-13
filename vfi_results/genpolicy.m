clear all
clc 
close all


%% Parameters that are used
% Bgrid = csvread('Bgridcuda.csv');
Kgrid = csvread('Kgrid.csv');
Zgrid = csvread('Zgrid.csv');
XXIgrid = csvread('XXIgrid.csv');
P = csvread('P.csv');

nk = length(Kgrid); % nb = length(Bgrid);
nz = length(Zgrid); nxi = length(XXIgrid);
cd ../MATLAB
vfi_para;
cd ../vfi_results;
%% Read in shadow value and policy functions
V = csvread('V.csv');

koptcuda = csvread('kopt.csv');
coptcuda = csvread('copt.csv');
% Rcuda = csvread('../R.csv');
%dcuda = csvread('dopt.csv');
ncuda = csvread('nopt.csv');
mmucuda = csvread('mmuopt.csv');
%wagecuda = csvread('wopt.csv');

% flagcuda = csvread('flag.csv');
% EEerrorcuda = csvread('EEerror.csv');

%% Reshape
V = reshape(V,[nk,nz,nxi]);

koptcuda = reshape(koptcuda,[nk,nz,nxi]);
coptcuda = reshape(coptcuda,[nk,nz,nxi]);
%Rcuda = reshape(Rcuda,[nk,nz,nxi]);
%dcuda = reshape(dcuda,[nk,nz,nxi]);
ncuda = reshape(ncuda,[nk,nz,nxi]);
%wagecuda = reshape(wagecuda,[nk,nz,nxi]);
mmucuda = reshape(mmucuda,[nk,nz,nxi]);

P = reshape(P,[nz*nxi nz*nxi]);

% flagcuda = reshape(flagcuda,[nk,nz,nxi]);

%%
% Linxi: Paint Kopt for i_z=5, i_xxi=5

% One slice of captial policy
figure
tightness = 1;
tfp = 5;
subplot(3,3,1)
plot(Kgrid,koptcuda(:,tfp,tightness)-Kgrid);
xlabel('Capital'); ylabel('Change in Capital')
subplot(3,3,2)
plot(Kgrid,coptcuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Consumption')
subplot(3,3,3)
%plot(Kgrid,wagecuda(:,tfp,tightness));
xlabel('Capital'); ylabel('Wage');
subplot(3,3,4)
%plot(Kgrid,dcuda(:,tfp,tightness));
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
plot(Kgrid,V(:,tfp,tightness));
xlabel('Capital'); ylabel('Shadow Value')

leg=legend('Lower Shadow Value','High Shadow Value');
set(leg,'Orientation','horizontal','Position',[0.3144 0.0123 0.4014 0.02725],'FontSize',10);
legend boxoff;

Kgrid_occa = Kgrid;
Zgrid_occa = Zgrid;
XXIgrid_occa = XXIgrid;
koptcuda_occa = koptcuda;
coptcuda_occa = coptcuda;
%wagecuda_occa = wagecuda;
%dcuda_occa = dcuda;
ncuda_occa = ncuda;
mmucuda_occa = mmucuda;

%% Find Euler equation error
burnin = 1000;
T = 10000+burnin;
uz = rand(1,T);
uxxi = rand(1,T);



