clear all
clc 
close all

%% Parameters that are used
% Bgrid = csvread('Bgridcuda.csv');
Kgrid = csvread('Kgrid.csv');
Zgrid = csvread('Zgrid.csv');
XXIgrid = csvread('XXIgrid.csv');
coeff = csvread('coeff.csv');
P = csvread('P.csv');

nk = length(Kgrid); % nb = length(Bgrid);
nz = length(Zgrid); nxxi = length(XXIgrid);
cd ../MATLAB
fpiter_para;
cd ../fpiter_results;
accuracy;

%% Read in shadow value and policy functions
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
koptcuda = reshape(koptcuda,[nk,nz,nxxi]);
coptcuda = reshape(coptcuda,[nk,nz,nxxi]);

%Rcuda = reshape(Rcuda,[nk,nz,nxi]);
%dcuda = reshape(dcuda,[nk,nz,nxi]);
ncuda = reshape(ncuda,[nk,nz,nxxi]);
%wagecuda = reshape(wagecuda,[nk,nz,nxi]);
mmucuda = reshape(mmucuda,[nk,nz,nxxi]);
coeff = reshape(coeff,[pk+1 pz+1 pxxi+1]);

P = reshape(P,[nz*nxxi nz*nxxi]);

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
%plot(Kgrid,moptcuda(:,tfp,tightness));
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
u = rand(1,T);
kindex = ones(1,T);
zindex = ones(1,T);
xxiindex = ones(1,T);
k = kss*ones(1,T);

for t = 1:T-1
    k(t+1) = koptcuda(kindex(t),zindex(t),xxiindex(t));
    cdf = cumsum(P(sub2ind([nz nxxi],zindex(t),xxiindex(t)),:));
    agg_shock = find(cdf>u(t),1,'first');
    [zindex(t+1),xxiindex(t+1)] = ind2sub([nz nxxi],agg_shock);
    [~,kindex(t+1)] = min(abs(Kgrid-k(t+1))); 
end

eee = zeros(1,T-burnin);
for i = burnin+1:T
    c = coptcuda(kindex(i),zindex(i),xxiindex(i));
    kplus = koptcuda(kindex(i),zindex(i),xxiindex(i));
    [~,i_kplus] = min(abs(Kgrid-kplus)); 
    sum = 0;
    for i_zplus = 1:nz
        for i_xxiplus = 1:nxxi
            sum = sum + bbeta*P(sub2ind([nz nxxi],zindex(i),xxiindex(i)),sub2ind([nz nxxi],i_zplus,i_xxiplus))*(1-ddelta+(1-mmucuda(i_kplus,i_zplus,i_xxiplus))*ttheta*Zgrid(i_zplus)*Kgrid(i_kplus)^(ttheta-1)*ncuda(i_kplus,i_zplus,i_xxiplus)^(1-ttheta))/(coptcuda(i_kplus,i_zplus,i_xxiplus));
        end
    end
    ctilde = (1-mmucuda(kindex(i),zindex(i),xxiindex(i))*XXIgrid(xxiindex(i)))/sum;
    eee(i) = abs(c/ctilde-1);
end

eulererror = mean(eee(burnin+1:end))


