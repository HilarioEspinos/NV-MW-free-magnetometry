% Estimate magnetic-field related parameters by brute-force Bayesian
% evaluation of a log-likelihood provided by `log_likelihood_derivatives`.
% Produces 1D marginals and selected 2D marginal maps.

% - Angles are in radians internally. Degrees are used only for axis labels.
% - Magnetic fields are treated in Tesla internally. Axis labels show mT.
% - A calibration factor (1.25) is applied to Bz to correct a known
% experimental miscalibration.
%
% Author: Hilario Espinós Martínez  

clear;clc

%% Load experimental data

% The PL map has already been normalized in a separate file; only the PL
% values are stored in the _NORMALIZED file
filename = 'random_axis_2/1mT_with_4AmpYaxis_2.48AmpXaxis_transverse_rotation_scan_merged_with_time_NORMALIZED.csv';
PL_exp = cell2mat(readcell(filename));   

% Raw file used to extract Bz axis (first column)
filename_raw = "random_axis_2/1mT_with_4AmpYaxis_2.48AmpXaxis_transverse_rotation_scan_merged_with_time.csv";
df = cell2mat(table2cell(readtable(filename_raw)));

Bz = 1.25*df(2:end,1)*1e-4; % Bz in Gauss => convert to Tesla with 1e-4, then apply calibration factor 1.25.
phi = linspace(0,355,size(PL_exp,1))*pi/180; % radians

% Subsampling & noise
step = 50; % Reduce data size if needed to speed up the brute-force sweep.
Bz = Bz(1:step:end);
PL_exp = PL_exp(:,1:step:end);
sigma = 0.0005; % Experimental noise scale for the Gaussian likelihood


%% Brute force bayesian estimation

% Define the grid to sweep for the posterior.
bz = 0 * 1e-3; % Tesla
Bperp = 1 *1e-3; % Tesla
alpha = linspace(0 , 360, 10) * pi/180; % Radians
beta =  acos(1/sqrt(3))-acos(linspace(1/sqrt(3),1, 50)); % Radians
zeta = linspace(0,120,50) * pi/180; % Radians

na = length(alpha);
nb = length(beta);
nz = length(zeta);
nbz = length(bz);
nbp = length(Bperp);

log_posterior = zeros(nbz,nbp,na,nb,nz);

tic
parfor ia = 1:na
    ia
    for ib = 1:nb
        for iz = 1:nz
            for ibz = 1:nbz
                for ibp = 1:nbp
                    log_posterior(ibz,ibp,ia,ib,iz) = log_likelihood_derivatives(Bz,phi,bz(ibz), Bperp(ibp), alpha(ia), beta(ib), zeta(iz), PL_exp,sigma);
                end
            end
        end
    end
end
toc

% Replace non-finite with a very negative number
log_posterior(~isfinite(log_posterior)) = -1e300;

% Substract the maximum to avoid overflows
log_posterior = log_posterior - max(log_posterior(:));
posterior_exp = exp(log_posterior);

% Marginalization
posterior_bz = squeeze(sum(sum(sum(sum(posterior_exp, 2), 3), 4), 5)); % sum over bperp, alpha, beta and zeta
posterior_Bperp = squeeze(sum(sum(sum(sum(posterior_exp, 1), 3), 4), 5));  % sum over bz, alpha, beta and zeta
posterior_alpha = squeeze(sum(sum(sum(sum(posterior_exp, 1), 2), 4), 5));  % sum over bz, bperp, beta and zeta
posterior_beta = squeeze(sum(sum(sum(sum(posterior_exp, 1), 2), 3), 5));  % sum over bz, bperp, alpha and zeta
posterior_zeta = squeeze(sum(sum(sum(sum(posterior_exp, 1), 2), 3), 4));    % sum over bz, bperp, alpha and beta


posterior_alphabeta = squeeze(sum(sum(sum(posterior_exp,1), 2), 5)); % Marginal over alpha and beta
posterior_alphazeta = squeeze(sum(sum(sum(posterior_exp,1), 2), 4)); % Marginal over alpha and zeta

% Marginals can be normalized to unit area for comparability; not needed
% for the purpose of this work.


%% Representation
figure
plot(bz*1e3, posterior_bz,'LineWidth',1.5)
xlabel('b_z (mT)')

figure
plot(Bperp*1e3, posterior_Bperp,'LineWidth',1.5)
xlabel('b_\perp (mT)')

figure
plot(alpha, posterior_alpha,'LineWidth',1.5)
xlabel('\alpha (rad)')

figure
plot(beta, posterior_beta,'LineWidth',1.5)
xlabel('\beta (rad)')

figure
plot(zeta, posterior_zeta,'LineWidth',1.5)
xlabel('\zeta (rad)')


figure
pcolor(alpha*180/pi, beta*180/pi, posterior_alphabeta')
shading flat
xlabel('\alpha (deg)'); ylabel('\beta (deg)')


figure
pcolor(alpha*180/pi, zeta*180/pi, posterior_alphazeta')
shading flat
xlabel('\alpha (deg)'); ylabel('\zeta (deg)')