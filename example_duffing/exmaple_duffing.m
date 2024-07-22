% An example of BEEPSIS workflow using the experimental data 
% of a nanoprticle trapped in the optical tweezers in vacuum. 
% The optical force in this case correspond to the Duffing oscillator with 
% softening nonlinearity. For the deatils of the experiemntal setup see
% Flajsmanova Scientific Reports 2020. 
% The procedure described here is a basis for data analysis described in
% Section IV.A and figure 3 of BEEPSIS paper Siler Phys Rev Applied 2023, 
% https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.19.064059
% 
% The supplied data file '2MHz_20180806_150419.txt' contains one meaurement
% performed at pressure 10 Pa (0.1 mbar) and is ~1.7% of all data used for
% analysis presented in the published paper (at given pressure)
% 
% This file is part of the BEEPSIS toolbox.
% See LICENSE.md for information about using/distributing this file.

%%
clear

%% definitions
rho        = 2000; % particle density kg/m^3
a0         = 170e-9/2; % particle radius [m]
MM         = rho*4/3*pi*a0^3; % particle mass
T0         = 293; % environment temperature
kT_M       = 1.38e-23*T0/MM; % k_B T over mass
pressure   = 10; %medium pressure [Pa]
GammaTh    = 437.6804; % theoretical value of damping rate \Gamma [s] (see Li NaturePhysics 2011 and https://www.engineeringtoolbox.com/nitrogen-N2-dynamic-kinematic-viscosity-temperature-pressure-d_2067.html?vA=25&degree=C)

fsample    = 1764705.882353; % sampling frequency
dt         = 1/fsample; % time step
freq_block = 500; % number of frequencies for PSD blocking (see Berg-Sørensen, RevSciInstruments 2004)

calibration_factor = 1.9019; % calibration factor from PSD voltage to micrometers (see Flajsmanova SciRep 2020)
bandpass           = [80 100]*1e3; % frequencies for PSD bandpass filtering

%% load experimental data (voltages from quadrant photodiode)
datafile = '2MHz_20180806_150419.txt';
rawdata  = load(datafile);
data     = rawdata(:,1) * calibration_factor;


%% power spectraldensity and data filtering
T        =   (length(data)-1)*dt;
f        =   (0 : length(data)-1)' / T;
data     = data - mean(data);
FT       =   dt*fft(data);
P        =   FT .* conj(FT) / T;
ind      =   f <= fsample/2; % only to the Nyquist f 
f        =   f(ind);
P        =   P(ind);
% frequency blocking (see Berg-Sørensen RevSciInstruments 2004, Tolić-Nørrelykke Comp. Phys. Commun. 2004)
nbin     =   floor(length(f)/freq_block);
[fb, Pb] = deal(zeros(nbin, 1));
for kk = 1 : nbin
    fb(kk)   = mean(f((kk-1)*freq_block+1 : kk*freq_block));
    Pb(kk)   = mean(P((kk-1)*freq_block+1 : kk*freq_block));
end

% bandpass filter in frequency domain
FTfilt      = zeros(size(FT));
ind         = f>=bandpass(1) & f<=bandpass(2);
FTfilt(ind) = FT(ind);
data_filt   = ifft(FTfilt,'symmetric')/dt; 

%% PSD figure
figure
semilogy(fb/1e3, Pb, 'LineWidth',2);
xlim([50 150])
ax = gca;
hold on
plot(bandpass(1) * [1 1]/1e3, ax.YLim, 'k--', 'LineWidth',2)
plot(bandpass(2) * [1 1]/1e3, ax.YLim, 'k--', 'LineWidth',2)
hold off
grid on
xlabel('\omega/(2\pi) [Hz]', 'Interpreter','tex')
ylabel('{\itP}_{xx} [\mum^2Hz^{-1}]', 'Interpreter','tex')
clear ax
%% apply non-parametric BEEPSIS and check force profile 

% non-parametric BEEPSIS
bin_res = beepsis_bin(data_filt*1e-6,dt, 31, GammaTh, 'c');
% find average linear force
stiff = nlinfit(bin_res.bincenter,bin_res.force, @(p,x) p*x, -1);

%% plot obtained force and check deviations from linear profile
figure
tiledlayout(2,1, "TileSpacing","tight",'Padding','tight')
nexttile
% non-parametic force andits linear fit
plot(bin_res.bincenter*1e6,bin_res.force, 'LineWidth', 2);
hold on
plot(bin_res.bincenter*1e6,stiff * bin_res.bincenter, 'LineWidth', 2);
hold off
grid on
ylabel('{\itF_x/m} [ms^{-2}]','Interpreter','tex')
legend('non-parm. BEEPSIS', 'linear fit')
nexttile
% difference showing the cubic 
plot(bin_res.bincenter*1e6,bin_res.force -stiff * bin_res.bincenter, 'LineWidth', 2);
ylabel('({\itF_x/m})^{(BEEPSIS)} -  ({\itF_x/m})^{(Linear)}]','Interpreter','tex')
xlabel('{\itx} [\mum]','Interpreter','tex')

%% define force function and use BEEPSIS to also obtain damping and relative temperature
% We will consider a combination linear harmonic restoring force combined
% with softening cubic non-linear term (Duffing oscillator).
% We will search for frequency of Harmonic trap in kHz (note that mass
%   normalized stiffness is (2*pi*f)^2 and the strength of duffing
%   weakening in [\mum^{-2}]. We use these units in order to avoid
%   big difference in orders of the searched parameters which may influence some
%   of the minimization algorithms.

%force function
ffun  = @(p, x) -(2*pi*p(1)*1e3)^2 * x .* (1-p(2)*1e12*x.^2);

% we use force function and non-parametric BEEPSIS results to obtain
%  initial values
init = nlinfit(bin_res.bincenter, bin_res.force, ffun, [sqrt(-stiff)/2/pi/1e3, 1]);

% BEEPSIS, 
% first forces than all parameters, for fun use fminseach and fminunc
bres = beepsis(data_filt*1e-6, dt,ffun, init, GammaTh,kT_M, '1cA|2cA');

% check results
disp(['Fharm, '])
disp([bres.force,bres.Gamma, bres.effTr, bres.Svalue/1e7])
% we caan see 
%  - line 1+3, colums 1+2 - both fminsearch and fminunct gave same trap
%    frequency and coefficient of Duffing non-linearity. Gamma and effective
%    temperature were fixed on initial values. Value of minimized
%    function is the same in both cases 
% - line 2 (fminsearch, minimization of all four parameters). Damping decreased,
%    effective temperature increased, force remained same as on the previous
%    line. Value of minimized function decreased
% - line 4 (fminunc, minimization of all four parameters). Damping,
%   effective temperature and value of minimized function did not change.
%   This indicated that minimization failed and one should use results from
%   line 2


%% final figure plotting force profiles
iBest = 2;
figure
xd = linspace(min(data_filt),max(data_filt), 100);
plot(xd, ffun(bres.force(iBest,:),xd*1e-6),'LineWidth', 2)
hold on;
plot(bin_res.bincenter*1e6,bin_res.force, 'LineWidth', 2);
hold off
grid on
ylabel('{\itF_x/m} [ms^{-2}]','Interpreter','tex')
xlabel('{\itx} [\mum]','Interpreter','tex')
legend('BEEPSIS', 'non-parm. BEEPSIS')
