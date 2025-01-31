% Coded up by Mustafa Yilmaz - 11/1/2024
% This code calculates the elevation angle and slant range between WiFi
% access points (APs) and satellite. To calculate these two, 'aer function'
% which is a part of satellite communication toolbox is utilized. This
% function requires two inputs: satellite and ground station data.
%
% The spreadsheet inlcuding lat/long information of 74,002 census tracks 
% was provided by Cbale Labs. These lat/long values are considered as the 
% locations of APs. The name of the file is 'apLocations'.
%
% As satellite system, the TLE file (saved as txt file) of SES-2 satellite 
% is used since CableLabs performed the simulation with this satellite. 
% TLE file is downloaded from 'www.space-track.org'.
% The mean altitude of SES-2 satellite is 35,793.5 km according to 
% https://in-the-sky.org/spacecraft.php?id=37809

clear,clc
apLocations = table2array(readtable("apLocations.csv"));    % Upload AP data from census track spreadsheet
% Form a satellite scenario with satcom toolbox.
startTime  = datetime(2024,6,4,14,27,0,TimeZone="America/New_York");
stopTime   = startTime;
sampleTime = 1;      
sc = satelliteScenario(startTime,stopTime,sampleTime);
apName = string(1:size(apLocations,1));
% Add WiFi ground stations to satellite scenario
grS = groundStation(sc,apLocations(:,1),apLocations(:,2),Name=apName);
% Add satellite data to satellite scenario
sat = satellite(sc,"ses2.txt");
% When the number of ground stations or satellites is high, Matlab throws
% memory error. Therefore, it is better to process all ground stations or
% satellites part by part and save the elevation angle and slant range
% information of all links to use them in another piece of code to
% calculate the propagation, clutter, atmospheric and building entry loss.
intervals = 0:1e4:size(apLocations,1);
intervals = [intervals,size(apLocations,1)];
allAPs = [];
allTic = tic;
for i = 1:length(intervals)-1
    startT = tic;
    [~,el,range] = aer(grS(1+intervals(i):intervals(i+1)),sat);
    allAPs = [allAPs;[el,range]];
    stopT = toc(startT);
    fprintf('Elapsed Time for %d.Interval: %d\n',i,stopT)
end
save("er_APtoSat.mat","allAPs")
allToc = toc(allTic);
fprintf('Total Elapsed Time: %d\n',allToc)

%% SECOND PART OF THE CODE
%% This code calculates 
%   * Atmospheric loss (loss_a) with P.676
%   * Building entry loss (loss_bel) with P.2109
%   * Clutter loss (loss_c) with P.2108
%   * Propagation loss (PL) with free-space equation
% Finally INR is calculated and 'INR vs. CDF' is plotted
% The parameter values are taken from Cable Labs slides (CL 7 GHz NTIA
% v0.8). Next to the definition of each parameter, slide number where the
% parameter values are obtained from is provided.

% load("er_APtoSat.mat")                          % Once this data is saved above, the first part of the code can be commented out.
f_GHz   = 8.140;                                  % Center frequency [GHz] - Based on Roy's answer via email
rx_antG = 26.1;                                   % Satellite antenna gain [dB] - [Slide 28]
hGS     = randi([3 9],1,size(allAPs,1))/1e3;      % Ground terminal height [km]
hS      = 100;                                    % Top of atmosphere [km]
beta1   = (pi/180)*allAPs(:,1);                   % Elevation angle referenced to zenith [radian]
lambda  = physconst('LightSpeed')/(f_GHz*1e9);    % Wavelength [m]
type    = 1;                                      % 1 for traditional, 2 for thermally efficient construction
% Randomly select indoor APs and assign EIRP based on EIRP table - [Slide 17]
numAPs   = 47e6;                                  % Number of APs in Cable Labs simulation
satBands_MHz = 8122:8158;                         % The frequency bands used by SES2 satellite [MHz]
psdLimit = [8,23];                                % PSD limits [indoor outdoor] APs - [Slide 25]
numAP_indoor = numAPs*.98;                        % Number of indoor APs - 98% of the total APs - [Slide 25]
[eirps_data_indoor,eirps_beacon_indoor] = ueSelection(numAP_indoor,satBands_MHz,psdLimit(1));   % EIRPs of indoor APs
numAP_outdoor = numAPs*.02;                       % Number of outdoor APs - 2% of the total APs - [Slide 25]
[eirps_data_outdoor,eirps_beacon_outdoor] = ueSelection(numAP_outdoor,satBands_MHz,psdLimit(2));   % EIRPs of outdoor APs
all_eirps = 10*log10([eirps_data_indoor eirps_beacon_indoor eirps_data_outdoor eirps_beacon_outdoor]);
numAP_all = length(all_eirps);                    % Total number of indoor and outdoor APs selected 
monteC    = 1000;                                 % Number of Monte Carlo iterations for CL and BEL
iStart    = tic;
I_aggregate = zeros(100,monteC);
for i = 1:10
    rng("shuffle")
    selectAPs = randi([1 size(allAPs,1)],1,numAP_all);% Randomly select 'numAP_all' number of indoor and outdoor APs    
    rx_dBm    = zeros(numAP_all,monteC);
    Start     = tic;
    for ue = 1:numAP_all
        ap = tic;
        loss_a   = p676_AL(f_GHz,hGS(selectAPs(ue)),hS,beta1(selectAPs(ue)));
        loss_fs  = 10*log10((4*pi*allAPs(selectAPs(ue),2)/lambda)^2);
        prct     = rand(1,monteC)*100;                % Percentage of locations [%]
        prob     = rand(1,monteC)*100;                % Probability that loss is not exceeded
        loss_bel = zeros(1,monteC);
        loss_cl  = zeros(1,monteC);
        if ue < length([eirps_data_indoor eirps_beacon_indoor])+1
            for bel = 1:monteC
                loss_bel(bel) = p2109_BEL(f_GHz,allAPs(selectAPs(ue),1),prob(bel),type);
            end
        end
        for cl = 1:monteC
            loss_cl(cl)  = p2108_CL(f_GHz,allAPs(selectAPs(ue),1),prct(cl));
        end
        % Polarization mismatch loss of 3dB will be added according to
        % CableLabs presentation
        rx_dBm(ue,:) = all_eirps(ue) + rx_antG - loss_fs - loss_a - loss_bel - loss_cl - 3;
        apStop = toc(ap);
        % fprintf('Elapsed Time: %d - (%d)\n',apStop,ue)
    end
    Stop = toc(Start);
    fprintf('Total Elapsed Time - %d.Iteration: %d\n',i,Stop)
    % Aggregate interference calculation
    I_aggregate(i,:) = 10*log10(abs(sum(10.^(rx_dBm/10))));
end
iStop = toc(iStart);
fprintf('Total Elapsed Time: %d\n',iStop)

%% Noise calculation
BW     = length(satBands_MHz)*1e6-1;                   % Transponder bandwidth [Hz]
T      = 257;                                     % Receiver noise temperature of IntelSat 35e [Kelvin]
NF     = 0;                                       % Noise figure [dBm]
noiseP = 10*log10(physconst('Boltzmann')*T*BW) + NF + 30;   % dBm
INR    = I_aggregate - noiseP;                    % Interference-to-Noise Ratio

%% Plot CDF of INR
h = figure;
[bins,edges] = histcounts(reshape(INR,1,[]),'Normalization','cdf');
plot(edges,[0 bins],'b','LineWidth',2);
xlabel('INR (dB)'), ylabel('CDF'),grid on,ylim([0,1])

%% UE selection algorithm
% This algorithm assumes that data and beacon signals are transmitted from
% the same channel used by satellite, i.e., co-channel interference is
% considered for all transmissions.
% For WiFi AP case, the satellite is considered to use 8122-8158 MHz.
% The parameter values are taken from Cable Labs slides (CL 7 GHz NTIA
% v0.8). Next to the definition of each parameter, slide number where the
% parameter values are obtained from is provided.
% Total number of randomly selected UEs will be determined based on the
% Cable Labs assumption that 47M UEs are distributed in CONUS.

function [allEirps_data,allEirps_beacon] = ueSelection(numUE,satBands_MHz,psdLimit)
satBW_MHz         = length(satBands_MHz)-1;                     % BW of satellite
cbw               = [20,40,80,160,320];                         % Channel bandwidth - [Slide 17]
num_bands         = [63 32 16 8 4];                             % Number of bands for each CBW - [Slide 27]
cbw_percentage    = [.1 .1 .3 .3 .2];                           % The percentage of each CBW usage - [Slide 19]
data_percentage   = [.71 .36 .18 .09 .09]/100;                  % The percentage of data transmission of each CBW - [Slide 19]
beacon_percentage = [.9929 .9964 .9982 .9991 .9991]*.0058;      % The percentage of beacon transmission of each CBW - [Slide 19]
%% Calculate the number of UEs transmitting data and beacon for each CBW case
numUE_perCbw      = round(numUE*cbw_percentage);                % # of UE in each CBW case
numUE_perBand     = round(numUE_perCbw./num_bands);             % # of UE per band in each CBW case
numUE_perData     = round(numUE_perBand.*data_percentage);      % # of UE transmitting data per CBW
numUE_perBeacon   = round(numUE_perBand.*beacon_percentage);    % # of UE transmitting beacon per CBW
%% Assign EIRPs
eirp_values       = [1000 250 100 50 13 1];                     % Max EIRP values [mW]- [Slide 17]
if psdLimit == 23                                               % Increase EIRP values by 6dB if UEs are outdoor - [Slide 25]
    eirp_values   = 10.^((10*log10([1000 250 100 50 13 1])+6)/10);
end
eirp_percentage   = [.01 .12 .06 .19 .54 .08];                  % The percentage of usage of max EIRP values - [Slide 17]
numUE_perEirpDtemp= round(numUE_perData'.*eirp_percentage);     % Number of UE for each EIRP percentage transmitting data 
numUE_perEirpBeac = round(numUE_perBeacon'.*eirp_percentage);   % Number of UE for each EIRP percentage transmitting beacon
maxEirp_beacon    = 10^((psdLimit + 10*log10(20))/10);          % Max EIRP of beacon signal. PSD is from [Slide 25]
if satBW_MHz < 20
    maxEirp_beacon= 10^((psdLimit + 10*log10(20) - 10*log10(20/satBW_MHz))/10);
end
% Adjust beacon EIRPs based on the bandwidth of beacon [20 MHz fixed]
eirp_values_upd = eirp_values;
eirp_values_upd(eirp_values_upd > maxEirp_beacon) = maxEirp_beacon;
eirp_table_beacon = eirp_values_upd;
eirp_table_beacon = repmat(eirp_table_beacon,size(numUE_perEirpBeac,1),1);
% Adjust data EIRPs based on CBWs and satellite BW
eirp_table_data   = [];
numUE_perEirpData = [];
for i = 1:length(cbw)
    bands = ueBands(cbw(i));
    for b = 1:size(bands,1)
        maxEirp_data = maxEirp(bands(b,:),satBands_MHz,psdLimit);
        eirp_values_upd = eirp_values;
        eirp_values_upd(eirp_values_upd > maxEirp_data) = maxEirp_data;
        eirp_table_data = [eirp_table_data;eirp_values_upd];
        numUE_perEirpData = [numUE_perEirpData;numUE_perEirpDtemp(i,:)]; % Increase the number of UEs by the number of overlapping bands with satellite bands
    end
end
% Multiplex the EIRP values as many as the number of UEs in 'numUE_perEirp..'
% variable and put them in a variable orderly.
allEirps_data   = [];
for i = 1:size(numUE_perEirpData,1)
    for k = 1:size(numUE_perEirpData,2)
        allEirps_data   = [allEirps_data,repmat(eirp_table_data(i,k),1,numUE_perEirpData(i,k))];
    end
end
allEirps_beacon = [];
for i = 1:size(numUE_perEirpBeac,1)
    for k = 1:size(numUE_perEirpBeac,2)
        allEirps_beacon = [allEirps_beacon,repmat(eirp_table_beacon(i,k),1,numUE_perEirpBeac(i,k))];
    end
end
end

%% Calculate max EIRP for data and beacon signals
function maxEirp_data = maxEirp(bands,satBands_MHz,psdLimit)
for i = 1:size(bands,1)
    overlap_bw = length(intersect(satBands_MHz,bands(i,:)))-1;      % Find overlapping portion of UE bands with satellite bands
    maxEirp_data = 10^((psdLimit + 10*log10(overlap_bw))/10);
end
end
%% UE bands overlapping with satellite bands
function bands = ueBands(cbw)
if cbw == 20
    bands = [8105:8125; 8125:8145; 8145:8165];
elseif cbw == 40
    bands = [8105:8145; 8145:8185];
elseif cbw == 80
    bands = 8105:8185;
elseif cbw == 160
    bands = 8025:8185;
elseif cbw == 320
    bands = 7865:8185;
end
end

%% ITU-R P.2108 Implementation
% Coded up by Mustafa Yilmaz - 10/3/2024
% Clutter loss as defined in ITU P.2108 Section 3.3.2, Equation-7
% Input parameters
% theta is the elevation angle which is the angle of the airborne platform 
% or satellite as seen from the terminal - 0 to 90 [degree]
% f_GHz is the center frequency - 10 to 100 [GHz]
% prct is percentage of locations - 0 to 100 [%]
function loss_c = p2108_CL(f_GHz,theta,prct)        % [dB]
f       = [10,15];          % Frequency values to use for extrapolation [GHz]
A_1     = 0.05;
K_1     = 93*f.^0.175;
l_p2108 = (-K_1.*(log(1-prct/100)).*cot(A_1*(1-theta/90)+pi*theta/180)).^(0.5*(90-theta)/90)-1+0.6*norminv(prct/100);
% This model is for f>10GHz. To get exact value for 8 GHz, we need to use extrapolation
slope = (l_p2108(2)-l_p2108(1))/(f(2)-f(1));
loss_c = l_p2108(2)+slope*(f_GHz - f(2));
end

%% ITU-R P.2109 Implementation
% Coded up by Mustafa Yilmaz - 10/3/2024
% Building entry loss as defined in ITU P.2109 Section 3, Equation-1
% Input parameters
% theta is the elevation angle of the path at the building façade [degrees]
% f_GHz is the center frequency [GHz]
% prob is the probability that loss is not exceeded (0.0 ≤ P ≤ 1.0)
% type - 1 for traditional, 2 for thermally efficient
function loss_bel = p2109_BEL(f_GHz,theta,prob,type)        % [dB]
r      = [12.64 28.19];
s      = [3.72 -3];
t      = [0.96 8.48];
u      = [9.6 13.5];
v      = [2 3.8];
w      = [9.1 27.8];
x      = [-3 -2.9];
y      = [4.5 9.4];
z      = [-2 -2.1];
sigma1 = u(type) + v(type)*log10(f_GHz);
sigma2 = y(type) + z(type)*log10(f_GHz);
mu2    = w(type) + x(type)*log10(f_GHz);
L_e    = 0.212*abs(theta); % Correction for elevation angle of the path at the building facade                                          
L_h    = r(type) + s(type)*log10(f_GHz) + t(type)*log10(f_GHz)^2; % Median loss for horizontal paths
mu1    = L_h + L_e;
C      = -3;
A_P    = sigma1*norminv(prob/100) + mu1;
B_P    = sigma2*norminv(prob/100) + mu2;
loss_bel = 10*log10(10^(0.1*A_P)+10^(0.1*B_P)+10^(0.1*C));
end

%% ITU-R P.676 Implementation
% Code was downloaded from https://github.com/NTIA/propagation/wiki/P676 - 10/3/2024
% Atmospheric loss as defined in ITU P.676
% This code allows the user to select from any of the six reference 
% atmospheres defined in Recommendation ITU-R P.835 (Reference standard atmospheres).
% Input parameters
% beta1 is elevation angle referenced to zenith [radian]
% h1 is the ground terminal height [km]
% h2 is the top of atmosphere [km]
% f_GHz is the center frequency [GHz]
% P676ReferenceAtmosphere.MAGRA refers to global average atmospheric model
% Other atmospheric models:
%  HighLatitudeSummer,HighLatitudeWinter,MidLatitudeSummer,MidLatitudeWinter,LowLatitude
function loss_a = p676_AL(f_GHz,h1,h2,beta1)        % [dB]
atmosphere = P676ReferenceAtmosphere.MAGRA;
[~, result] = P676.SlantPathAttenuation(f_GHz, h1, h2, beta1, atmosphere);
loss_a = result.A_gas__db;
end