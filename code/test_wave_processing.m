%% load data
file_number = 20;
load(['C:/Users/hjohn/Documents/work/utqiagvik_mooring/Matlab_Format_with_coord_transforms/S106174A002_let_s_go_' num2str(file_number) '.mat']);

% set plotting limits
datenum_lims = [min(Data.Average_Time) max(Data.Average_Time)];

% Find indices for each burst
firstindex = find(Data.Burst_EnsembleCount==1);
spacing = mean(diff(firstindex));
if isnan(spacing) % exception for only one burst in file
    spacing=length(Data.Burst_EnsembleCount)-firstindex; 
end 
nb = length(firstindex);

% initialize burst counter
bcounter = 1;

% initialize structure
sigBurst = struct;

% Loop through bursts and process
for b = 1:nb

    % select burst data
    burstInd = firstindex(b) + [1:spacing] - 1;
    burstInd( burstInd > length(Data.Burst_Time) ) = [];

    burst_sampling_rate = double(Config.Burst_SamplingRate);
    burst_pressure = Data.Burst_Pressure(burstInd) - 40;
    
    bursttime = Data.Burst_Time(burstInd);

    sigBurst(bcounter).time = bursttime(1);

    [ Hs, Tp, Hig, Tig, E, f ] = Pwaves(burst_pressure,burst_sampling_rate);
    sigBurst(bcounter).sigwaveheight = Hs;
    sigBurst(bcounter).peakwaveperiod = Tp;

    % increment burst counter;
    bcounter = bcounter + 1;
end

%% create plots
% plot significant wave height from pressure
figure(3);
plot([sigBurst.time],[sigBurst.sigwaveheight]);
ylabel('Significant wave height (m)');
datetick('x');

% plot peak wave period from pressure
figure(2);
plot([sigBurst.time],[sigBurst.peakwaveperiod]);
ylabel('Peak wave period (s)');
datetick('x');