function data = Preprocessing(mouseName,garrFile,intanFolder)
% Preprocessing(mouseName,garrFile,intanFolder)
% reads garr and intan and saves them as useful variables aligned with the
% right time stamps
%   mouseName: name of mouse (string)
%   garrFile: name of file where garr mat data is stored (string)
%   intanFolder: name of file where intan data is stored (string)
% RL 11/2024

%% load variables

% mouseName = "MS13B";
% garrFile = 'dis_vel_rew_2024-7-25_19_14_21_MS13B_REC_BASE.mat';
% intanFolder = "W:\mouse_chemogenetics\MS13B\MS13B_240725_175414";
% mouseName = "MS14A";
% garrFile = 'dis_vel_rew_2024-9-6_18_11_35_MS14_REC.mat';
% intanFolder = "W:\mouse_chemogenetics\MS14A\MS14A_240906_162727";

mkdir Preprocess_figures % To store saved figures

load(garrFile, 'garr', 'adata');

if garr(end,1) == 0
    garr(end,:) = [];
end

% Parameters
samplerate = 30;
nBin = 100;

% Extract date from file name
[~,fdate] = fileparts(garrFile);
fdate = strtok(fdate(13:22), '_');

environment = garr(:,7);
%envnum = environment(floor(end/2));

% Only one environment for now
[envs,~,ic] = unique(environment);
for ii = 1:length(envs)
    [~,VRradius(ii)] = findEnvironment(envs(ii));
end
vrrad_eachbin = VRradius(ic);  % radius of environment at each time point, seems to be very accurate even at env transitions

pos_raw = garr(:,5); % position in laps (based on VR)
poscm_raw = pos_raw*VRradius*2*pi;
pos = mod(pos_raw,1);
poscm = pos*VRradius*2*pi;
posBin = floor(pos*nBin) + 1;

samples = [0:size(garr,1)-1]';
vel = [0; diff(poscm_raw)] * samplerate; % Velocity (cm/s)
ts = samples / samplerate; % Timestamp (s)
lap = ones(length(samples),1);
reward = [0; diff(garr(:,3)) < 0]; % Reward
licks = garr(:,4); % Lick
objectChanges = garr(:,6);

% Smooth velocity
vel(vel > 120) = 0; % Remove unusually high velocity values
vel_smoothed = smoothdata(vel,"gaussian",3);

% Redundant variables (for now)
% pos_enc = garr(:,1); % position in revolutions of wheel (based on rotary encoder)
% pos_incm_enc = pos_enc.*12.*pi; % position in cm (based on rotary encoder) -- in movie mode this will be different from pos
% inst_vel = garr(:,2);  %cm/s (based on VR)
% vel_enc = [0;diff(pos_incm_enc)];   %cm/sample (based on encoder)
% vel_enc = vel_enc.*samplerate; %cm/s (based on encoder)

% Split dataset into laps
lapStart_idx = [1; find(diff(pos) < -0.95) + 1]; % Threshold difference at 95% of track distance
lapEnd_idx = find(diff(pos) < -0.95);
nLap = length(lapEnd_idx);
lapStart_idx = lapStart_idx(1:nLap);

for lap_no = 1:nLap
    lap(lapStart_idx(lap_no):lapEnd_idx(lap_no)) = lap_no;
end
lap(lapEnd_idx(end)+1:end) = NaN;

dataLap_all = [lap posBin ts pos poscm vel vel_smoothed reward licks objectChanges environment];

% Need movement data only (zero velocity)
posChange = (diff(dataLap_all(:,4)) ~= 0);
% nonzeroVel_idx = logical([1; posChange]);
nonzeroVel_idx = logical([posChange; 0]);
dataLap = dataLap_all(nonzeroVel_idx,:);

% Align timestamps with Intan
ori = pwd;
cd(intanFolder);
digIn = readNPY(mouseName + "-digIn.npy");
digInts = readNPY(mouseName + "-digInts.npy");
% lfp = readNPY(mouseName + "-lfp.npy");
% lfpts = readNPY(mouseName + "-lfpts.npy");
intants = loadintants(mouseName + "-digIn.npy");
cd(ori);

figure;
subplot(2,1,1)
plot(diff(intants),'b');
hold on
plot(diff(dataLap(:,3)),'r');
hold off
title("Timestamp Alignment (Before)");
xlabel("Timestamps"); ylabel("Magnitude of Difference (s)");
subplot(2,1,2)
% check timestamp alignment
dt = diff(dataLap(:,3));
dint = diff(intants);
[c,lag] = xcorr(dint,dt);
[~,ix] = max(c);
offset = lag(ix);
% [~,idx1] = max(diff(intants));
% [~,idx2] = max(diff(dataLap(:,3)));
% offset = idx1 - idx2;
if offset > 0 % Intan starts first
    time_offset = intants(abs(offset));
    offset_shift = find(intants > time_offset,1,'first');
    % intants = intants(abs(offset)+1:end);
    intants = intants(offset_shift:end);
    
    % offset_shift = find(lfpts > time_offset,1,'first');
    % lfpts = lfpts(offset_shift:end);
    % lfp = lfp(:,offset_shift:end);
else          % Smoothwalk starts first
    time_offset = dataLap(abs(offset),3);
    offset_shift = find(dataLap(:,3) > time_offset,1,'first');
    % dataLap = dataLap(abs(offset)+1:end,:);
    dataLap = dataLap(offset_shift:end,:);
end
hold all
plot(diff(intants),'b');
plot(diff(dataLap(:,3)),'r');
hold off
title("Timestamp Alignment (After)");
xlabel("Timestamps"); ylabel("Magnitude of Difference (s)");
saveas(gcf,"Preprocess_figures/" + mouseName +  "_Timestamp_Alignment.png")

dataLap(:,3) = intants;

% Remove fake laps (due to jittering of mice position around lap transition point)
fakeLaps = unique(dataLap(find(diff(dataLap(:,2)) > 90),1));
fakeLaps = [dataLap(1,1); fakeLaps]; % Also remove the first incomplete lap
fakeLaps_idx = zeros(size(dataLap,1),1);
fakeLapNo_correction = dataLap(:,1);
for i = 1:length(fakeLaps)
    temp_idx = find(dataLap(:,1) == fakeLaps(i));
    fakeLaps_idx(temp_idx) = 1;
    fakeLapNo_correction(temp_idx(1):end) = fakeLapNo_correction(temp_idx(1):end) - 1; % Fix lap no.
end
if dataLap(1,1) > 1 % Start lap is not labeled as 1
    fakeLapNo_correction = fakeLapNo_correction - dataLap(1,1) + 1;
end
fakeLaps_idx(find(isnan(dataLap(:,1)))) = 1; % Remove final incomplete lap
nLap = nLap - (dataLap(1,1) - 1) - length(fakeLaps);
dataLap(:,1) = fakeLapNo_correction;
dataLap(logical(fakeLaps_idx),:) = []; % Remove fake laps + first and final incomplete lap from behaviour data

% Bin velocity by position bin per lap
dataBinned = zeros(nLap*nBin,6);
for i = 1:size(dataBinned,1)
    lap_no = floor((i-1)/nBin) + 1;
    dataBinned(i,1) = lap_no; % Lap No.
    bin_no = mod(i,nBin);
    if mod(i,nBin) == 0
        bin_no = nBin;
    end
    dataBinned(i,2) = bin_no; % Bin No.
    temp_array = dataLap((dataLap(:,1) == lap_no) & (dataLap(:,2) == bin_no),:);
    if ~isempty(temp_array)
        dataBinned(i,3) = temp_array(1,3); % Start time
        dataBinned(i,4) = temp_array(end,3); % End time
        dataBinned(i,5) = temp_array(end,3) - temp_array(1,3) + (1/samplerate); % Duration
        dataBinned(i,6) = mean(temp_array(:,7),'omitnan'); % Use smoothed velocity
    end
end

vel_binned = reshape(dataBinned(:,6),[nBin,nLap])';

% Extract reward signals from intan
digInts_arr = digInts(find(digInts >= dataLap(1,3),1,'first'):find(digInts <= dataLap(end,3),1,'last'));
digIn_arr = digIn(2,find(digInts >= dataLap(1,3),1,'first'):find(digInts <= dataLap(end,3),1,'last'));
digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
intan_pos = interp1(dataLap(:,3),dataLap(:,5),digInts_arr);

reward_ts = dataLap(find(dataLap(:,8) == 1),3);
reward_pos = dataLap(find(dataLap(:,8) == 1),5) / (VRradius*2*pi);
intan_reward_ts = digInts_arr(find(digIn_arr_ == 1));
intan_reward_pos = intan_pos(find(digIn_arr_ == 1)) / (VRradius*2*pi);

% figure;
% % subplot(3,1,1)
% % hold on
% % yyaxis left
% % % plot(temp.garr(1:end-1,8),temp.garr(1:end-1,1))
% % plot(timestamp,pos)
% % title("Smoothwalk"); ylabel("Distance (cm)");
% % yline(0.83*VRradius*2*pi,'--');
% % yyaxis right
% % plot(timestamp,reward,'r')
% % hold off
% subplot(2,1,1)
% hold on
% yyaxis left
% plot(digInts_arr,intan_pos)
% title("Intan"); ylabel("Distance (cm)");
% yline(0.83*VRradius*2*pi,'--');
% yyaxis right
% plot(digInts_arr,digIn_arr_,'r')
% hold off
% subplot(2,1,2)
% plot(dataLap(:,3),dataLap(:,6))
% title("Velocity (cm/s)"); ylabel("Time (s)");
% ylim([-10,120])

figure;
subplot(2,1,1)
scatter(reward_ts,reward_pos,'bx'); yline(0.83,'k--');
% mean(reward_pos)
title("Smoothwalk"); xlabel("Time (s)"); ylabel("Reward Position");
ylim([0 1]);
subplot(2,1,2)
scatter(intan_reward_ts,intan_reward_pos,'bx'); yline(0.83,'k--');
% mean(intan_reward_pos)
title("Intan"); xlabel("Time (s)"); ylabel("Reward Position");
ylim([0 1]);
saveas(gcf,"Preprocess_figures/" + mouseName +  "_Lap_Reward_Position.png")

% figure;
% yyaxis left
% plot(dataLap(:,3), dataLap(:,8))
% yyaxis right
% plot(digInts_arr, digIn_arr_)

% start_trial = 1; % 310 - 360
% end_trial = 716;
% test_lap = dataLap(find(dataLap(:,1) == start_trial,1,'first'):find(dataLap(:,1) == end_trial,1,'last'),:);
% timestamp = test_lap(:,3);
% reward = test_lap(:,8);
% vel = test_lap(:,6);
% pos = test_lap(:,5);
% digInts_arr = digInts(find(digInts >= test_lap(1,3),1,'first'):find(digInts <= test_lap(end,3),1,'last'));
% digIn_arr = digIn(2,find(digInts >= test_lap(1,3),1,'first'):find(digInts <= test_lap(end,3),1,'last'));
% digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
% 
% intan_pos = interp1(timestamp,pos,digInts_arr);
% 
% reward_ts = timestamp(find(reward == 1));
% intan_reward_ts = digInts_arr(find(digIn_arr_ == 1));

% Extract lick signals from intan


data.dataLap_all = dataLap_all;
data.dataLap = dataLap;
data.dataBinned = dataBinned;
data.vel_binned = vel_binned;
data.nonzeroVel_idx = nonzeroVel_idx;
data.date = fdate;
% data.lfp = lfp;
% data.lfpts = lfpts;
data.offset = offset;
data.time_offset = time_offset;
data.reward = [intan_reward_ts intan_reward_pos];
save("sess.mat","data",'-v7.3');

% Plot behavioural data
figure;
set(gcf, 'Position',  [100, 100, 1500, 500])
hold all
% First 10 laps
start_lap = 1;
end_lap = 10;
plot_start_idx = find(dataLap(:,1) == start_lap,1,'first'); 
plot_end_idx = find(dataLap(:,1) == end_lap,1,'last');
subplot(2,2,1)
yyaxis left
plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,5))
title("Distance (cm)");
xlabel('Time (s)'); ylabel('Distance (cm)');
yline(0.83*VRradius*2*pi,'--');
yyaxis right
digInts_arr = digInts(find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
digIn_arr = digIn(2,find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
plot(digInts_arr,digIn_arr_,'r')
ylabel('Reward');
subplot(2,2,3)
plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,6))
title("Velocity (cm/s)");
xlabel('Time (s)'); ylabel('Velocity (cm/s)');
ylim([0,120])
hold off
hold all
% Last 10 laps
start_lap = dataLap(end,1) - 9;
end_lap = dataLap(end,1);
plot_start_idx = find(dataLap(:,1) == start_lap,1,'first'); 
plot_end_idx = find(dataLap(:,1) == end_lap,1,'last');
subplot(2,2,2)
yyaxis left
plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,5))
title("Distance (cm)");
xlabel('Time (s)'); ylabel('Distance (cm)');
yline(0.83*VRradius*2*pi,'--');
yyaxis right
digInts_arr = digInts(find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
digIn_arr = digIn(2,find(digInts >= dataLap(plot_start_idx,3),1,'first'):find(digInts <= dataLap(plot_end_idx,3),1,'last'));
digIn_arr_ = [0 (abs(diff(int8(digIn_arr))))];
plot(digInts_arr,digIn_arr_,'r')
ylabel('Reward');
subplot(2,2,4)
plot(dataLap(plot_start_idx:plot_end_idx,3),dataLap(plot_start_idx:plot_end_idx,6))
title("Velocity (cm/s)");
xlabel('Time (s)'); ylabel('Velocity (cm/s)');
ylim([0,120])
hold off
saveas(gcf,"Preprocess_figures/" + mouseName +  "_Distance_Velocity.png")

% Plot position-binned data
figure;
imagesc(vel_binned); cbar = colorbar; cbar.Label.String = 'Velocity (cm/s)';
title("Velocity (cm/s)");
xlabel('Position (Bins)'); ylabel('Lap');
saveas(gcf,"Preprocess_figures/" + mouseName +  "_Velocity_Binned.png")

figure;
xdat = (1:nBin)*0.01*VRradius*2*pi;
ydat = mean(vel_binned,'omitnan');
y_std = std(vel_binned,'omitnan');
plot(xdat,ydat,'k-','LineWidth', 2);
hold on
curve1 = ydat + y_std;
curve2 = ydat - y_std;
x2 = [xdat, fliplr(xdat)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.3);
title(mouseName);
xlabel('Position (cm)'); ylabel('Velocity (cm/s)');
hold off
saveas(gcf,"Preprocess_figures/" + mouseName +  "_Velocity.png")

% figure;
% dur_binned = reshape(dataBinned(:,5),[nBin,nLap])';
% imagesc(dur_binned); colorbar;














