%% 
%clear all; close all; clc;
%%
colorO = [0.69, 0.94, 1];
colorI = [1, 0.69, 0.69];
edgeC=[0, 0, 0];
colors = [0.6196 0.0039 0.2588 ; 0.9353 0.2431 0.3098 ; 0.9569 0.4275 0.2627 ; 0.9922 0.6824 0.3804 ; 0.9961 0.8784 0.5451 ; 0.9961 0.8784 0.5451 ; 0.6706 0.8667 0.6431 ; 0.4000 0.7608 0.6471 ; 0.1961 0.5333 0.7412 ; 0.3686 0.3098 0.6353];

%% Intensity changes compared with previous contact angles (not working well)

intensityChange = [];

for i = 1:length(contactAngles)
    if i == 1,
        intensityChange(i) = 0;
    else
    intensityChange(i) = histStorage(i) - histStorage(i-1);
    end
end

figure; hold on;
plot(intensityChange,contactAngles);

% [X,Y] = meshgrid(1:1:360,contactAngles);
% s = pcolor(X(:,1:360),Y(:,1:180),intensityChange(:,1:180));

titlestr = ['Intensity changes compared with previous contact angles'];
title(titlestr)
xlabel('Intensity changes(a.u.)')
ylabel('Contact Angles(\Theta\circ)')
%zlabel('Angle of Bin Rays(\circ)')
%colorbar
%colormap(viridis)

%% Ratiometric bottom readout vs contact angles written by katrina copyright 2023
    
intensityZero = [];
intensity45 = [];
ratiomatricIntensity = [];

%0&45
for i = 1:length(contactAngles)
    %intensityZero(i) = mean(intHistStorage(i,265:275))/mean(intHistStorage(1,260:280));
    intensityZero(i) = mean(histStorage(i,265:275));
    %intensity45(i) = mean(intHistStorage(i,310:320))/mean(intHistStorage(1,305:325)); 
    intensity45(i) = mean(histStorage(i,310:320));
end
figure;
subplot(2,1,1); hold on;
%subplot(1,1,1); hold on;
titlestr = 'Degrees of bottom readout intensity';
title(titlestr)
%titlestr = ['Degrees of bottom readout']
plot(contactAngles,intensityZero,'DisplayName','0-degrees-intensity')
plot(contactAngles,intensity45,'DisplayName','45-degrees-intensity')

xlabel('Contact Angles(\Theta\circ)')
ylabel('Intensity(a.u.)')

%ylim([0.75 1.5])
legend

%0/45
ratiomatricIntensity = intensityZero./intensity45;
%xlim([135 180])
subplot(2,1,2); hold on;
%subplot(1,1,1);
titlestr = 'Ratiometric bottom readout';
title(titlestr)
%titlestr = ['Ratiometric bottom readout']
plot(contactAngles,ratiomatricIntensity/ratiomatricIntensity(1))
xlabel('Contact Angles(\Theta\circ)')
%ylim([0.75 1.5])
ylabel('Ratiometric 0\circ/45\circ intensity (a.u.)')

%% bottom readout vs contact angles written by katrina copyright 2025
    
intensityZero = [];
intensity15 = [];
intensity30 = [];
intensity45 = [];
intensity60 = [];
intensity75 = [];
intensity90 = [];
ratiomatricIntensity = [];

%0&45
for i = 1:length(contactAngles)
%     %intensityZero(i) = mean(intHistStorage(i,265:275))/mean(intHistStorage(1,260:280));
%     intensityZero(i) = mean(histStorage(i,268:272)); %270
%     intensity15(i) = mean(histStorage(i,283:287)); %285
%     intensity30(i) = mean(histStorage(i,298:302)); %300
%     
%     %intensity45(i) = mean(intHistStorage(i,310:320))/mean(intHistStorage(1,305:325)); 
%     intensity45(i) = mean(histStorage(i,313:317)); %315
%     intensity60(i) = mean(histStorage(i,328:332)); %330
%     intensity75(i) = mean(histStorage(i,343:347)); %345
%     intensity90(i) = mean(histStorage(i,358:360)); %360
    
    %intensityZero(i) = mean(intHistStorage(i,265:275))/mean(intHistStorage(1,260:280));
    intensityZero(i) = mean(histStorage(i,269:271)); %270
    intensity15(i) = mean(histStorage(i,284:286)); %285
    intensity30(i) = mean(histStorage(i,299:301)); %300
    
    %intensity45(i) = mean(intHistStorage(i,310:320))/mean(intHistStorage(1,305:325)); 
    intensity45(i) = mean(histStorage(i,314:316)); %315
    intensity60(i) = mean(histStorage(i,329:331)); %330
    intensity75(i) = mean(histStorage(i,344:346)); %345
    intensity90(i) = mean(histStorage(i,359:360)); %360
end
figure;
subplot(2,1,1); hold on;
%subplot(1,1,1); hold on;
titlestr = 'Degrees of bottom readout intensity';
title(titlestr)
%titlestr = ['Degrees of bottom readout']
plot(contactAngles,intensityZero,'DisplayName','0-degrees-intensity')
plot(contactAngles,intensity15,'DisplayName','15-degrees-intensity')
plot(contactAngles,intensity30,'DisplayName','30-degrees-intensity')
plot(contactAngles,intensity45,'DisplayName','45-degrees-intensity')
plot(contactAngles,intensity60,'DisplayName','60-degrees-intensity')
plot(contactAngles,intensity75,'DisplayName','75-degrees-intensity')
plot(contactAngles,intensity90,'DisplayName','90-degrees-intensity')

xlabel('Contact Angles(\Theta\circ)')
ylabel('Intensity(a.u.)')

%ylim([0.75 1.5])
legend

%0/45
ratiomatricIntensity = intensityZero./intensity45;
%xlim([135 180])
subplot(2,1,2); hold on;
%subplot(1,1,1);
titlestr = 'Ratiometric bottom readout';
title(titlestr)
%titlestr = ['Ratiometric bottom readout']
plot(contactAngles,ratiomatricIntensity/ratiomatricIntensity(1))
xlabel('Contact Angles(\Theta\circ)')
%ylim([0.75 1.5])
ylabel('Ratiometric 0\circ/45\circ intensity (a.u.)')

%% Normalized bottom readout vs contact angles written by katrina copyright 2025

% Initialize intensity arrays
intensityZero = [];
intensity15 = [];
intensity30 = [];
intensity45 = [];
intensity60 = [];
intensity75 = [];
intensity90 = [];
ratiomatricIntensity = [];

% Compute intensities for different contact angles
for i = 1:length(contactAngles)
    intensityZero(i) = mean(histStorage(i,269:271)); %270
    intensity15(i) = mean(histStorage(i,284:286)); %285
    intensity30(i) = mean(histStorage(i,299:301)); %300
    intensity45(i) = mean(histStorage(i,314:316)); %315
    intensity60(i) = mean(histStorage(i,329:331)); %330
    intensity75(i) = mean(histStorage(i,344:346)); %345
    intensity90(i) = mean(histStorage(i,359:360)); %360
end

% Normalize intensities by their first values
intensityZeroNorm = intensityZero / intensityZero(1); % or normalize by end value intensityZeroNorm = intensityZero / intensityZero(end);
intensity15Norm = intensity15 / intensity15(1);
intensity30Norm = intensity30 / intensity30(1);
intensity45Norm = intensity45 / intensity45(1);
intensity60Norm = intensity60 / intensity60(1);
intensity75Norm = intensity75 / intensity75(1);
intensity90Norm = intensity90 / intensity90(1);

% Create figure and plot normalized intensities
figure;
subplot(2,1,1); hold on;
title('Normalized Degrees of Bottom Readout Intensity')

plot(contactAngles, intensityZeroNorm, 'DisplayName', '0°-intensity')
plot(contactAngles, intensity15Norm, 'DisplayName', '15°-intensity')
plot(contactAngles, intensity30Norm, 'DisplayName', '30°-intensity')
plot(contactAngles, intensity45Norm, 'DisplayName', '45°-intensity')
plot(contactAngles, intensity60Norm, 'DisplayName', '60°-intensity')
plot(contactAngles, intensity75Norm, 'DisplayName', '75°-intensity')
plot(contactAngles, intensity90Norm, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend

% Compute and plot normalized ratiometric intensity
ratiomatricIntensity = intensityZero ./ intensity45;
ratiomatricIntensityNorm = ratiomatricIntensity / ratiomatricIntensity(1);

subplot(2,1,2); hold on;
title('Normalized Ratiometric Bottom Readout')

plot(contactAngles, ratiomatricIntensityNorm, 'DisplayName', '0°/45° ratio')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/45° Intensity (a.u.)')
legend

%% Normalized bottom readout from CA90-180 vs contact angles written by katrina copyright 2025

% Define the range of contact angles to plot (modify these values as needed)
minAngle = 90;  % Set minimum contact angle (e.g., 90)
maxAngle = 180; % Set maximum contact angle (e.g., 180)

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Compute intensities for different contact angles
intensityZero = mean(histStorage(:, 269:271), 2); % 270
intensity15 = mean(histStorage(:, 284:286), 2); % 285
intensity30 = mean(histStorage(:, 299:301), 2); % 300
intensity45 = mean(histStorage(:, 314:316), 2); % 315
intensity60 = mean(histStorage(:, 329:331), 2); % 330
intensity75 = mean(histStorage(:, 344:346), 2); % 345
intensity90 = mean(histStorage(:, 359:360), 2); % 360

% Filter the intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity15 = intensity15(angleRange);
intensity30 = intensity30(angleRange);
intensity45 = intensity45(angleRange);
intensity60 = intensity60(angleRange);
intensity75 = intensity75(angleRange);
intensity90 = intensity90(angleRange);

% Normalize intensities by their last values in the selected range
intensityZeroNorm = intensityZero / intensityZero(1);
intensity15Norm = intensity15 / intensity15(1);
intensity30Norm = intensity30 / intensity30(1);
intensity45Norm = intensity45 / intensity45(1);
intensity60Norm = intensity60 / intensity60(1);
intensity75Norm = intensity75 / intensity75(1);
intensity90Norm = intensity90 / intensity90(1);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ')'])

plot(filteredAngles, intensityZeroNorm, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity15Norm, 'DisplayName', '15°-intensity')
plot(filteredAngles, intensity30Norm, 'DisplayName', '30°-intensity')
plot(filteredAngles, intensity45Norm, 'DisplayName', '45°-intensity')
plot(filteredAngles, intensity60Norm, 'DisplayName', '60°-intensity')
plot(filteredAngles, intensity75Norm, 'DisplayName', '75°-intensity')
plot(filteredAngles, intensity90Norm, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend

% Compute and plot normalized ratiometric intensity for the selected contact angle range
ratiomatricIntensity = intensityZero ./ intensity45;
ratiomatricIntensityNorm = ratiomatricIntensity / ratiomatricIntensity(1);

subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (CA ' num2str(minAngle) '-' num2str(maxAngle) ')'])

plot(filteredAngles, ratiomatricIntensityNorm, 'DisplayName', '0°/45° ratio')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/45° Intensity (a.u.)')
legend
grid on

%% Normalized Bottom Readout from CA90-180 vs Contact Angles - Written by Katrina, Copyright 2025, with NA=0.22 counted, +-12, normalize to CA=90 

% Define the range of contact angles to plot (modify these values as needed)
minAngle = 90;  % Set minimum contact angle (e.g., 90)
maxAngle = 180; % Set maximum contact angle (e.g., 180)

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±13 degrees around each contact angle
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute intensities for different contact angles with ±13-degree summation
intensityZero = sumIntensityRange(histStorage, 270); % 270 ± 12
intensity15 = sumIntensityRange(histStorage, 285); % 285 ± 12
intensity30 = sumIntensityRange(histStorage, 300); % 300 ± 12
intensity45 = sumIntensityRange(histStorage, 315); % 315 ± 12
intensity60 = sumIntensityRange(histStorage, 330); % 330 ± 12
intensity75 = sumIntensityRange(histStorage, 345); % 345 ± 12
intensity90 = sumIntensityRange(histStorage, 360); % 360 ± 12

% Filter the intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity15 = intensity15(angleRange);
intensity30 = intensity30(angleRange);
intensity45 = intensity45(angleRange);
intensity60 = intensity60(angleRange);
intensity75 = intensity75(angleRange);
intensity90 = intensity90(angleRange);

% Normalize intensities by their first values in the selected range
intensityZeroNorm = intensityZero / intensityZero(1);
intensity15Norm = intensity15 / intensity15(1);
intensity30Norm = intensity30 / intensity30(1);
intensity45Norm = intensity45 / intensity45(1);
intensity60Norm = intensity60 / intensity60(1);
intensity75Norm = intensity75 / intensity75(1);
intensity90Norm = intensity90 / intensity90(1);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ')'])

plot(filteredAngles, intensityZeroNorm, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity15Norm, 'DisplayName', '15°-intensity')
plot(filteredAngles, intensity30Norm, 'DisplayName', '30°-intensity')
plot(filteredAngles, intensity45Norm, 'DisplayName', '45°-intensity')
plot(filteredAngles, intensity60Norm, 'DisplayName', '60°-intensity')
plot(filteredAngles, intensity75Norm, 'DisplayName', '75°-intensity')
plot(filteredAngles, intensity90Norm, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend

% Compute and plot normalized ratiometric intensity for the selected contact angle range
ratiomatricIntensity = intensityZero ./ intensity45;
ratiomatricIntensityNorm = ratiomatricIntensity / ratiomatricIntensity(1);

subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (CA ' num2str(minAngle) '-' num2str(maxAngle) ')'])

plot(filteredAngles, ratiomatricIntensityNorm, 'DisplayName', '0°/45° ratio')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/45° Intensity (a.u.)')
legend
grid on

%% Normalized Bottom Readout from CA90-180 vs Contact Angles - Written by Katrina, Copyright 2025, with NA=0.22 counted, ±12°, normalized to CA=180

% Define the range of contact angles to plot (modify these values as needed)
minAngle = 90;  % Set minimum contact angle (e.g., 90)
maxAngle = 180; % Set maximum contact angle (e.g., 180)

% Automatically filter contact angles based on the selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around each contact angle
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute intensities for different contact angles with ±12-degree summation
intensityZero = sumIntensityRange(histStorage, 270); % 270 ± 12
intensity15 = sumIntensityRange(histStorage, 285); % 285 ± 12
intensity30 = sumIntensityRange(histStorage, 300); % 300 ± 12
intensity45 = sumIntensityRange(histStorage, 315); % 315 ± 12
intensity60 = sumIntensityRange(histStorage, 330); % 330 ± 12
intensity75 = sumIntensityRange(histStorage, 345); % 345 ± 12
intensity90 = sumIntensityRange(histStorage, 360); % 360 ± 12

% Filter the intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity15 = intensity15(angleRange);
intensity30 = intensity30(angleRange);
intensity45 = intensity45(angleRange);
intensity60 = intensity60(angleRange);
intensity75 = intensity75(angleRange);
intensity90 = intensity90(angleRange);

% Normalize intensities by their last values in the selected range (CA = 180°)
intensityZeroNorm = intensityZero / intensityZero(end);
intensity15Norm = intensity15 / intensity15(end);
intensity30Norm = intensity30 / intensity30(end);
intensity45Norm = intensity45 / intensity45(end);
intensity60Norm = intensity60 / intensity60(end);
intensity75Norm = intensity75 / intensity75(end);
intensity90Norm = intensity90 / intensity90(end);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12°)'])

plot(filteredAngles, intensityZeroNorm, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity15Norm, 'DisplayName', '15°-intensity')
plot(filteredAngles, intensity30Norm, 'DisplayName', '30°-intensity')
plot(filteredAngles, intensity45Norm, 'DisplayName', '45°-intensity')
plot(filteredAngles, intensity60Norm, 'DisplayName', '60°-intensity')
plot(filteredAngles, intensity75Norm, 'DisplayName', '75°-intensity')
plot(filteredAngles, intensity90Norm, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend
grid on

% Compute and plot normalized ratiometric intensity for the selected contact angle range
ratiomatricIntensity = intensityZero ./ intensity30;
ratiomatricIntensityNorm = ratiomatricIntensity / ratiomatricIntensity(end); % Normalize by last value (CA = 180°)

subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12°)'])

plot(filteredAngles, ratiomatricIntensityNorm, 'DisplayName', '0°/30° ratio')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/30° Intensity (a.u.)')
legend
grid on


%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2°, ±12°, and ±16°, 
% normalized to the first point in the selected range.
% for 0 degree emission only

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on the selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to Sum Intensity Over a Given Range
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for 0° at:
% - **Single-point (no summation)**
% - **Summed over ±2°**
% - **Summed over ±12°**
% - **Summed over ±16°**
intensityZero_single = histStorage(:, 270);  % Single intensity at 270°
intensityZero_2deg = sumIntensityRange(histStorage, 270, 2);  % 270 ±2°
intensityZero_12deg = sumIntensityRange(histStorage, 270, 12); % 270 ±12°
intensityZero_16deg = sumIntensityRange(histStorage, 270, 16); % 270 ±16°

% Filter intensity values based on the selected contact angle range
intensityZero_single = intensityZero_single(angleRange);
intensityZero_2deg = intensityZero_2deg(angleRange);
intensityZero_12deg = intensityZero_12deg(angleRange);
intensityZero_16deg = intensityZero_16deg(angleRange);

% Normalize all intensities to the first value in the selected range
normalize = @(data) data / data(1);  % Function for normalization

intensityZero_singleNorm = normalize(intensityZero_single);
intensityZero_2degNorm = normalize(intensityZero_2deg);
intensityZero_12degNorm = normalize(intensityZero_12deg);
intensityZero_16degNorm = normalize(intensityZero_16deg);

% Plot All Four Intensities in One Graph
figure; hold on;
title(['Normalized Bottom Readout Intensity at 0° (CA ' num2str(minAngle) '-' num2str(maxAngle) ', First-Point Normalized)'])

plot(filteredAngles, intensityZero_singleNorm, '-d', 'LineWidth', 1.5, 'DisplayName', '0° Single-Point Intensity')
plot(filteredAngles, intensityZero_2degNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0° ±2° Summed Intensity')
plot(filteredAngles, intensityZero_12degNorm, '-^', 'LineWidth', 1.5, 'DisplayName', '0° ±12° Summed Intensity')
plot(filteredAngles, intensityZero_16degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '0° ±16° Summed Intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2°, ±12°, and ±16°, 
% normalized to the first point in the selected range.
% for emission angle of 30 degrees

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on the selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to Sum Intensity Over a Given Range
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for **30°** at:
% - **Single-point (no summation)**
% - **Summed over ±2°**
% - **Summed over ±12°**
% - **Summed over ±16°**
intensity30_single = histStorage(:, 300);  % Single intensity at 300° (corresponds to 30°)
intensity30_2deg = sumIntensityRange(histStorage, 300, 2);  % 300 ±2°
intensity30_12deg = sumIntensityRange(histStorage, 300, 12); % 300 ±12°
intensity30_16deg = sumIntensityRange(histStorage, 300, 16); % 300 ±16°

% Filter intensity values based on the selected contact angle range
intensity30_single = intensity30_single(angleRange);
intensity30_2deg = intensity30_2deg(angleRange);
intensity30_12deg = intensity30_12deg(angleRange);
intensity30_16deg = intensity30_16deg(angleRange);

% Normalize all intensities to the first value in the selected range
normalize = @(data) data / data(1);  % Function for normalization

intensity30_singleNorm = normalize(intensity30_single);
intensity30_2degNorm = normalize(intensity30_2deg);
intensity30_12degNorm = normalize(intensity30_12deg);
intensity30_16degNorm = normalize(intensity30_16deg);

% Plot All Four Intensities in One Graph for **30°**
figure; hold on;
title(['Normalized Bottom Readout Intensity at 30° (CA ' num2str(minAngle) '-' num2str(maxAngle) ', First-Point Normalized)'])

plot(filteredAngles, intensity30_singleNorm, '-d', 'LineWidth', 1.5, 'DisplayName', '30° Single-Point Intensity')
plot(filteredAngles, intensity30_2degNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '30° ±2° Summed Intensity')
plot(filteredAngles, intensity30_12degNorm, '-^', 'LineWidth', 1.5, 'DisplayName', '30° ±12° Summed Intensity')
plot(filteredAngles, intensity30_16degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '30° ±16° Summed Intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Ratiometric Bottom Readout (0°/30°) from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Ratiometric Intensities (0°/30°) for ±2°, ±12°, and ±16°,
% normalized to the first point in the selected range.
% for 0/30 ratiometric 

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on the selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to Sum Intensity Over a Given Range
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for **0°** and **30°** at:
% - **Single-point (no summation)**
% - **Summed over ±2°**
% - **Summed over ±12°**
% - **Summed over ±16°**
intensity0_single = histStorage(:, 270);  % Single intensity at 0° (270° in histStorage)
intensity0_2deg = sumIntensityRange(histStorage, 270, 2);  % 0° ±2°
intensity0_12deg = sumIntensityRange(histStorage, 270, 12); % 0° ±12°
intensity0_16deg = sumIntensityRange(histStorage, 270, 16); % 0° ±16°

intensity30_single = histStorage(:, 300);  % Single intensity at 30° (300° in histStorage)
intensity30_2deg = sumIntensityRange(histStorage, 300, 2);  % 30° ±2°
intensity30_12deg = sumIntensityRange(histStorage, 300, 12); % 30° ±12°
intensity30_16deg = sumIntensityRange(histStorage, 300, 16); % 30° ±16°

% Filter intensity values based on the selected contact angle range
intensity0_single = intensity0_single(angleRange);
intensity0_2deg = intensity0_2deg(angleRange);
intensity0_12deg = intensity0_12deg(angleRange);
intensity0_16deg = intensity0_16deg(angleRange);

intensity30_single = intensity30_single(angleRange);
intensity30_2deg = intensity30_2deg(angleRange);
intensity30_12deg = intensity30_12deg(angleRange);
intensity30_16deg = intensity30_16deg(angleRange);

% Compute Ratiometric Intensities (0°/30°)
ratio_0_30_single = intensity0_single ./ intensity30_single;
ratio_0_30_2deg = intensity0_2deg ./ intensity30_2deg;
ratio_0_30_12deg = intensity0_12deg ./ intensity30_12deg;
ratio_0_30_16deg = intensity0_16deg ./ intensity30_16deg;

% Normalize all ratios to the first value in the selected range
normalize = @(data) data / data(1);  % Function for normalization

ratio_0_30_singleNorm = normalize(ratio_0_30_single);
ratio_0_30_2degNorm = normalize(ratio_0_30_2deg);
ratio_0_30_12degNorm = normalize(ratio_0_30_12deg);
ratio_0_30_16degNorm = normalize(ratio_0_30_16deg);

% Plot Normalized Ratiometric Intensities (0°/30°) for Different Summation Ranges
figure; hold on;
title(['Normalized Ratiometric Bottom Readout (0°/30°) (CA ' num2str(minAngle) '-' num2str(maxAngle) ', First-Point Normalized)'])

plot(filteredAngles, ratio_0_30_singleNorm, '-d', 'LineWidth', 1.5, 'DisplayName', '0°/30° Single-Point')
plot(filteredAngles, ratio_0_30_2degNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°/30° ±2° Summed')
plot(filteredAngles, ratio_0_30_12degNorm, '-^', 'LineWidth', 1.5, 'DisplayName', '0°/30° ±12° Summed')
plot(filteredAngles, ratio_0_30_16degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '0°/30° ±16° Summed')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Intensity (0°/30°) (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±12°, normalized to CA = 90°
% normalize to CA=90

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around each index
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute Intensities for Different Contact Angles (Summed over ±12°)
intensityZero = sumIntensityRange(histStorage, 270); % 0° ±12°
intensity15 = sumIntensityRange(histStorage, 285);  % 15° ±12°
intensity30 = sumIntensityRange(histStorage, 300);  % 30° ±12°
intensity45 = sumIntensityRange(histStorage, 315);  % 45° ±12°
intensity60 = sumIntensityRange(histStorage, 330);  % 60° ±12°
intensity75 = sumIntensityRange(histStorage, 345);  % 75° ±12°
intensity90 = sumIntensityRange(histStorage, 360);  % 90° ±12°

% Filter intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity15 = intensity15(angleRange);
intensity30 = intensity30(angleRange);
intensity45 = intensity45(angleRange);
intensity60 = intensity60(angleRange);
intensity75 = intensity75(angleRange);
intensity90 = intensity90(angleRange);

% Normalize intensities by their first values in the selected range (CA = 90°)
normalize = @(data) data / data(1);
intensityZeroNorm = normalize(intensityZero);
intensity15Norm = normalize(intensity15);
intensity30Norm = normalize(intensity30);
intensity45Norm = normalize(intensity45);
intensity60Norm = normalize(intensity60);
intensity75Norm = normalize(intensity75);
intensity90Norm = normalize(intensity90);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, intensityZeroNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity15Norm, '-s', 'LineWidth', 1.5, 'DisplayName', '15°-intensity')
plot(filteredAngles, intensity30Norm, '-^', 'LineWidth', 1.5, 'DisplayName', '30°-intensity')
plot(filteredAngles, intensity45Norm, '-d', 'LineWidth', 1.5, 'DisplayName', '45°-intensity')
plot(filteredAngles, intensity60Norm, '-p', 'LineWidth', 1.5, 'DisplayName', '60°-intensity')
plot(filteredAngles, intensity75Norm, '-h', 'LineWidth', 1.5, 'DisplayName', '75°-intensity')
plot(filteredAngles, intensity90Norm, '-x', 'LineWidth', 1.5, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

% Compute and plot normalized ratiometric intensity (0°/30°)
ratiometricIntensity = intensityZero ./ intensity30;
ratiometricIntensityNorm = normalize(ratiometricIntensity);

subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (0°/30°) (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, ratiometricIntensityNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°/30° ratio')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/30° Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±12°, normalized to CA = 180°
% normalize to CA=180

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around each index
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute Intensities for Different Contact Angles (Summed over ±12°)
intensityZero = sumIntensityRange(histStorage, 270); % 0° ±12°
intensity15 = sumIntensityRange(histStorage, 285);  % 15° ±12°
intensity30 = sumIntensityRange(histStorage, 300);  % 30° ±12°
intensity45 = sumIntensityRange(histStorage, 315);  % 45° ±12°
intensity60 = sumIntensityRange(histStorage, 330);  % 60° ±12°
intensity75 = sumIntensityRange(histStorage, 345);  % 75° ±12°
intensity90 = sumIntensityRange(histStorage, 360);  % 90° ±12°

% Filter intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity15 = intensity15(angleRange);
intensity30 = intensity30(angleRange);
intensity45 = intensity45(angleRange);
intensity60 = intensity60(angleRange);
intensity75 = intensity75(angleRange);
intensity90 = intensity90(angleRange);

% Normalize intensities by their last value in the selected range (CA = 180°)
normalize = @(data) data / data(end);
intensityZeroNorm = normalize(intensityZero);
intensity15Norm = normalize(intensity15);
intensity30Norm = normalize(intensity30);
intensity45Norm = normalize(intensity45);
intensity60Norm = normalize(intensity60);
intensity75Norm = normalize(intensity75);
intensity90Norm = normalize(intensity90);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, intensityZeroNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity15Norm, '-s', 'LineWidth', 1.5, 'DisplayName', '15°-intensity')
plot(filteredAngles, intensity30Norm, '-^', 'LineWidth', 1.5, 'DisplayName', '30°-intensity')
plot(filteredAngles, intensity45Norm, '-d', 'LineWidth', 1.5, 'DisplayName', '45°-intensity')
plot(filteredAngles, intensity60Norm, '-p', 'LineWidth', 1.5, 'DisplayName', '60°-intensity')
plot(filteredAngles, intensity75Norm, '-h', 'LineWidth', 1.5, 'DisplayName', '75°-intensity')
plot(filteredAngles, intensity90Norm, '-x', 'LineWidth', 1.5, 'DisplayName', '90°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

% Compute and plot normalized ratiometric intensity (0°/30°)
ratiometricIntensity = intensityZero ./ intensity30;
ratiometricIntensityNorm = normalize(ratiometricIntensity);

subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (0°/30°) (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, ratiometricIntensityNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°/30° ratio')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric 0°/30° Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±12°, normalized to CA = 180°

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around each index
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute Intensities for Different Contact Angles (Summed over ±12°)
intensityZero = sumIntensityRange(histStorage, 270); % 0° ±12°
intensity30 = sumIntensityRange(histStorage, 300);  % 30° ±12°

% Filter intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity30 = intensity30(angleRange);

% Normalize intensities by their last value in the selected range (CA = 180°)
normalize = @(data) data / data(end);
intensityZeroNorm = normalize(intensityZero);
intensity30Norm = normalize(intensity30);

% Compute and normalize ratiometric intensities
ratiometric_ZeroOver30 = intensityZero ./ intensity30;
ratiometric_30OverZero = intensity30 ./ intensityZero;
ratiometric_ZeroOver30Norm = normalize(ratiometric_ZeroOver30);
ratiometric_30OverZeroNorm = normalize(ratiometric_30OverZero);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, intensityZeroNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity30Norm, '-^', 'LineWidth', 1.5, 'DisplayName', '30°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

% Plot 0°/30° and 30°/0° normalized ratios on the same graph
subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (0°/30° vs 30°/0°) (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, ratiometric_ZeroOver30Norm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°/30° ratio')
plot(filteredAngles, ratiometric_30OverZeroNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '30°/0° ratio')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles
% Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±12°, normalized to CA = 90°

% Define Contact Angle Range
minAngle = 90;  % Minimum contact angle
maxAngle = 180; % Maximum contact angle

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around each index
sumIntensityRange = @(data, idx) sum(data(:, max(1, idx-12) : min(size(data, 2), idx+12)), 2);

% Compute Intensities for Different Contact Angles (Summed over ±12°)
intensityZero = sumIntensityRange(histStorage, 270); % 0° ±12°
intensity30 = sumIntensityRange(histStorage, 300);  % 30° ±12°

% Filter intensity values based on the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity30 = intensity30(angleRange);

% Normalize intensities by their first value in the selected range (CA = 90°)
normalize = @(data) data / data(1);
intensityZeroNorm = normalize(intensityZero);
intensity30Norm = normalize(intensity30);

% Compute and normalize ratiometric intensities
ratiometric_ZeroOver30 = intensityZero ./ intensity30;
ratiometric_30OverZero = intensity30 ./ intensityZero;
ratiometric_ZeroOver30Norm = normalize(ratiometric_ZeroOver30);
ratiometric_30OverZeroNorm = normalize(ratiometric_30OverZero);

% Create figure and plot normalized intensities for the selected contact angle range
figure;
subplot(2,1,1); hold on;
title(['Normalized Bottom Readout Intensity (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, intensityZeroNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°-intensity')
plot(filteredAngles, intensity30Norm, '-^', 'LineWidth', 1.5, 'DisplayName', '30°-intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

% Plot 0°/30° and 30°/0° normalized ratios on the same graph
subplot(2,1,2); hold on;
title(['Normalized Ratiometric Bottom Readout (0°/30° vs 30°/0°) (CA ' num2str(minAngle) '-' num2str(maxAngle) ', ±12° Summed)'])

plot(filteredAngles, ratiometric_ZeroOver30Norm, '-o', 'LineWidth', 1.5, 'DisplayName', '0°/30° ratio')
plot(filteredAngles, ratiometric_30OverZeroNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '30°/0° ratio')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Normalized Bottom Readout from CA90-180 vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2°, ±12°, ±16°
% Normalized to CA = 90°

% Define the range of contact angles to plot
minAngle = 90;  
maxAngle = 180; 

% Automatically filter contact angles based on selected range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over a given range
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute intensities at 90° (360° in `histStorage`) for different summation ranges
intensity90_2deg = sumIntensityRange(histStorage, 360, 2);  % 90° ±2°
intensity90_12deg = sumIntensityRange(histStorage, 360, 12); % 90° ±12°
intensity90_16deg = sumIntensityRange(histStorage, 360, 16); % 90° ±16°

% Filter intensity values based on the selected contact angle range
intensity90_2deg = intensity90_2deg(angleRange);
intensity90_12deg = intensity90_12deg(angleRange);
intensity90_16deg = intensity90_16deg(angleRange);

% Normalize intensities by their first values in the selected range (CA = 90°)
intensity90_2degNorm = intensity90_2deg / intensity90_2deg(1);
intensity90_12degNorm = intensity90_12deg / intensity90_12deg(1);
intensity90_16degNorm = intensity90_16deg / intensity90_16deg(1);

% Plot Normalized Intensities for 90° Overlaid in One Graph
figure; hold on;
title(['Normalized Bottom Readout Intensity at 90° (CA ' num2str(minAngle) '-' num2str(maxAngle) ')'])

plot(filteredAngles, intensity90_2degNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '90° ±2° Intensity')
plot(filteredAngles, intensity90_12degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '90° ±12° Intensity')
plot(filteredAngles, intensity90_16degNorm, '-d', 'LineWidth', 1.5, 'DisplayName', '90° ±16° Intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025

% Initialize intensity arrays
intensityZero = [];
intensity15 = [];
intensity30 = [];
intensity45 = [];
intensity60 = [];
intensity75 = [];
intensity90 = [];
ratiomatricIntensity = [];

% Function to sum intensity over ±12 degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute intensities (Original Data)
for i = 1:length(contactAngles)
    intensityZero(i) = mean(histStorage(i, 269:271)); % 270
    intensity15(i) = mean(histStorage(i, 284:286)); % 285
    intensity30(i) = mean(histStorage(i, 299:301)); % 300
    intensity45(i) = mean(histStorage(i, 314:316)); % 315
    intensity60(i) = mean(histStorage(i, 329:331)); % 330
    intensity75(i) = mean(histStorage(i, 344:346)); % 345
    intensity90(i) = mean(histStorage(i, 359:360)); % 360
end

% Compute Intensities Summed Over ±12° Range
intensityZero_12deg = sumIntensityRange(histStorage, 270, 12);
intensity15_12deg = sumIntensityRange(histStorage, 285, 12);
intensity30_12deg = sumIntensityRange(histStorage, 300, 12);
intensity45_12deg = sumIntensityRange(histStorage, 315, 12);
intensity60_12deg = sumIntensityRange(histStorage, 330, 12);
intensity75_12deg = sumIntensityRange(histStorage, 345, 12);
intensity90_12deg = sumIntensityRange(histStorage, 360, 12);

% Plot Original Intensities
figure;
subplot(2,1,1); hold on;
title('Degrees of Bottom Readout Intensity')

plot(contactAngles, intensityZero, '-o', 'DisplayName', '0° Intensity (Original)')
plot(contactAngles, intensity15, '-o', 'DisplayName', '15° Intensity (Original)')
plot(contactAngles, intensity30, '-o', 'DisplayName', '30° Intensity (Original)')
plot(contactAngles, intensity45, '-o', 'DisplayName', '45° Intensity (Original)')
plot(contactAngles, intensity60, '-o', 'DisplayName', '60° Intensity (Original)')
plot(contactAngles, intensity75, '-o', 'DisplayName', '75° Intensity (Original)')
plot(contactAngles, intensity90, '-o', 'DisplayName', '90° Intensity (Original)')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)')
legend('Location', 'best')
grid on

% Plot Intensities Summed Over ±12° for Comparison
subplot(2,1,2); hold on;
title('Summed Intensities (±12°) for Bottom Readout')

plot(contactAngles, intensityZero_12deg, '-s', 'DisplayName', '0° ±12° Intensity')
plot(contactAngles, intensity15_12deg, '-s', 'DisplayName', '15° ±12° Intensity')
plot(contactAngles, intensity30_12deg, '-s', 'DisplayName', '30° ±12° Intensity')
plot(contactAngles, intensity45_12deg, '-s', 'DisplayName', '45° ±12° Intensity')
plot(contactAngles, intensity60_12deg, '-s', 'DisplayName', '60° ±12° Intensity')
plot(contactAngles, intensity75_12deg, '-s', 'DisplayName', '75° ±12° Intensity')
plot(contactAngles, intensity90_12deg, '-s', 'DisplayName', '90° ±12° Intensity')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Summed Intensity (a.u.)')
legend('Location', 'best')
grid on

% Compute Ratiometric Intensity for Original Data
ratiomatricIntensity = intensityZero ./ intensity45;

% Compute Ratiometric Intensity for ±12° Summed Data
ratiomatricIntensity_12deg = intensityZero_12deg ./ intensity45_12deg;

% Plot Ratiometric Comparison
figure; hold on;
title('Ratiometric Bottom Readout (0° / 45°)')

plot(contactAngles, ratiomatricIntensity / ratiomatricIntensity(1), '-o', 'LineWidth', 1.5, 'DisplayName', '0° / 45° Ratio (Original)')
plot(contactAngles, ratiomatricIntensity_12deg / ratiomatricIntensity_12deg(1), '-s', 'LineWidth', 1.5, 'DisplayName', '0° / 45° Ratio (±12°)')

xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Intensity (a.u.)')
legend('Location', 'best')
grid on
hold off;

%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±12°, Normalized to CA=180°
% Contact Angle Range Limited to 90° - 180°

% Define Contact Angle Range
minAngle = 90;  
maxAngle = 180; 

% Filter contact angles within the specified range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Function to sum intensity over ±12 degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for 0° and 30° Readout (Original Data)
intensityZero = mean(histStorage(:, 269:271), 2); % 270° (0° readout)
intensity30 = mean(histStorage(:, 299:301), 2); % 300° (30° readout)

% Compute Intensities Summed Over ±12° Range
intensityZero_12deg = sumIntensityRange(histStorage, 270, 12);
intensity30_12deg = sumIntensityRange(histStorage, 300, 12);

% Filter intensity values within the selected contact angle range
intensityZero = intensityZero(angleRange);
intensity30 = intensity30(angleRange);
intensityZero_12deg = intensityZero_12deg(angleRange);
intensity30_12deg = intensity30_12deg(angleRange);

% Normalize intensities to last value (CA = 180°)
intensityZeroNorm = intensityZero / intensityZero(end);
intensity30Norm = intensity30 / intensity30(end);
intensityZero_12degNorm = intensityZero_12deg / intensityZero_12deg(end);
intensity30_12degNorm = intensity30_12deg / intensity30_12deg(end);

% Compute Ratiometric Intensity for 0°/30°
ratiometricIntensity = intensityZero ./ intensity30;
ratiometricIntensity_12deg = intensityZero_12deg ./ intensity30_12deg;

% Normalize ratiometric intensity to CA = 180°
ratiometricIntensityNorm = ratiometricIntensity / ratiometricIntensity(end);
ratiometricIntensity_12degNorm = ratiometricIntensity_12deg / ratiometricIntensity_12deg(end);

% Create figure with three subplots
figure;

% Subplot 1: Original Intensities
subplot(3,1,1); hold on;
title('Normalized Bottom Readout Intensity (Original, 90°-180°)')
plot(filteredAngles, intensityZeroNorm, '-o', 'LineWidth', 1.5, 'DisplayName', '0° Intensity')
plot(filteredAngles, intensity30Norm, '-o', 'LineWidth', 1.5, 'DisplayName', '30° Intensity')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: ±12° Summed Intensities
subplot(3,1,2); hold on;
title('Normalized Bottom Readout Intensity (Summed ±12°, 90°-180°)')
plot(filteredAngles, intensityZero_12degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '0° ±12° Intensity')
plot(filteredAngles, intensity30_12degNorm, '-s', 'LineWidth', 1.5, 'DisplayName', '30° ±12° Intensity')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 3: Ratiometric Intensities (0°/30°)
subplot(3,1,3); hold on;
title('Normalized Ratiometric 0°/30° Readout (90°-180°)')
plot(filteredAngles, ratiometricIntensityNorm, '--', 'LineWidth', 2, 'DisplayName', '0°/30° Ratio (Original)')
plot(filteredAngles, ratiometricIntensity_12degNorm, '--', 'LineWidth', 2, 'DisplayName', '0°/30° Ratio (±12°)')
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;


%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2° and ±12°, Normalized to CA=180°
% Contact Angle Range Limited to 90° - 180°

% Define Contact Angle Range
minAngle = 90;  
maxAngle = 180; 

% Filter contact angles within the specified range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Define Angles for Computation
angleIndices = [270, 285, 300, 315, 330, 345, 360]; % Corresponding indices for 0°, 15°, 30°, etc.
angleLabels = {'0°', '15°', '30°', '45°', '60°', '75°', '90°'};

% Function to sum intensity over ±X degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for All Angles (±2°)
intensities_2deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_2deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 2);
end

% Compute Intensities for All Angles (±12°)
intensities_12deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_12deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 12);
end

% Filter intensity values within the selected contact angle range
intensities_2deg = intensities_2deg(angleRange, :);
intensities_12deg = intensities_12deg(angleRange, :);

% Normalize intensities to CA = 180°
intensities_2degNorm = intensities_2deg ./ intensities_2deg(end, :);
intensities_12degNorm = intensities_12deg ./ intensities_12deg(end, :);

% Compute Ratiometric Intensities for 0° over All Angles (±2° and ±12°)
ratiometric_2deg = intensities_2deg(:, 1) ./ intensities_2deg(:, 2:end);
ratiometric_12deg = intensities_12deg(:, 1) ./ intensities_12deg(:, 2:end);

% Normalize ratiometric intensity to CA = 180°
ratiometric_2degNorm = ratiometric_2deg ./ ratiometric_2deg(end, :);
ratiometric_12degNorm = ratiometric_12deg ./ ratiometric_12deg(end, :);

% Create figure with four subplots
figure;

% Subplot 1: Normalized Intensities for All Angles (±2°)
subplot(4,1,1); hold on;
title('Normalized Bottom Readout Intensity (Summed ±2°, 90°-180°)')
for i = 1:size(intensities_2degNorm,2)
    plot(filteredAngles, intensities_2degNorm(:, i), '-o', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Normalized Intensities for All Angles (±12°)
subplot(4,1,2); hold on;
title('Normalized Bottom Readout Intensity (Summed ±12°, 90°-180°)')
for i = 1:size(intensities_12degNorm,2)
    plot(filteredAngles, intensities_12degNorm(:, i), '-s', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 3: Ratiometric Intensities for ±2°
subplot(4,1,3); hold on;
title('Normalized Ratiometric 0° Readout over All Angles (±2°, 90°-180°)')
for i = 1:size(ratiometric_2degNorm,2)
    plot(filteredAngles, ratiometric_2degNorm(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±2°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

% Subplot 4: Ratiometric Intensities for ±12°
subplot(4,1,4); hold on;
title('Normalized Ratiometric 0° Readout over All Angles (±12°, 90°-180°)')
for i = 1:size(ratiometric_12degNorm,2)
    plot(filteredAngles, ratiometric_12degNorm(:, i), ':', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±12°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2° and ±12°, Normalized to CA=180°
% Contact Angle Range Limited to 90° - 180°

% Define Contact Angle Range
minAngle = 90;  
maxAngle = 180; 

% Filter contact angles within the specified range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Define Angles for Computation
angleIndices = [270, 285, 300, 315, 330, 345, 360]; % Corresponding indices for 0°, 15°, 30°, etc.
angleLabels = {'0°', '15°', '30°', '45°', '60°', '75°', '90°'};

% Function to sum intensity over ±X degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for All Angles (±2°)
intensities_2deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_2deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 2);
end

% Compute Intensities for All Angles (±12°)
intensities_12deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_12deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 12);
end

% Filter intensity values within the selected contact angle range
intensities_2deg = intensities_2deg(angleRange, :);
intensities_12deg = intensities_12deg(angleRange, :);

% Normalize intensities to CA = 180°
intensities_2degNorm = intensities_2deg ./ intensities_2deg(end, :);
intensities_12degNorm = intensities_12deg ./ intensities_12deg(end, :);

% Compute Ratiometric Intensities for 0° over all other angles (±2° and ±12°)
ratiometric_2deg = intensities_2deg(:, 1) ./ intensities_2deg(:, 2:end);
ratiometric_12deg = intensities_12deg(:, 1) ./ intensities_12deg(:, 2:end);

% Normalize ratiometric intensity to CA = 180°
ratiometric_2degNorm = ratiometric_2deg ./ ratiometric_2deg(end, :);
ratiometric_12degNorm = ratiometric_12deg ./ ratiometric_12deg(end, :);

% **Figure 1: All Data for ±2°**
figure;

% Subplot 1: Normalized Intensities for All Angles (±2°)
subplot(2,1,1); hold on;
title('Normalized Bottom Readout Intensity (Summed ±2°, 90°-180°)')
for i = 1:size(intensities_2degNorm,2)
    plot(filteredAngles, intensities_2degNorm(:, i), '-o', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Ratiometric Intensities 0° over All Angles (±2°)
subplot(2,1,2); hold on;
title('Normalized Ratiometric 0° Readout over All Angles (±2°)')
for i = 1:size(ratiometric_2degNorm,2)
    plot(filteredAngles, ratiometric_2degNorm(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±2°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

% **Figure 2: All Data for ±12°**
figure;

% Subplot 1: Normalized Intensities for All Angles (±12°)
subplot(2,1,1); hold on;
title('Normalized Bottom Readout Intensity (Summed ±12°, 90°-180°)')
for i = 1:size(intensities_12degNorm,2)
    plot(filteredAngles, intensities_12degNorm(:, i), '-s', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Ratiometric Intensities 0° over All Angles (±12°)
subplot(2,1,2); hold on;
title('Normalized Ratiometric 0° Readout over All Angles (±12°)')
for i = 1:size(ratiometric_12degNorm,2)
    plot(filteredAngles, ratiometric_12degNorm(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±12°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2°, ±5°, and ±12°, No Normalization
% Contact Angle Range Limited to 90° - 180°

% Define Contact Angle Range
minAngle = 90;  
maxAngle = 180; 

% Filter contact angles within the specified range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Define Angles for Computation
angleIndices = [270, 285, 300, 315, 330, 345, 360]; % Corresponding indices for 0°, 15°, 30°, etc.
angleLabels = {'0°', '15°', '30°', '45°', '60°', '75°', '90°'};

% Function to sum intensity over ±X degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for All Angles
intensities_2deg = zeros(length(contactAngles), length(angleIndices));
intensities_5deg = zeros(length(contactAngles), length(angleIndices));
intensities_12deg = zeros(length(contactAngles), length(angleIndices));

for i = 1:length(angleIndices)
    intensities_2deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 2);
    intensities_5deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 5);
    intensities_12deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 12);
end

% Filter intensity values within the selected contact angle range
intensities_2deg = intensities_2deg(angleRange, :);
intensities_5deg = intensities_5deg(angleRange, :);
intensities_12deg = intensities_12deg(angleRange, :);

% Compute Ratiometric Intensities for 0° over all other angles
ratiometric_2deg = intensities_2deg(:, 1) ./ intensities_2deg(:, 2:end);
ratiometric_5deg = intensities_5deg(:, 1) ./ intensities_5deg(:, 2:end);
ratiometric_12deg = intensities_12deg(:, 1) ./ intensities_12deg(:, 2:end);

% **Figure 1: Original Values & Ratiometric Intensities for ±2°**
figure;

% Subplot 1: Original Intensities for All Angles (±2°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±2°, 90°-180°)')
for i = 1:size(intensities_2deg,2)
    plot(filteredAngles, intensities_2deg(:, i), '-o', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Ratiometric Intensities 0° over All Angles (±2°)
subplot(2,1,2); hold on;
title('Ratiometric 0° Readout over All Angles (±2°)')
for i = 1:size(ratiometric_2deg,2)
    plot(filteredAngles, ratiometric_2deg(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±2°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

% **Figure 2: Original Values & Ratiometric Intensities for ±5°**
figure;

% Subplot 1: Original Intensities for All Angles (±5°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±5°, 90°-180°)')
for i = 1:size(intensities_5deg,2)
    plot(filteredAngles, intensities_5deg(:, i), '-d', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Ratiometric Intensities 0° over All Angles (±5°)
subplot(2,1,2); hold on;
title('Ratiometric 0° Readout over All Angles (±5°)')
for i = 1:size(ratiometric_5deg,2)
    plot(filteredAngles, ratiometric_5deg(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±5°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

% **Figure 3: Original Values & Ratiometric Intensities for ±12°**
figure;

% Subplot 1: Original Intensities for All Angles (±12°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±12°, 90°-180°)')
for i = 1:size(intensities_12deg,2)
    plot(filteredAngles, intensities_12deg(:, i), '-s', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)')
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Ratiometric Intensities 0° over All Angles (±12°)
subplot(2,1,2); hold on;
title('Ratiometric 0° Readout over All Angles (±12°)')
for i = 1:size(ratiometric_12deg,2)
    plot(filteredAngles, ratiometric_12deg(:, i), '--', 'LineWidth', 2, 'DisplayName', ['0°/' angleLabels{i+1} ' Ratio (±12°)']);
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Ratiometric Value')
legend('Location', 'best')
grid on;
hold off;

%% Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025
% NA = 0.22 considered, Summed Intensities for CA ±2°, ±5°, and ±12°, No Normalization
% Contact Angle Range Limited to 90° - 180°

% Define Contact Angle Range
minAngle = 90;  
maxAngle = 180; 

% Filter contact angles within the specified range
angleRange = (contactAngles >= minAngle) & (contactAngles <= maxAngle);
filteredAngles = contactAngles(angleRange);

% Define Angles for Computation
angleIndices = [270, 285, 300, 315, 330, 345, 360]; % Corresponding indices for 0°, 15°, 30°, etc.
angleLabels = {'0°', '15°', '30°', '45°', '60°', '75°', '90°'};

% Function to sum intensity over ±X degrees around a given index
sumIntensityRange = @(data, idx, range) sum(data(:, max(1, idx-range) : min(size(data, 2), idx+range)), 2);

% Compute Intensities for All Angles (±2°)
intensities_2deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_2deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 2);
end

% Compute Intensities for All Angles (±5°)
intensities_5deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_5deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 5);
end

% Compute Intensities for All Angles (±12°)
intensities_12deg = zeros(length(contactAngles), length(angleIndices));
for i = 1:length(angleIndices)
    intensities_12deg(:, i) = sumIntensityRange(histStorage, angleIndices(i), 12);
end

% Filter intensity values within the selected contact angle range
intensities_2deg = intensities_2deg(angleRange, :);
intensities_5deg = intensities_5deg(angleRange, :);
intensities_12deg = intensities_12deg(angleRange, :);

% **Figure 1: Original Values for ±2°**
figure;

% Subplot 1: Original Intensities for All Angles (±2°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±2°, 90°-180°)')
for i = 1:size(intensities_2deg,2)
    plot(filteredAngles, intensities_2deg(:, i), '-o', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Duplicate of the first plot for comparison
subplot(2,1,2); hold on;
title('Original Bottom Readout Intensity (Summed ±2°, 90°-180°) (All Angles)')
for i = 1:size(intensities_2deg,2)
    plot(filteredAngles, intensities_2deg(:, i), '-o', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;

% **Figure 2: Original Values for ±5°**
figure;

% Subplot 1: Original Intensities for All Angles (±5°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±5°, 90°-180°)')
for i = 1:size(intensities_5deg,2)
    plot(filteredAngles, intensities_5deg(:, i), '-d', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Duplicate of the first plot for comparison
subplot(2,1,2); hold on;
title('Original Bottom Readout Intensity (Summed ±5°, 90°-180°) (All Angles)')
for i = 1:size(intensities_5deg,2)
    plot(filteredAngles, intensities_5deg(:, i), '-d', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;

% **Figure 3: Original Values for ±12°**
figure;

% Subplot 1: Original Intensities for All Angles (±12°)
subplot(2,1,1); hold on;
title('Original Bottom Readout Intensity (Summed ±12°, 90°-180°)')
for i = 1:size(intensities_12deg,2)
    plot(filteredAngles, intensities_12deg(:, i), '-s', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;

% Subplot 2: Duplicate of the first plot for comparison
subplot(2,1,2); hold on;
title('Original Bottom Readout Intensity (Summed ±12°, 90°-180°) (All Angles)')
for i = 1:size(intensities_12deg,2)
    plot(filteredAngles, intensities_12deg(:, i), '-s', 'LineWidth', 1.5, 'DisplayName', angleLabels{i});
end
xlabel('Contact Angles (\Theta\circ)')
ylabel('Intensity (a.u.)') % No normalization
legend('Location', 'best')
grid on;
hold off;


%% Normalized Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025, 0deg only

% Compute and normalize only the 0-degree intensity
intensityZero = arrayfun(@(i) mean(histStorage(i,269:271)), 1:length(contactAngles));
%intensityZeroNorm = intensityZero / intensityZero(1); normalize with the
%first value
intensityZeroNorm = intensityZero / intensityZero(end); % normalize wtih the last value

% Plot normalized 0-degree intensity
figure; hold on;
title('Normalized 0° Bottom Readout Intensity')

plot(contactAngles, intensityZeroNorm, 'DisplayName', '0°-intensity', 'LineWidth', 2)
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend
grid on

%% Normalized Bottom Readout vs Contact Angles - Written by Katrina, Copyright 2025

% Extract and normalize the 0-degree intensity directly from histStorage
intensityZero = histStorage(:, 270); % Extract the 270th column directly
intensityZeroNorm = intensityZero / intensityZero(1); % Normalize by first value

% Plot normalized 0-degree intensity
figure; hold on;
title('Normalized 0° Bottom Readout Intensity')

plot(contactAngles, intensityZeroNorm, 'DisplayName', '0°-intensity', 'LineWidth', 2)
xlabel('Contact Angles (\Theta\circ)')
ylabel('Normalized Intensity (a.u.)')
legend
grid on

%% Scanning intensity profile of a droplet from bottom 0 to 45(BinAngles 270 to 360)

variousAnglesIntensity = [];
figure; hold on;
titlestr = ['Scanning intensity profile of a droplet from bottom 0\circ to 70\circ' '-VR=' num2str(vr)];
title(titlestr)
xlabel('Detection angles(\circ)')
ylabel('Intensity(a.u.)')
%selectDroplets = [51  63  75  87 99 111 123];
selectDroplets = [90 100 110 120 130 140 150 160 170 180];
%selectDroplets = [20 32 40 52 60 72];

normalization = (histStorage/max(histStorage));
%normalization = (histStorage2/max(histStorage2)); i dont remember what is
%this 2 now, its 2024 now, by Katrina; ah i know now, histStorage is for
%incident light shining intensity, histStorage2 is for fluorescent
%intensity

%noBaselineHistStorage = histStorage %- mean(histStorage); %this is hacky and incorrect
noBaselineHistStorage = histStorage./normalization; 


for i = 1:length(contactAngles) %contactAnglesCorrected(45)
    if ismember(contactAngles(i),selectDroplets)
        for j = 270:340
        %for j = 268:356
            variousAnglesIntensity(j) = mean(noBaselineHistStorage(i,(j-4):(j+4)));
        end
        [M,maxValue(i)] = max(variousAnglesIntensity);
        %[pks,locs] = findpeaks(variousAnglesIntensity(270:340))
        %STORAGE(i,:) = variousAnglesIntensity; 
        plot(0:70,variousAnglesIntensity(270:340),'DisplayName',num2str(contactAngles(i)),'LineWidth',1.5)
        %plot(-2:86,variousAnglesIntensity(268:356),'DisplayName',num2str(contactAngles(i)),'LineWidth',1.5)
        legend
        
        plot(maxValue(i)-270,M,'o','DisplayName',num2str(contactAngles(i)))
        %plot(maxValue(i)-268,M,'or')
    end
end

% ylim([0.75 1.5])

%% Scanning intensity profile of a droplet from bottom 0 to 45(BinAngles 270 to 360) katrina copyright 2025

variousAnglesIntensity = [];
figure; hold on;
titlestr = ['Scanning intensity profile of a droplet from bottom 0\circ to 90\circ' '-VR=' num2str(vr)];
title(titlestr)
xlabel('Detection angles(\circ)')
ylabel('Intensity(a.u.)')
%selectDroplets = [51  63  75  87 99 111 123];
selectDroplets = [90 100 110 120 130 140 150 160 170 180];
%selectDroplets = [20 32 40 52 60 72];

normalization = (histStorage/max(histStorage));
%normalization = (histStorage2/max(histStorage2)); i dont remember what is
%this 2 now, its 2024 now, by Katrina; ah i know now, histStorage is for
%incident light shining intensity, histStorage2 is for fluorescent
%intensity

%noBaselineHistStorage = histStorage %- mean(histStorage); %this is hacky and incorrect
noBaselineHistStorage = histStorage./normalization; 


for i = 1:length(contactAngles) %contactAnglesCorrected(45)
    if ismember(contactAngles(i),selectDroplets)
        for j = 270:359
        %for j = 268:356
            variousAnglesIntensity(j) = mean(noBaselineHistStorage(i,(j-1):(j+1)));
        end
        [M,maxValue(i)] = max(variousAnglesIntensity);
        %[pks,locs] = findpeaks(variousAnglesIntensity(270:340))
        %STORAGE(i,:) = variousAnglesIntensity; 
        plot(0:89,variousAnglesIntensity(270:359),'DisplayName',num2str(contactAngles(i)),'LineWidth',1.5)
        %plot(-2:86,variousAnglesIntensity(268:356),'DisplayName',num2str(contactAngles(i)),'LineWidth',1.5)
        legend
        
        plot(maxValue(i)-270,M,'o','DisplayName',num2str(contactAngles(i)))
        %plot(maxValue(i)-268,M,'or')
    end
end

% ylim([0.75 1.5])

%% Scanning intensity profile of a droplet from bottom 0° to 90° (BinAngles 270 to 360)
% Katrina, Copyright 2025

figure; hold on;
titlestr = ['Normalized Scanning Intensity Profile (0° to 90°) - VR=' num2str(vr)];
title(titlestr)
xlabel('Detection angles (\circ)')
ylabel('Normalized Intensity (a.u.)')

% Selected contact angles for plotting
selectDroplets = [90 125 155 165 175];

% Normalize intensity profile
normalization = (histStorage / max(histStorage)); % Normalize to max intensity

% Remove baseline intensity
noBaselineHistStorage = histStorage ./ normalization;

% Loop through selected contact angles
for i = 1:length(contactAngles)
    if ismember(contactAngles(i), selectDroplets)
        % Extract intensity values for angles 270 to 359
        variousAnglesIntensity = zeros(1, 90); % Preallocate
        for j = 270:359
            variousAnglesIntensity(j-269) = mean(noBaselineHistStorage(i, (j-1):(j+1))); % 3-point smoothing
        end
        
        % Normalize intensity to 90° (BinAngle = 360)
        variousAnglesIntensityNorm = variousAnglesIntensity / variousAnglesIntensity(end);
        
        % Find peak intensity value
        [M, maxValueIdx] = max(variousAnglesIntensityNorm);
        
        % Plot normalized intensity
        plot(0:89, variousAnglesIntensityNorm, 'LineWidth', 1.5, 'DisplayName', ['CA = ' num2str(contactAngles(i)) '°']);
        
        % Mark peak intensity
        plot(maxValueIdx-1, M, 'o', 'MarkerSize', 6, 'DisplayName', ['Peak: ' num2str(contactAngles(i)) '°']);
    end
end

legend('Location', 'best')
grid on;
hold off;

%% surface plotting
histBinNorm = [];

%this is for after water-air or whatever external interface
figure; hold on;
titlestr = ['VR=' num2str(vr) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) '), water plane'];
title(titlestr)
[X,Y] = meshgrid(1:1:360,contactAngles);
%[X,Y] = meshgrid(1:1:360,contactAnglesCorrected);
s = pcolor(X(:,1:360),Y(:,1:360),histStorage(:,1:360));
xlim([180 360])
%ylim([90 180])
colorbar
%caxis([0 5000])
s.EdgeColor = 'none';

xlabel('Angle of Bin Rays(\circ)')
ylabel('Contact Angle of Droplet(\circ)')
zlabel('Intensity of Rays(a.u.)')
colormap(viridis)
title(titlestr)
%  rectangle('position',[260 0 20 180])
%  rectangle('position',[215 0 20 180])

%%
%this is just as light leaves the droplet enter water before air

%interest = histStorage;
%intHisStorage is before Snell's law, hisStorage is after Snell's law
interest = intHistStorage;

figure;
titlestr = ['VR=' num2str(vr) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) '),Rays Leaving Droplet, Color is Intensity'];
%[X,Y] = meshgrid(1:1:360,contactAngles);
[X,Y] = meshgrid(1:1:360,90:1:180);
s = surf(X(:,1:360),Y(:,1:360),intHistStorage(:,1:360)); %use pcolor instead of surf for 2d
s.EdgeColor = 'none';
colormap(viridis)
xlim([180 360])
xlabel('Angle Ray Is Heading After Leaving Droplet(\circ)')
ylabel('Contact Angle of Droplet(\circ)')
zlabel('Intensity of Rays(a.u.)')
title(titlestr)
colorbar

%% playaround with normalization
%this is just as light leaves the droplet enter water before air

%interest = histStorage;
%intHisStorage is before Snell's law, hisStorage is after Snell's law
interest = intHistStorage;
normalization = (histStorage2/max(histStorage2));
figure;
titlestr = ['VR=' num2str(vr) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) '),Rays Leaving Droplet, Color is Intensity'];
[X,Y] = meshgrid(1:1:360,15:165);
s = surf(X(:,190:350),Y(:,190:350),histStorage(15:165,190:350)./normalization(190:350)); %use pcolor instead of surf for 2d
%[X,Y] = meshgrid(1:1:360,contactAngles);
%s = surf(X(:,1:360),Y(:,1:360),histStorage(:,1:360)); %use pcolor instead of surf for 2d
s.EdgeColor = 'none';
colormap(viridis)
xlim([180 360])
xlabel('Angle Ray Is Heading After Leaving Droplet(\circ)')
ylabel('Contact Angle of Droplet(\circ)')
zlabel('Intensity of Rays(a.u.)')
title(titlestr)
colorbar

%%
%now histStorage is an array containing all of the points where rays hit, now
%we place it into a histogram for each ray
data = []; degIntensityZero = []; degIntensityInterest = []; degIntensityZeroSA = []; degIntensityInterestSA = [];
for k = 1:size(histStorage,2)
    data(:,k) = intHistStorage(:,k);%.*(sind(k)^2);
end

for i = 1:numDrops %reset to numdrops if you want everything
    degIntensityZero(i,:) = mean(data(i,265:275)); %+- 5deg around direct north
    degIntensityInterest(i,:) =  mean(data(i,40:50)); %+- 5deg around 45deg off axis
    thisCA(i,:) = contactAnglesCorrected(i);
end

degIntensityZeroSA = (degIntensityZero);
degIntensityInterestSA = (degIntensityInterest);
intensityRatio = (degIntensityZeroSA./(degIntensityInterestSA));

%zeroDegIntensityNorm = normalize(degIntensityZero,'range',[1 2]); %(degIntensityZero - min(degIntensityZero)) / ( max(degIntensityZero) - min(degIntensityZero) );
%degIntensityInterestNorm = normalize(degIntensityInterest,'range',[1 2]); %(degIntensityInterest - min(degIntensityInterest)) / ( max(degIntensityInterest) - min(degIntensityInterest) );

%[M,idx] = min(normThenDiv(1:20));
%M = mean(normThenDiv(80:90));
skip = 4;

h = figure;
axes1 = axes('Parent',h);
hold on;
title(['HC:',num2str(ni),', FC:',num2str(no),', cont:',num2str(environment(1)),' W/A Interface']);
xlabel('Contact Angle of Droplet');
ylabel('Relative intensity at probe (a.u.)');

%plot(contactAnglesCorrected,zeroDegIntensityNorm,'lineWidth',2)
%plot(contactAnglesCorrected,degIntensityInterestNorm,'lineWidth',2)
plot(thisCA,degIntensityZeroSA/degIntensityZeroSA(90),'-o','lineWidth',2,'Color',[colors(9,:)])
plot(thisCA,degIntensityInterestSA/degIntensityInterestSA(90),'-o','lineWidth',2,'Color',[colors(3,:)])
plot(thisCA,intensityRatio,'-o','lineWidth',2,'Color',[colors(8,:)])

%scatter(dephwatx-11,normalize(dephwaty,'range',[1 2]))

xlim([0,90])
%ylim([0,1])
%legend('0{\circ}','45{\circ}','0/45{\circ}')
%legend('85-95^{\circ}','40-50^{\circ}','0^{\circ}/45^{\circ}','45 Degree/Zero Degree','location','northeast')

%legend('HFEDEBWAT 0deg','HFEDEBWAT 45deg','C1.334-HC1.495-FC1.324 0deg','C1.334-HC1.495-FC1.324 45deg','C1.334-HC1.498-FC1.332 0deg','C1.334-HC1.498-FC1.332 45deg')
legend('0deg','45deg','ratio')
%legend('DEB:HFE','TOL:HFE/FC43','DEPH:MPFB')
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);
grid on;

%%
data = []; degIntensityZero = []; degIntensityInterest = []; degIntensityZeroSA = []; degIntensityInterestSA = [];
for k = 1:size(histStorage,2)
    data(:,k) = histStorage(:,k);%.*(sind(k)^2);
end

for i = 1:numDrops %reset to numdrops if you want everything
    degIntensityZero(i,:) = mean(data(i,86:96)); %+- 5deg around direct north
    degIntensityInterest(i,:) =  mean(data(i,41:51)); %+- 5deg around 45deg off axis
    thisCA(i,:) = contactAnglesCorrected(i);
end

degIntensityZeroSA = (degIntensityZero);
degIntensityInterestSA = (degIntensityInterest);
intensityRatio = (degIntensityZeroSA./(degIntensityInterestSA));
plot(thisCA,intensityRatio/intensityRatio(61),'-o','lineWidth',2,'Color',[colors(10,:)])
legend('DEB:HFE','TOL:HFE/FC43','DEPH:MPFB')

%% Intensity of the emitting (step 1) phase 
h = figure;
axes1 = axes('Parent',h);
hold on;
%title('Sum of intensity of step 1 emission, 10rd, 10deg cone, different widths')
xlabel('Contact Angle of Droplet');
ylabel('Sum intensity');

%% plotting the intensity of incoming light
for i = 1:size(intensityStorage,2)
    interest = intensityStorage{i};
    ms(i) = sum(sum(interest));
end

plot(contactAnglesCorrected,ms/ms(end),'lineWidth',2,'Color',[colors(10,:)])
xlim([0,90])
legend('DEB:HFE','TOL:HFE/FC43','DEPH:MPFB')

box(axes1,'on');
set(axes1,'FontSize',14);
grid on;
% Set the remaining axes properties
set(axes1,'FontSize',14);

%% Plotting individual droplets
useScreen = logic(1);
useCoverslip = logic(2);
useWaterPlane = logic(3);
returnRays = logic(4);
parLimit = logic(5);
forIntensity = logic(6);
nm = environment(1);
na = environment(2);
nc = environment(3);
ri = environment(4);
minimumAmplitude = environment(5);
propagationDistance = environment(6);
screenHeight = environment(7);
resolution = environment(8);
imgDim = [environment(9) environment(10)];
cTop = planes(1,:);
cBottom = planes(2,:);
waterPlane = planes(3,:);
startPlane = planes(4,:);

for i = 12:12
        figure; hold on;
        Ri = collectedRi(i);
        theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
%         for Ray=rayStorage{i}
%             drawRay(Ray, [0, 0, 0.5]);
%             hold on
%         end
        if useCoverslip
             plot([-4*Rd, 4*Rd], [cTop(1), cTop(1)], 'c');
             plot([-4*Rd, 4*Rd], [cBottom(1), cBottom(1)], 'c');
        end
        if useWaterPlane
            plot([-5*Rd, 5*Rd], [3*Rd, 3*Rd], 'c')
        end
        xlabel('x')
        ylabel('y')
        title(strcat('CA: ',num2str(contactAnglesCorrected(i)),' - Rays:',num2str(size(nRays))))
        set(gcf, 'color', 'white')
        drawDrop(theDrop,100,edgeC,colorO,colorI)        
        xlim([-4*Rd 4*Rd])
        ylim([-4*Rd 4*Rd])
        axis equal

        %scatter(randomRays(1:rayCut,1),randomRays(1:rayCut,2))
        hold off
    end


%%
figure;
hold on
txt = [];
for ii = 1:3
    histStorage = histStorageAll{ii}
    histCut = histStorage(:,2:89)
    dropAngle = dropAngleAll{ii};
    
        [row, col] = find(ismember(histStorage, max(histCut(:))))
        
        %[M, I] = max(histStorage(:,2:89));
        dropData(ii,:) = [contactAnglesCorrected(row), dropAngle];
        
    %plot(dropData(:,1),dropData(:,2),'lineWidth',2)
    %plot(dropDataOnly(:,1),dropDataOnly(:,2),'lineWidth',2)
    txt = [txt num2str(dropAngle)];
end
scatter(dropData(:,1),dropData(:,2))

%legend('0 deg','45deg','90deg')
titlestr = ['Angle of maximum droplet intensity per droplet (ABOVE)']
ylabel('Angle of droplet rotation');
xlabel('Contact angle of droplet');
title(titlestr)

%% Droplet Generator
figure; hold on;
checkAngles = [90];
for i = 1:length(checkAngles)
    for j = 1:length(contactAnglesCorrected)
    if contactAnglesCorrected(j) == checkAngles(i)
            checkStorage(i,:) = intHistStorage(j,:);
            Ri = collectedRi(j);
            theDrop=drop([i*2.2 0], Rd, Ri, no, ni, 1, dropAngle);
            drawDrop(theDrop,100,edgeC,colorO,colorI);
        else
        end
    end
end
xlim([0,length(checkAngles)*2.5])

axis equal

%% circle ...?
checkAngles = [51];
skip = 3;
for i = 1:length(checkAngles)
    for j = 1:length(contactAnglesCorrected)
    if contactAnglesCorrected(j) == checkAngles(i)
            checkStorage(i,:) = intHistStorage(j,:);
            Ri = collectedRi(j);
        else
        end
    end
end
h = figure; hold on; 
[num ang] = size(checkStorage)
angles = linspace(0,359,360)
legendentries = [];
for i = 1:1
    x = angles(1:skip:181);
    y = checkStorage(i,1:skip:181)/sum(checkStorage(i,1:skip:181));
    [xC yC] = pol2cart(x*pi/180,2+(y-min(y))*125)
    plot(xC,yC,'lineWidth',3,'Color',colors(colorselect(3),:),'MarkerIndices',1:5:length(y))
    
    theDrop=drop([0 0], Rd*1.9, Ri*1.9, no, ni, vr, dropAngle);
    drawDrop(theDrop,100,edgeC,colorO,colorI);
end
xlim([-3,3])
ylim([-3,3])
%legend('13', '27', '63', '90')
box on; axis equal; axis off;

%% Checkin Checkin Checkin, droplets over their intensity

checkAngles = [11 26 51 61];
skip = 4;
figure; hold on;
for i = 1:length(checkAngles)
    for j = 1:length(contactAnglesCorrected)
    if contactAnglesCorrected(j) == checkAngles(i)
            checkStorage(i,:) = intHistStorage(j,:);
            Ri = collectedRi(j);
            theDrop=drop([0 -i*2.1], Rd, Ri, no, ni, vr, dropAngle);
            drawDrop(theDrop,100,edgeC,colorO,colorI);
        else
        end
    end
end
axis equal;
hold off;
h = figure; hold on; 
[num ang] = size(checkStorage)
angles = linspace(0,359,360)
legendentries = [];
colorselect = [3,1,7,9]
for i = 1:num
    x = angles(1:skip:181) - 90;
    y = checkStorage(i,1:skip:181)/sum(checkStorage(i,1:skip:181));
    plot(x,y,'lineWidth',2,'Color',colors(colorselect(i),:),'MarkerIndices',1:5:length(y))
end
xlim([-90,90]); set(gca,'FontSize',18);
legend('11', '26', '51', '61')
box on;