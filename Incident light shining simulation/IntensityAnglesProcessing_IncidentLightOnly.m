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

%% Ratiometric bottom readout vs contact angles written by katrina copyright 2025
    
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

%% Scanning intensity profile of a droplet from bottom 0 to 45(BinAngles 270 to 360)

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
 rectangle('position',[260 0 20 180])
 rectangle('position',[215 0 20 180])

%%
%this is just as light leaves the droplet enter water before air

%interest = histStorage;
%intHisStorage is before Snell's law, hisStorage is after Snell's law
interest = intHistStorage;

figure;
titlestr = ['VR=' num2str(vr) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) '),Rays Leaving Droplet, Color is Intensity'];
[X,Y] = meshgrid(1:1:360,contactAngles);
s = surf(X(:,1:360),Y(:,1:360),intHistStorage(:,1:360)); %use pcolor instead of surf for 2d
s.EdgeColor = 'none';
colormap(viridis)
xlim([180 360])
xlabel('Angle Ray Is Heading After Leaving Droplet(\circ)')
ylabel('Contact Angle of Droplet(\circ)')
zlabel('Intensity of Rays(a.u.)')
title(titlestr)
colorbar
