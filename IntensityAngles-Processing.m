%% 
clear all; close all; clc;
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
title(titlestr)
titlestr = ['Degrees of bottom readout']
plot(contactAngles,intensityZero,'DisplayName','0-degrees-intensity')
plot(contactAngles,intensity45,'DisplayName','45-degrees-intensity')

xlabel('Contact Angles(\Theta\circ)')
ylabel('Intensity(a.u.)')

%ylim([0.75 1.5])
legend

%0/45
ratiomatricIntensity = intensityZero./intensity45;
%xlim([135 180])
subplot(2,1,2);
%subplot(1,1,1);
title(titlestr)
titlestr = ['Ratiometric bottom readout']
plot(contactAngles,ratiomatricIntensity/ratiomatricIntensity(1))
xlabel('Contact Angles(\Theta\circ)')
%ylim([0.75 1.5])
ylabel('Ratiometric 0\circ/45\circ intensity (a.u.)')

%% Scanning intensity profile of a droplet from bottom 0 to 45(BinAngles 270 to 315)

variousAnglesIntensity = [];
figure; hold on;
titlestr = ['Scanning intensity profile of a droplet from bottom 0\circ to 70\circ' '-VR=' num2str(vr)];
title(titlestr)
xlabel('Detection angles(\circ)')
ylabel('Intensity(a.u.)')
%selectDroplets = [51  63  75  87 99 111 123];
selectDroplets = [40 60 80 100 120 140 160 180];
%selectDroplets = [20 32 40 52 60 72];

normalization = (histStorage2/max(histStorage2));

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