%%
output = addpath('output/');
colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];
colors = [0.6196 0.0039 0.2588 ; 0.9353 0.2431 0.3098 ; 0.9569 0.4275 0.2627 ; 0.9922 0.6824 0.3804 ; 0.9961 0.8784 0.5451 ; 0.9961 0.8784 0.5451 ; 0.6706 0.8667 0.6431 ; 0.4000 0.7608 0.6471 ; 0.1961 0.5333 0.7412 ; 0.3686 0.3098 0.6353];

%% Processing custom INTER DROPLET EMISSION MAPS!

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];

intensityStorage = emIntensityStorage;
for i =  1:length(emIntensityStorage)
    intCheck(i) = max(max(intensityStorage{i}));
end
max(intCheck)

for i = 1:length(intensityStorage)
    interest = intensityStorage{i};
    if sum(sum(interest)) > 0
        figure; hold on; axis image; axis off;

        Ri = collectedRi(i);
        xlim([0 imgDim(1)])
        ylim([0 imgDim(2)])
        imagesc(interest);
        colorbar;
        caxis([0 100])
       
        titlestr = ['EM VR=' num2str(vr) '-CA=' num2str(contactAngles(i)) '-cont(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) ')'] ;
        title(titlestr);
        colormap(inferno)

        %Intensity of emission
        x = linspace(1,180,180)*(500/180);
        checkStorage = intHistStorage(i,:);
        y = (smooth(checkStorage(:,1:180),10)/mean(checkStorage(:,1:180)))*250-200;
        %export_fig(titlestr,'-png','-r300')
        %plot(x,y,'lineWidth',2,'Color',colors(7,:))

        xlim([0,500]); set(gca,'FontSize',18);  
        xlabel('Intensity around top hemisphere of droplet');
        theDrop = drop(dropLocation + imgDim/2, Rd*imgDim(1)/2.5, Ri*imgDim(1)/2.5, no, ni, vr,dropAngle);
        %drawDrop(theDrop,1000,edgeC,colorO,colorI);

    end
end

%% Processing custom FLUORESCENCE DROPLETS!

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];

for i =  1:length(fIntensityStorage)
    intCheck(i) = max(max(fIntensityStorage{i}));
end

for i = 1:length(fIntensityStorage)
    interest = fIntensityStorage{i};
    if sum(sum(interest)) > 0
        figure; hold on; axis image; axis off;

        Ri = collectedRi(i);
        xlim([0 imgDim(1)])
        ylim([0 imgDim(2)])
        imagesc(interest);
        colorbar;
        %caxis([0 max(intCheck)])
        %caxis([0 10000])
       
        titlestr = ['FEM VR=' num2str(vr) '-CA=' num2str(contactAngles(i)) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) ')'] ;
        title(titlestr);
        colormap(viridis)

        %Intensity of emission
        %x = linspace(1,180,180)*(500/180);
        %checkStorage = fIntensityStorage(i,:);
        %y = (smooth(checkStorage(:,1:180),10)/mean(checkStorage(:,1:180)))*250-200;
        %export_fig(titlestr,'-png','-r300')
        %plot(x,y,'lineWidth',2,'Color',colors(7,:))

        xlim([0,500]); %set(gca,'FontSize',18);  
        xlabel('Intensity around top hemisphere of droplet');
        theDrop = drop(dropLocation + imgDim/2, Rd*imgDim(1)/2.5, Ri*imgDim(1)/2.5, no, ni, vr,dropAngle);
        drawDrop(theDrop,1000,edgeC,colorO,colorI);

    end
end

%% Processing custom FLUORESCENCE Control Droplet!

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];

    intCheck = max(max(fIntensityStorage2));
    %luckydrop = 90;

    if sum(sum(fIntensityStorage2)) > 0
        figure; hold on; axis image; axis off;

        Ri = collectedRi(luckydrop);
        xlim([0 imgDim(1)])
        ylim([0 imgDim(2)])
        imagesc(fIntensityStorage2);
        colorbar;
        %caxis([0 max(intCheck)])
        caxis([0 10000])
       
        titlestr = ['FEM VR=' num2str(vr) '-CA=' num2str(contactAngles(luckydrop)) '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ')-FC(n=' num2str(no) ')'] ;
        title(titlestr);
        colormap(viridis)

        %Intensity of emission
        %x = linspace(1,180,180)*(500/180);
        %checkStorage = fIntensityStorage(i,:);
        %y = (smooth(checkStorage(:,1:180),10)/mean(checkStorage(:,1:180)))*250-200;
        %export_fig(titlestr,'-png','-r300')
        %plot(x,y,'lineWidth',2,'Color',colors(7,:))

        xlim([0,500]); %set(gca,'FontSize',18);  
        xlabel('Intensity around top hemisphere of droplet');
        theDrop = drop(dropLocation + imgDim/2, Rd*imgDim(1)/2.5, Ri*imgDim(1)/2.5, no, ni, vr,dropAngle);
        %drawDrop(theDrop,1000,edgeC,colorO,colorI);

    end

%% Emitter position ploting

figure;
u = uv(:,1);
v = uv(:,2);
scatter(u,v);

%% emIntensity && Intensity on top of each other

for i = 1:10:(length(intensityStorage))
    interest = intensityStorage{i};
    if sum(sum(interest)) > 0

        Ri = collectedRi(i);
        theDrop = drop(dropLocation + imgDim/2, Rd*imgDim(1)/5, Ri*imgDim(1)/5, no, ni, vr,dropAngle);
        
        interestInt = intensityStorage{i};
        intAlpha = zeros(imgDim(1),imgDim(2));
        intAlpha(interestInt < 50) = 1;        
       
        hf = figure;
        h1 = axes;
        colormap(inferno);
        p1=imagesc(interest);
        set(h1,'ydir','normal');
        
        h2 = axes;
        p2 = imagesc(interestInt);
        set(p2, 'Alphadata', intAlpha);
        colormap(inferno);
        set(h2,'ydir','normal');
                linkaxes([h1 h2])
        colorbar;
        caxis([0 6])
        %xlim([0 imgDim(1)])
        %ylim([0 imgDim(2)])
        %hold on;
         
        %axis image
        titlestr = ['CA_d=' num2str(contactAngles(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
        title(titlestr);

    end
end

%% Processing custom intensity maps for droplet tracing
%colorO = [0.5, 0.5, 0.5];
colorO = [0.94, 0.94, 0.94];
colorI = [0.94, 0.94, 0.94];

for i = 1:length(intensityStorage)
    if sum(sum(intensityStorage{i})) > 0
        interest = intensityStorage{i};
        thismap = interest(interest >= 6000);
        %interest = interest/min(min(interestnonNull));
        figure;
        %imgDim = [500 500];
        hold on
        Ri = collectedRi(i);
        theDrop = drop(dropLocation + imgDim/2, Rd*imgDim(1)/2.5, Ri*imgDim(1)/2.5, no, ni, vr,dropAngle);
        %xlim([0 imgDim(1)])
        %ylim([0 imgDim(2)])
        axis image
        
        titlestr = ['CA_d=' num2str(contactAngle(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
        %title(titlestr);
        colorbar
        colormap(inferno)
        %caxis([0 max(max(intensityStorage{:,:}))])
        %caxis ([0 4.5])
        %set(gca,'ColorScale','log')
        %export_fig(titlestr,'-png','-r300','-p0.02')
        g = drawDrop(theDrop,1000,edgeC,colorO,colorI);
        h=imagesc(thismap)
%         norm = thismap/mean((mean(thismap(1:10,1:10))));
%         h = imagesc(norm);
        %g = drawDrop(theDrop,1000,edgeC,colorO,colorI)

        %h.AlphaData = 0.8;

        %xlim([0 500])
        %ylim([0 500])
        %axis off
    end
end

%% A little section to output some sort of confocal crossection of a traced droplet
h = figure; hold on;
% Set the remaining axes properties

%% HERE IS WHERE WE DO THE DROPLET TRACING
h = figure; hold on;
set(gca,'FontSize',18); box on; grid on;
droplet = intensityStorage{1};
normal = mean(mean(droplet(245:255,320:330)));
interest = droplet(51:449,240:260)';
interestMean = mean(interest,1);
x = linspace(0,1,size(interestMean,2));
y = interestMean/normal;
%y = interestMean/(mean(droplet(500,1:50)))
%y = (interestMean)/(mean(droplet(500,:)))
plot(x,y,'Color',colors(9,:),'lineWidth',3,'Color',[colors(1,:)])


%axes1 = axes('Parent',h);
ylim([0.5 2])
ylabel('Intensity (a.u.)');
xlabel('Droplet Cross Section');
%legend('{\theta} = 1{\circ}','{\theta} = 10{\circ}','{\theta} = 26{\circ}','{\theta} = 45{\circ}','{\theta} = 90{\circ}')

%% Bottom!
h = figure; hold on;
set(gca,'FontSize',18); box on; grid on;
droplet = intensityStorage{4};
normal = mean(mean(droplet(245:255,320:330)));
interest = droplet(240:250,1:500);
interestMean = mean(interest,1);
x = linspace(0,1,size(interestMean,2));
y = interest/normal;
%y = interestMean/(mean(droplet(500,1:50)))
%y = (interestMean)/(mean(droplet(500,:)))
plot(x,y,'Color',colors(9,:),'lineWidth',3,'Color',[colors(1,:)])
%scatter(x,y)

%axes1 = axes('Parent',h);
%ylim([0.5 2])
ylabel('Intensity (a.u.)');
xlabel('Droplet Cross Section');
%legend('{\theta} = 1{\circ}','{\theta} = 10{\circ}','{\theta} = 26{\circ}','{\theta} = 45{\circ}','{\theta} = 90{\circ}')


%%
h = figure; hold on;
%%
set(gca,'FontSize',18); box on; grid on;
droplet = intensityStorage{end};
interest = droplet(300:320,:);
interestMean = mean(interest,1);
x = linspace(0,1,size(interestMean,2));
y = interestMean/mean(interestMean(245:255));
%y = interestMean/(mean(droplet(500,1:50)))
%y = (interestMean)/(mean(droplet(500,:)))
scatter(contactAngles(1),mean(mean(intensityStorage{1})))

%axes1 = axes('Parent',h);
ylabel('Intensity (a.u.)');
xlabel('Droplet Cross Section');
%legend('{\theta} = 1{\circ}','{\theta} = 10{\circ}','{\theta} = 26{\circ}','{\theta} = 45{\circ}','{\theta} = 90{\circ}')