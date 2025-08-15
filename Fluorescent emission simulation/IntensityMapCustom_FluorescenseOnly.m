%%
output = addpath('output/');
colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];
colors = [0.6196 0.0039 0.2588 ; 0.9353 0.2431 0.3098 ; 0.9569 0.4275 0.2627 ; 0.9922 0.6824 0.3804 ; 0.9961 0.8784 0.5451 ; 0.9961 0.8784 0.5451 ; 0.6706 0.8667 0.6431 ; 0.4000 0.7608 0.6471 ; 0.1961 0.5333 0.7412 ; 0.3686 0.3098 0.6353];


%% Processing custom FLUORESCENCE DROPLETS! 2025 Katrina: for fluorescent dye emission from droplets with fluorescent dye in HC

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
        caxis([0 max(intCheck)])
        %caxis([0 3000]) % caxis([0 10000]), caxis([0 4500])
       
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

%% modified custom FLUORESCENCE DROPLETS intenisty maps with same scale bar nomralized to the maximum scale bar

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];

% Determine the global maximum intensity across all intensity maps
maxIntensity = max(cellfun(@(x) max(x(:)), fIntensityStorage));
minIntensity = min(cellfun(@(x) min(x(:)), fIntensityStorage));

for i = 1:length(fIntensityStorage)
    interest = fIntensityStorage{i};
    if sum(sum(interest)) > 0
        figure; hold on; axis image; axis off;

        Ri = collectedRi(i);
        xlim([0 imgDim(1)]);
        ylim([0 imgDim(2)]);
        imagesc(interest);
        colorbar;
        caxis([minIntensity, maxIntensity]); % Apply consistent scale to all plots
        
        titlestr = ['FEM VR=' num2str(vr) '-CA=' num2str(contactAngles(i)) ...
                    '-Solution(n=' num2str(environment(1)) ')-HC(n=' num2str(ni) ...
                    ')-FC(n=' num2str(no) ')'];
        title(titlestr);
        colormap(viridis);

        % Intensity of emission (optional code for visualization or export)
        % x = linspace(1,180,180)*(500/180);
        % checkStorage = fIntensityStorage(i,:);
        % y = (smooth(checkStorage(:,1:180),10)/mean(checkStorage(:,1:180)))*250-200;
        % export_fig(titlestr,'-png','-r300')
        % plot(x, y, 'lineWidth', 2, 'Color', colors(7,:))

        xlim([0, 500]);
        xlabel('Intensity around top hemisphere of droplet');
        
        % Draw the drop
        theDrop = drop(dropLocation + imgDim / 2, Rd * imgDim(1) / 2.5, ...
                       Ri * imgDim(1) / 2.5, no, ni, vr, dropAngle);
        drawDrop(theDrop, 1000, edgeC, colorO, colorI);
    end
end

%% Processing custom FLUORESCENCE Control Droplet! 

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];

    intCheck = max(max(fIntensityStorage2)); % KATRINA 2025: what is fIntensityStorage2, only for fluorescent dye containing droplets emission result
    luckydrop = 90;

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
