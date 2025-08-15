clear all; %close all; clc;
options = optimset('Display','off');
%parpool('local',16,'IdleTimeout', Inf) %Enable if distributed 
%%

% Droplet parameters & generation of contact angles:
    vr=1; %volume ratio.  V_internal/v_external
    Rd=1; %droplet radius
    dropAngle = 0;
    dropLocation=[0,0];
    % Indexes of refraction:
    nm=1.334; % Surrounding Medium
    ni=1.462; %inner phase, hydrocarbon
    no=1.361; % outer phase, fluorocarbon
    ri=nm; %This is the phase the rays start in. VERY IMPORTANT.
    angleRadii = 1;

% How many droplets
    numPoints = 4; %In this current version, this is the amount of drops between 0 and 90, which gets doubled.
    contactAngles = [90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180]; %This must ALWAYS run with a 90 degree drop, as the baseline
    %contactAngles = [10 27 35 50 65 85];
    interestAngles = [90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180]; %These are the droplet contact angles we gather intensity maps for
    reverseAngles=true; %Make true if you want to look at the other side of the droplet conformations

% Raytracing Parameters:
    nRays = 10000; %total rays to trace per droplet
    rayPerCore = nRays/5; %rays sent to each core [always less than total rays]
    
%Per-experiment logic
    plotDrops = false; %Do you really want to plot all of the droplets?
    useScreen = false; %This will limit rays to those that hit the screen. Currently unused 
    useCoverslip = false;
    useWaterPlane = true; 
    forIntensity = false; %This makes it so only droplet-and external phase matter for intensity mapping. In this case we only care about direction

%Coverslip (or any glass substrate)
    cThickness=Rd; %thickness of coverslip
    nc=1.52; %refractive index of coverslip
    na=1; % refractive index outside of coverslip (air probably)
    %coverslip data (planes that define top and bottom):
    cTop=[-Rd*1.01, 0, 1];
    cBottom=[-Rd*1.01-cThickness, 0, 1];
    %if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip

%Screen info
    screenDistance=-10;
    screenSize=20;
    
    %plane to propagate rays to after interaction with the droplet
    %if rays don't hit this plane they are thrown away:
    %equation of plane is p(1)=p(2)x+p(3)y.
    screen=[ -screenDistance, 0, 1]; 
    propagationDistance=1000; %instead propagate all rays after they leave the drop by this amount

%Intensity Settings (For intensity mapping, defines the size of the map you
%generate
    imgDim = [500 500];
    screenHeight = 3;
    startPlane=[-3, 1, 0];
    distance = 10;
    resolution=200;
    numSlices=resolution*distance;
    minimumAmplitude=0.01;
    
%Water Plane
waterPlane = [8*Rd,0,1];

%Misc Settings
    %Droplet color
    colorO = [0.69, 0.94, 1];
    colorI = [1, 0.69, 0.69];
    edgeC  = [0, 0, 0];

%RAY SETTINGS

% %%45 degree collimated light source from left
% uv(:,1)=linspace(-5*Rd, 0, nRays); 
% uv(:,2)=linspace(0, 5*Rd, nRays);
% rDir=repmat([1, -1],nRays ,1 ); 

%%45 degree collimated light source from right
% uv(:,1)=linspace(0, 5*Rd, nRays); %distribution of x coordinates for start
% uv(:,2)=linspace(5*Rd, 0*Rd, nRays); %distribution of y coordinates for start
% rDir=repmat([-1, -1],nRays ,1 ); %directions

%%0 degree collimated light source top to bottom
% uv(:,1)=linspace(-5*Rd, 5*Rd, nRays); 
% uv(:,2)=linspace(5*Rd, 5*Rd, nRays);
% rDir=repmat([0, -1],nRays ,1 ); 

%%0 degree collimated light source bottom to top
uv(:,1)=linspace(-5*Rd, 5*Rd, nRays); 
uv(:,2)=linspace(-5*Rd, -5*Rd, nRays);
rDir=repmat([0, 1],nRays ,1 );

%This would be collimated rays, which
%generally, we don't have
%newDir = 258 + 24*rand(emRays,1); %New ray direction generation
%[x,y] = pol2cart(newDir(:,:)*pi/180,1);
%rDir(:,1) = x;
%rDir(:,2) = y;

%%Solving for the radii to use for the ray tracer
collectedRi = zeros(1,numPoints);
for m = 1 : numPoints
    contactAngle=contactAngles(m)*pi/180;
    [solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, angleRadii, vr, cos(contactAngle)), [0.8*angleRadii, 1*angleRadii] );
    if (exitflag<0)
        continue
    end
    d=solu(1);
    Ri=solu(2);
    collectedRi(1,m) = Ri;
end %Using the solver to find the Ri's for each desired angle

if reverseAngles == true
    positiveRi = collectedRi*Rd;
    negativeRi = -flip(collectedRi)*Rd;
    positiveAngles = contactAngles;
    negativeAngles = -flip(contactAngles);
    contactAnglesCorrected = [contactAngles, 90+contactAngles];
    collectedRi = [positiveRi, negativeRi];%;, -flip(collectedRi)*Rd];
else
    collectedRi = collectedRi*Rd;
    contactAnglesCorrected = contactAngles;
end
numDrops = length(contactAnglesCorrected);

% Generate Emitters (stage 1)
histStorage = zeros(numDrops,360);
intHistStorage = zeros(numDrops,360);
histsqStorage = zeros(numDrops,360);
intHistsqStorage = zeros(numDrops,360);

%You can kind of ignore this, it is only setting the information for
%parralellization
logic = [useScreen useCoverslip useWaterPlane forIntensity];
planes = [cTop;cBottom;waterPlane;startPlane;screen];
hist = zeros(1,360);
intHist = zeros(1,360);
histsq = zeros(1,360);
intHistsq = zeros(1,360);
cint = zeros(imgDim(2),imgDim(1));
emInt = zeros(imgDim(2),imgDim(1));
cIntHist = zeros(screenHeight*2*resolution, numSlices);
eint = zeros(imgDim(2),imgDim(1));
environment = [nm na nc ri minimumAmplitude propagationDistance screenHeight resolution imgDim no];
name = strcat(num2str(contactAnglesCorrected(1)),'-',num2str(contactAnglesCorrected(length(contactAnglesCorrected))),'-',num2str(nRays),'rd.mat');

%%    
tic;
for ii = 1:numDrops

    hist = zeros(1,360);
    intHist = zeros(1,360);
    emInt = zeros(imgDim(2),imgDim(1));
    if ismember(contactAnglesCorrected(ii),interestAngles)
        logic(4) = true; %sets forintensity to be true sometimes, so we can run only a few int maps
    else 
        logic(4) = false;
    end
    
    Ri = collectedRi(ii);
    theDrop = [];
    
%     if contactAnglesCorrected(ii) > 90
%         no_corrected = ni;
%         ni_corrected = no;
%     else
        no_corrected = no;
        ni_corrected = ni;
    %end
    theDrop = drop(dropLocation, Rd, Ri, no_corrected, ni_corrected, vr, dropAngle);
    rayNum = [];
	
    locsCut = uv;
    rayChop = nRays/rayPerCore;

    for l = 1:rayChop
        if l == rayChop
            choppedRays = locsCut((l-1)*rayPerCore+1:nRays,:);
            choppedDirs = rDir((l-1)*rayPerCore+1:nRays,:);
            rayNum = length(choppedDirs);
        else
            choppedRays = locsCut((l-1)*rayPerCore+1:(l*rayPerCore+1),:);
            choppedDirs = rDir((l-1)*rayPerCore+1:(l*rayPerCore+1),:);
            rayNum = length(choppedDirs);
        end
%         [chist,cIntHist,cint] = multiDropParTrace(rayNum,choppedRays,theDrop,choppedDirs,logic,environment,planes);
%         hist = hist+chist;
%         intHist = intHist+cIntHist;
        %Question, how does this work? and can i make a program which will
        %get me all of the intensity data (water, hydrocarbon,
        %fluorocarbon) into separate arrays from the same sequence?
        
        %suggestion from brad: make my own 'multidroppartrace' which will
        %take the ray data and make me this information so i can plot it
        %myself
        
        %a separate project: is make this work with the 'spatial emitter
        %generation' which simulations a fluorocarbon dye
        
        %OUTSIDE DROPLET SPECIAL MADE
        emIntTemp = emitterGeneration(rayNum,choppedRays,theDrop,choppedDirs,environment); %INSIDE DROPLET SPECIAL MADE BY 
        %Question: How does this work? Can i make it also plot intensity
        %differences in the fluorocarbon phase?

        emInt = emInt + emIntTemp;
    end
	dirStorage{ii} = rDir;
    randomRays = []; tempRays = [];
    histStorage(ii,:) = hist;
    emIntensityStorage{ii} = emInt;
    intHistStorage(ii,:) = intHist;
    
    save(name,'contactAnglesCorrected','collectedRi','dropLocation','Rd','Ri','no','ni','vr','dropAngle','histStorage','numDrops','environment','logic','planes','nRays','imgDim','intHistStorage','emIntensityStorage')
    forIntensity = true;
     ii
end
toc;


%%
output = addpath('output/');
colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];
edgeC=[0, 0, 0];
colors = [0.6196 0.0039 0.2588 ; 0.9353 0.2431 0.3098 ; 0.9569 0.4275 0.2627 ; 0.9922 0.6824 0.3804 ; 0.9961 0.8784 0.5451 ; 0.9961 0.8784 0.5451 ; 0.6706 0.8667 0.6431 ; 0.4000 0.7608 0.6471 ; 0.1961 0.5333 0.7412 ; 0.3686 0.3098 0.6353];

%% Processing custom emission intensity maps

colorO = [0.74, 0.72, 0.8];
colorI = [1, 0.69, 0.69];

current = emIntensityStorage;
for i =  1:length(emIntensityStorage)
    intCheck(i) = max(max(current{i}));
end
max(intCheck)
total = length(current);
for i =  1:total
    interest = current{i};
    if sum(sum(interest)) > 0
        figure; hold on; axis image; axis off;

        Ri = collectedRi(i);
        xlim([0 imgDim(1)])
        ylim([0 imgDim(2)])
        imagesc(interest);
         colorbar;
        caxis([0 6])
       
        titlestr = ['EM CA_d=' num2str(contactAnglesCorrected(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
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
        theDrop=drop(dropLocation + imgDim/2, Rd*imgDim(1)/2.5, Ri*imgDim(1)/2.5, no, ni, vr,dropAngle);
        %drawDrop(theDrop,1000,edgeC,colorO,colorI);

    end
end


%% emIntensity && Intensity on top of eachother

for i = 1:10%(length(intensityStorage))
    interest = intensityStorage{i};
    if sum(sum(interest)) > 0

        Ri = collectedRi(i);
        theDrop=drop(dropLocation + imgDim/2, Rd*imgDim(1)/5, Ri*imgDim(1)/5, no, ni, vr,dropAngle);
        
        interestInt = intensityStorage{i};
        intAlpha = zeros(imgDim(1),imgDim(2));
        intAlpha(interestInt < 50) = 1;
        
       
        hf = figure;
        h1 = axes;
        colormap(h1,'parula');
        p1=imagesc(interest);
        set(h1,'ydir','normal');
        
        h2 = axes;
        p2 = imagesc(interestInt);
        set(p2, 'Alphadata', intAlpha);
        colormap(h2,'jet');
        set(h2,'ydir','normal');
                linkaxes([h1 h2])
        %xlim([0 imgDim(1)])
        %ylim([0 imgDim(2)])
        %hold on;
         
        %axis image
        titlestr = ['CA_d=' num2str(contactAnglesCorrected(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
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
        
        titlestr = ['CA_d=' num2str(contactAnglesCorrected(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
        %title(titlestr);
        colorbar
        colormap(cmappablo)
        %caxis([0 max(max(intensityStorage{:,:}))])
        %caxis ([0 4.5])
        %set(gca,'ColorScale','log')
        %export_fig(titlestr,'-png','-r300','-p0.02')
        g = drawDrop(theDrop,1000,edgeC,colorO,colorI);
        norm = thismap/mean((mean(thismap(1:10,1:10))));
        h = imagesc(norm);
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
scatter(contactAnglesCorrected(1),mean(mean(intensityStorage{1})))

%axes1 = axes('Parent',h);
ylabel('Intensity (a.u.)');
xlabel('Droplet Cross Section');
%legend('{\theta} = 1{\circ}','{\theta} = 10{\circ}','{\theta} = 26{\circ}','{\theta} = 45{\circ}','{\theta} = 90{\circ}')