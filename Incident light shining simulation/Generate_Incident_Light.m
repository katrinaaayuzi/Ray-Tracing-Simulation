%https://github.com/snnagel/MultiPhaseDropletRaytracer

clear all; %close all; clc;
options = optimset('Display','off');
%parpool('local',16,'IdleTimeout', Inf) %Enable if distributed 
%%

% Droplet parameters & generation of contact angles:
    vr = 1; %volume ratio.  V_internal/v_external
    Rd = 1; %droplet radius
    dropAngle = 0;
    dropLocation = [0,0];
    % Indexes of refraction:
    nm = 1.334; % Surrounding Medium
    ni = 1.462; %inner phase, hydrocarbon
    no = 1.361; % outer phase, fluorocarbon
    ri = nm; %This is the phase the rays start in. VERY IMPORTANT.
    angleRadii = 1;

% How many droplets
%     numPoints = 6; %In this current version, this is the amount of drops between 0 and 90, which gets doubled.
%     contactAngles = [5 10 15 20 25 90]; %This must ALWAYS run with a 90 degree drop, as the baseline
%     %contactAngles = [0 30 60 90];
%     interestAngles = [155 160 165 170 175]; %These are the droplet contact angles we gather intensity maps for
%     reverseAngles = true; %Make true if you want to look at the other side of the droplet conformations    
    CA_start = 90;
    CA_end = 180;
    numPoints = 19; 
    %Modify here
    interestAngles = [90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180]; %These are the droplet contact angles we gather intensity maps for
    reverseAngles = true; %Make true if you want to look at the other side of the droplet conformations

% Raytracing Parameters:
%     emRays = 50000;
    nRays = 100000; %total rays to trace per droplet
    rayPerCore = nRays/5; %rays sent to each core [always less than total rays]
    
%Per-experiment logic
    plotDrops = false; %Do you really want to plot all of the droplets?
    useScreen = false; %This will limit rays to those that hit the screen. Currently unused 
    useCoverslip = true;
    useWaterPlane = false; 
    forIntensity = false; %This makes it so only droplet-and external phase matter for intensity mapping. In this case we only care about direction

%Coverslip (or any glass substrate)
    cThickness = Rd; %thickness of coverslip
    nc = 1.52; %refractive index of coverslip
    na = 1; % refractive index outside of coverslip (air probably)
    %coverslip data (planes that define top and bottom):
    cTop = [-Rd*1.01, 0, 1];
    cBottom = [-Rd*1.01-cThickness, 0, 1];
    %if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip

%Screen info
    screenDistance = -10;
    screenSize = 20;
    
    %plane to propagate rays to after interaction with the droplet
    %if rays don't hit this plane they are thrown away:
    %equation of plane is p(1)=p(2)x+p(3)y.
    screen = [ -screenDistance, 0, 1]; 
    propagationDistance = 1000; %instead propagate all rays after they leave the drop by this amount

%Intensity Settings (For intensity mapping, defines the size of the map you
%generate
    imgDim = [500 500];
    screenHeight = 3;
    startPlane=[-3, 1, 0];
    distance = 10;
    resolution = 200;
    numSlices = resolution*distance;
    minimumAmplitude = 0.01;
    
%Water Plane
    waterPlane = [8*Rd,0,1];

%Misc Settings
    %Droplet color
    colorO = [0.69, 0.94, 1];
    colorI = [1, 0.69, 0.69];
    edgeC  = [0, 0, 0];

%RAY SETTINGS
    
% %%45 degree collimated light source from upper left
% uv(:,1) = linspace(-5*Rd, 0, nRays); 
% uv(:,2) = linspace(0, 5*Rd, nRays);
% rDir = repmat([1, -1],nRays ,1 ); 

%%45 degree collimated light source from upper right
% uv(:,1) = linspace(0, 5*Rd, nRays); %distribution of x coordinates for start
% uv(:,2) = linspace(5*Rd, 0*Rd, nRays); %distribution of y coordinates for start
% rDir = repmat([-1, -1],nRays ,1 ); %directions

%%45 degree collimated light source from down left
% uv(:,1) = linspace(-5*Rd, 0, nRays); 
% uv(:,2) = linspace(0, -5*Rd, nRays);
% rDir = repmat([1, 1],nRays ,1 ); 

%%0 degree collimated light source top to bottom
% uv(:,1) = linspace(-5*Rd, 5*Rd, nRays); 
% uv(:,2) = linspace(5*Rd, 5*Rd, nRays);
% rDir = repmat([0, -1],nRays ,1 ); 

%0 degree collimated light source bottom to top
uv(:,1) = linspace(-5*Rd, 5*Rd, nRays); 
uv(:,2) = linspace(-5*Rd, -5*Rd, nRays);
rDir = repmat([0, 1],nRays ,1 );

%This would be collimated rays, which
%generally, we don't have
%newDir = 258 + 24*rand(emRays,1); %New ray direction generation
%[x,y] = pol2cart(newDir(:,:)*pi/180,1);
%rDir(:,1) = x;
%rDir(:,2) = y;

if rDir == [1, -1]
    lightdir = '45topleft'
elseif rDir == [0, -1]
       lightdir = '0top'
elseif rDir == [0, 1]
       lightdir = '0bottom'
elseif rDir == [1, 1]
       lightdir = '45bottomleft'       
end   

%%Solving for the radii to use for the ray tracer
% collectedRi = zeros(1,numPoints);
% for m = 1 : numPoints
%     contactAngle = contactAngles(m)*pi/180;
%     [solu,fval,exitflag,output] = fsolve(@(x)dropShapeSolverContactAngle(x, angleRadii, vr, cos(contactAngle)), [0.8*angleRadii, 1*angleRadii] );
%     if (exitflag<0)
%         continue
%     end
%     d = solu(1);
%     Ri = solu(2);
%     collectedRi(1,m) = Ri;
% end %Using the solver to find the Ri's for each desired angle
% 
% if reverseAngles == true
%     positiveRi = collectedRi*Rd;
%     negativeRi = -flip(collectedRi)*Rd;
%     positiveAngles = contactAngles;
%     negativeAngles = -flip(contactAngles);
%     contactAnglesCorrected = [contactAngles, 90+contactAngles];
%     collectedRi = [positiveRi, negativeRi];%;, -flip(collectedRi)*Rd];
% else
%     collectedRi = collectedRi*Rd;
%     contactAnglesCorrected = contactAngles;
% end
% numDrops = length(contactAnglesCorrected);

%New Contact Angles generation method

[collectedRi, contactAngles] = generate(CA_start, CA_end, numPoints, vr, Rd);

% Generate Emitters (stage 1)
histStorage = zeros(numPoints,360);
intHistStorage = zeros(numPoints,360);
histsqStorage = zeros(numPoints,360);
intHistsqStorage = zeros(numPoints,360);

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
name = strcat(lightdir,'-',num2str(contactAngles(1)),'-',num2str(contactAngles(length(contactAngles))),'-',num2str(nRays),'rd.mat');

%%    
tic;
for ii = 1:numPoints

    hist = zeros(1,360);
    intHist = zeros(1,360);
    fint = zeros(imgDim(2),imgDim(1));
    emInt = zeros(imgDim(2),imgDim(1));
    if ismember(contactAngles(ii),interestAngles)
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
%         no_corrected = no;
%         ni_corrected = ni;
    %end
    %theDrop = drop(dropLocation, Rd, Ri, no_corrected, ni_corrected, vr, dropAngle);
    theDrop = drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
    rayNum = [];
	
    locsCut = uv;
    rayChop = nRays/rayPerCore;

    parfor l = 1:rayChop
        if l == rayChop
            choppedRays = locsCut((l-1)*rayPerCore+1:nRays,:);
            choppedDirs = rDir((l-1)*rayPerCore+1:nRays,:);
            rayNum = length(choppedDirs);
        else
            choppedRays = locsCut((l-1)*rayPerCore+1:(l*rayPerCore+1),:);
            choppedDirs = rDir((l-1)*rayPerCore+1:(l*rayPerCore+1),:);
            rayNum = length(choppedDirs);
        end               
        
        %Question, how does this work? and can i make a program which will
        %get me all of the intensity data (water, hydrocarbon,
        %fluorocarbon) into separate arrays from the same sequence?
        
        %suggestion from brad: make my own 'multidroppartrace' which will
        %take the ray data and make me this information so i can plot it
        %myself
        
        %a separate project: is make this work with the 'spatial emitter
        %generation' which simulations a fluorocarbon dye
        
        %OUTSIDE DROPLET SPECIAL MADE
        %Question: How does this work? Can i make it also plot intensity
        %differences in the fluorocarbon phase?
        
        cintTemp = emitterGeneration(rayNum,choppedRays,theDrop,choppedDirs,environment); %INSIDE DROPLET SPECIAL MADE BY 
        emInt = emInt + cintTemp;
        
        [chist,cIntHist,cint] = multiDropParTrace(rayNum,choppedRays,theDrop,choppedDirs,logic,environment,planes);
        hist = hist+chist; %this is the 'leaving every other interface we care about' version
        intHist = intHist+cIntHist; %this is the 'just leaving the droplet' version
        fint = fint + cint; %these are maps 
        
    end
    
	dirStorage{ii} = rDir;
    randomRays = []; tempRays = [];
    histStorage(ii,:) = hist;
    emIntensityStorage{ii} = emInt;
    fIntensityStorage{ii} = fint;
    intHistStorage(ii,:) = intHist;
    
    save(name,'contactAngles','collectedRi','dropLocation','Rd','Ri','no','ni','vr','dropAngle','histStorage','numPoints','environment','logic','planes','nRays','imgDim','intHistStorage','emIntensityStorage','fIntensityStorage')
    forIntensity = false;
     ii
end
toc;
%%
% %Plotting our test maps how fun1!!
% 
% for i = 1:length(emIntensityStorage)
% figure; hold on;
% title([contactAnglesCorrected(i)])
% titlestr = ['CA_d=' num2str(contactAnglesCorrected(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
% Ri = collectedRi(i);
% theDrop=drop(dropLocation, Rd, Ri, no, ni, vr);
% %drawDrop(theDrop)
% %xlim([0 500])
% %ylim([0 500])
% %caxis([0 5])
% colormap(inferno)
% imagesc(emIntensityStorage{i})
% axis image;
% 
% end
% 
% 
%%
% for i = 1:size(intensityStorage,2)
    % point(1,i) = sum(sum(emIntensityStorage{i}))
% end

%%
% figure; axis image;
% drawDrop(theDrop)
% scatter(uv(:,1),uv(:,2))
