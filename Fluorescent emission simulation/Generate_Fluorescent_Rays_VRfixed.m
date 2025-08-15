%https://github.com/snnagel/MultiPhaseDropletRaytracer

clear all; close all; clc;
options = optimset('Display','off');
%parpool('local',16,'IdleTimeout', Inf) %Enable if distributed 
%%

% Droplet parameters & generation of contact angles:
    vr = 0.8; %volume ratio.  V_internal/v_external
    Rd = 1; %droplet radius
    dropAngle = 0;
    dropLocation = [0,0];
    % Indexes of refraction:
    nm = 1.334; % Surrounding Medium
    ni = 1.462; %inner phase, hydrocarbon
    no = 1.361; % outer phase, fluorocarbon
    ri = ni; %This is the phase the rays start in. VERY IMPORTANT.
    angleRadii = 1;

% How many droplets
    CA_start = 90;
    CA_end = 180;
    numPoints = 91; %In this current version, this is the amount of drops between 0 and 90, which gets doubled.
    
    %Modify here
    interestAngles = [90 100 110 120 130 140 150 160 170 180]; %These are the droplet contact angles we gather intensity maps for
    reverseAngles = true; %Make true if you want to look at the other side of the droplet conformations

% Raytracing Parameters:
%     emRays = 50000;
    nRays = 1000000; %total rays to trace per droplet
    rayPerCore = 10000; %rays sent to each core [always less than total rays]
    
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

%heres where we generate all of the droplets
[collectedRi, contactAngles] = generate(CA_start, CA_end, numPoints, vr, Rd);
%%
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
name = strcat('VR=',num2str(vr),'-',num2str(numPoints),'-numPoints',num2str(contactAngles(1)),'-',num2str(contactAngles(length(contactAngles))),'-',num2str(nRays),'rd.mat');

%%    
tic;

%figure; hold on;
%checkAngles = [1];
%for i = 1:length(checkAngles)

    
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
    
    theDrop = drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
    rayNum = [];
    
    %This is the code for generating equally spaced emitters within the
    %hydrocarbon phase of the droplet which is necessary for simulating a
    %'1 step' fluorescence, where incident light total internal reflect is
    %not the vibe
    
        tempFig = figure('visible','off'); %We take the droplet and get the polygon out of it
        [fO, fI] = drawDrop(theDrop);
        if Ri >= 0
            f = fO;
        else
            f = fI;
        end
        close(tempFig)
        polyGen = polyshape({f(1,:),f(2,:)},{f(3,:),f(4,:)}); %Generate polygon from points
        triPoly = triangulation(polyGen); %Generate triangles from this polygon
        a = zeros(1,length(triPoly.ConnectivityList));
        for j = 1:length(triPoly.ConnectivityList) %Generate a list of triangle areas
            [p] = triPoly.ConnectivityList(j,:);
            a(j) = 1/2*abs(det([triPoly.Points(p(1),1),triPoly.Points(p(1),2),1;triPoly.Points(p(2),1),triPoly.Points(p(2),2),1;triPoly.Points(p(3),1),triPoly.Points(p(3),2),1]));
        end
        a = a/sum(a); %normalize triangles
        [~,~,simpind] = histcounts(rand(nRays,1),cumsum([0,a])); %Next we find out how many points each triangle should get,
        r1 = rand(nRays,1); %we randomly assign points to each triangle
        uv = triPoly.Points(triPoly.ConnectivityList(simpind,1),:).*r1 + triPoly.Points(triPoly.ConnectivityList(simpind,2),:).*(1-r1);
        r2 = sqrt(rand(nRays,1));
        uv = uv.*r2 + triPoly.Points(triPoly.ConnectivityList(simpind,3),:).*(1-r2); %and we generate points from this
        %rDir = -1 + (1+1)*rand(nRays,2); %Generate random ray directions, combined with rays later
        
        newDir = 360*rand(nRays,1); %New ray direction generation
        [x,y] = pol2cart(newDir(:,:)*pi/180,1);
        rDir(:,1) = x;
        rDir(:,2) = y;

    
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
        
        [chist,cIntHist,cint] = multiDropParTrace(rayNum,choppedRays,theDrop,choppedDirs,logic,environment,planes);
        hist = hist+chist; %this is the 'leaving every other interface we care about' version
        intHist = intHist+cIntHist; %this is the 'just leaving the droplet' version
        fint = fint + cint; %these are maps 
        
        
        %Question, how does this work? and can i make a program which will
        %get me all of the intensity data (water, hydrocarbon,
        %fluorocarbon) into separate arrays from the same sequence?
        
        %suggestion from brad: make my own 'multidroppartrace' which will
        %take the ray data and make me this information so i can plot it
        %myself
        
        %a separate project: is make this work with the 'spatial emitter
        %generation' which simulations a fluorocarbon dye
        
        %Question: How does this work? Can i make it also plot intensity
        %differences in the fluorocarbon phase?
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
% title([contactAngles(i)])
% titlestr = ['CA_d=' num2str(contactAngles(i)) '-cont_n=' num2str(environment(1)) '-HC_n=' num2str(ni) '-FC_n=' num2str(no)] ;
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
