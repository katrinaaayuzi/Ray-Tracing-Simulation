%https://github.com/snnagel/MultiPhaseDropletRaytracer

clear all; close all; clc;
options = optimset('Display', 'off');

%% Droplet parameters
vr = 1; % Volume ratio (V_internal/V_external)
Rd = 1; % Droplet radius
nm = 1.334; % Surrounding medium index of refraction
ni = 1.462; % Inner phase index of refraction
no = 1.361; % Outer phase index of refraction
ri = nm; % Starting phase (VERY IMPORTANT)

dropAngle = 0;
dropLocation = [0, 0];
angleRadii = 1;

% Contact angle settings
CA_start = 90;
CA_end = 180;
numPoints = 19;
interestAngles = 90:5:180;
reverseAngles = true;

% Ray tracing parameters
nRays = 10000;
rayPerCore = nRays / 5; % Rays sent to each core

% Miscellaneous settings
plotDrops = false;
useScreen = false;
useCoverslip = true;
useWaterPlane = false;
forIntensity = false;

% Coverslip settings
cThickness = Rd;
nc = 1.52; % Coverslip refractive index
na = 1; % Outside coverslip refractive index (air)
cTop = [-Rd * 1.01, 0, 1]; %coverslip data (planes that define top and bottom)
cBottom = [-Rd * 1.01 - cThickness, 0, 1];

% Screen settings
screenDistance = -10;
screenSize = 20;
screen = [-screenDistance, 0, 1];
propagationDistance = 1000;

% Intensity mapping settings
imgDim = [500, 500];
screenHeight = 3;
startPlane = [-3, 1, 0];
distance = 10;
resolution = 200;
numSlices = resolution * distance;
minimumAmplitude = 0.01;

% Water plane
waterPlane = [8 * Rd, 0, 1];

% Droplet appearance
colorO = [0.69, 0.94, 1];
colorI = [1, 0.69, 0.69];
edgeC = [0, 0, 0];

%% Angle loop
angles = 0:5:90;
for angle = angles
    theta = deg2rad(angle);
    rDir = repmat([-cos(theta), sin(theta)], nRays, 1);
    
    % Set starting positions (uv)
    uv = zeros(nRays, 2);
    if angle == 0
        uv(:, 1) = linspace(-5 * Rd, 5 * Rd, nRays);
        uv(:, 2) = linspace(-5 * Rd, -5 * Rd, nRays);
    else
        uv(:, 1) = linspace(-5 * Rd * cos(theta), 5 * Rd, nRays);
        uv(:, 2) = linspace(-5 * Rd, -5 * Rd * sin(theta), nRays);
    end
    
    % Generate light direction name
    lightdir = sprintf('angle_%d_degrees_bottomleft', angle);
    
    % Display for debugging
    fprintf('Angle: %d degrees\n', angle);
    disp('Starting positions (uv):');
    disp(uv);
    disp('Direction vectors (rDir):');
    disp(rDir);
    
    % Ray tracing logic here
    %end
    
    % Contact angle generation
    [collectedRi, contactAngles] = generate(CA_start, CA_end, numPoints, vr, Rd);
    
    % Storage initialization
    histStorage = zeros(numPoints, 360);
    intHistStorage = zeros(numPoints, 360);
    histsqStorage = zeros(numPoints, 360);
    intHistsqStorage = zeros(numPoints, 360);
    hist = zeros(1,360);
    intHist = zeros(1,360);
    histsq = zeros(1,360);
    intHistsq = zeros(1,360);
    cint = zeros(imgDim(2), imgDim(1));
    emInt = zeros(imgDim(2), imgDim(1));
    cIntHist = zeros(screenHeight * 2 * resolution, numSlices);
    eint = zeros(imgDim(2),imgDim(1));
    environment = [nm, na, nc, ri, minimumAmplitude, propagationDistance, screenHeight, resolution, imgDim, no];
    name = strcat(lightdir,'-',num2str(contactAngles(1)),'-',num2str(contactAngles(length(contactAngles))),'-',num2str(nRays),'rd.mat');
    
    % Main processing loop
    tic;
    for ii = 1:numPoints
        hist = zeros(1, 360);
        intHist = zeros(1, 360);
        fint = zeros(imgDim(2), imgDim(1));
        emInt = zeros(imgDim(2), imgDim(1));
        
        if ismember(contactAngles(ii), interestAngles)
            logic(4) = true;
        else
            logic(4) = false;
        end
        
        Ri = collectedRi(ii);
        theDrop = drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
        rayNum = [];
        locsCut = uv;
        rayChop = nRays / rayPerCore;
        
        for l = 1:rayChop
            if l == rayChop
                choppedRays = locsCut((l - 1) * rayPerCore + 1:end, :);
                choppedDirs = rDir((l - 1) * rayPerCore + 1:end, :);
                rayNum = length(choppedDirs);
            else
                choppedRays = locsCut((l - 1) * rayPerCore + 1:l * rayPerCore, :);
                choppedDirs = rDir((l - 1) * rayPerCore + 1:l * rayPerCore, :);
                rayNum = length(choppedDirs);
            end
            
            cintTemp = emitterGeneration(length(choppedDirs), choppedRays, theDrop, choppedDirs, environment);
            emInt = emInt + cintTemp;
            
            [chist, cIntHist, cint] = multiDropParTrace(length(choppedDirs), choppedRays, theDrop, choppedDirs, logic, environment, [cTop; cBottom; waterPlane; startPlane; screen]);
            hist = hist + chist;
            intHist = intHist + cIntHist;
            fint = fint + cint;
        end
        
        dirStorage{ii} = rDir;
        histStorage(ii, :) = hist;
        emIntensityStorage{ii} = emInt;
        fIntensityStorage{ii} = fint;
        intHistStorage(ii, :) = intHist;
        
        save(name, 'contactAngles', 'collectedRi', 'dropLocation', 'Rd', 'Ri', 'no', 'ni', 'vr', 'dropAngle', 'histStorage', 'numPoints', 'environment', 'logic', 'nRays', 'imgDim', 'intHistStorage', 'emIntensityStorage', 'fIntensityStorage');
        %save(sprintf('angle_%d_direction_%s_nRays_%d_contactAngle_%d.mat', angle, 'collimated', nRays, contactAngles(ii)), 'contactAngles', 'collectedRi', 'dropLocation', 'Rd', 'Ri', 'no', 'ni', 'vr', 'dropAngle', 'histStorage', 'numPoints', 'environment', 'logic', 'nRays', 'imgDim', 'intHistStorage', 'emIntensityStorage', 'fIntensityStorage');
        fprintf('Processed point %d/%d\n', ii, numPoints);
    end
end
toc;
