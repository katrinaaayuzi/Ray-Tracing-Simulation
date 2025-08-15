clear all; clc;
options = optimset('Display', 'off');

% Droplet parameters
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
contactAngles = linspace(CA_start, CA_end, numPoints);
interestAngles = 90:5:180;
reverseAngles = true;

% Ray tracing parameters
nRays = 1000;
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
cTop = [-Rd * 1.01, 0, 1];
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

% Incident light angles
lightAngles = 0:5:90;

% Initialize storage for each light angle
lightDirData = struct();

% Main processing loop
tic;
for angle = lightAngles
    theta = deg2rad(angle);
    rDir = repmat([-cos(theta), sin(theta)], nRays, 1);

    % Generate light direction name
    lightdir = sprintf('angle_%d_degrees_bottomleft', angle);

    % Initialize storage for this light angle
    lightDirData(angle + 1).lightAngle = angle;
    lightDirData(angle + 1).lightDirName = lightdir;
    lightDirData(angle + 1).contactAngles = contactAngles;
    lightDirData(angle + 1).histStorage = zeros(numPoints, 360);
    lightDirData(angle + 1).emIntensityStorage = cell(1, numPoints);
    lightDirData(angle + 1).fIntensityStorage = cell(1, numPoints);
    lightDirData(angle + 1).intHistStorage = zeros(numPoints, 360);
    
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
    
    for ii = 1:numPoints
        % Set starting positions (uv)
        uv = zeros(nRays, 2);
        if angle == 0
            uv(:, 1) = linspace(-5 * Rd, 5 * Rd, nRays);
            uv(:, 2) = linspace(-5 * Rd, -5 * Rd, nRays);
        else
            uv(:, 1) = linspace(-5 * Rd * cos(theta), 5 * Rd, nRays);
            uv(:, 2) = linspace(-5 * Rd, -5 * Rd * sin(theta), nRays);
        end

        hist = zeros(1, 360);
        intHist = zeros(1, 360);
        fint = zeros(imgDim(2), imgDim(1));
        emInt = zeros(imgDim(2), imgDim(1));

        if ismember(contactAngles(ii), interestAngles)
            logic(4) = true;
        else
            logic(4) = false;
        end

        Ri = contactAngles(ii); % Updated to use current contact angle
        theDrop = drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);
        locsCut = uv;
        rayChop = nRays / rayPerCore;

        for l = 1:rayChop
            if l == rayChop
                choppedRays = locsCut((l - 1) * rayPerCore + 1:end, :);
                choppedDirs = rDir((l - 1) * rayPerCore + 1:end, :);
            else
                choppedRays = locsCut((l - 1) * rayPerCore + 1:l * rayPerCore, :);
                choppedDirs = rDir((l - 1) * rayPerCore + 1:l * rayPerCore, :);
            end

            cintTemp = emitterGeneration(length(choppedDirs), choppedRays, theDrop, choppedDirs, environment);
            emInt = emInt + cintTemp;

            [chist, cIntHist, cint] = multiDropParTrace(length(choppedDirs), choppedRays, theDrop, choppedDirs, logic, environment, [cTop; cBottom; waterPlane; startPlane; screen]);
            hist = hist + chist;
            intHist = intHist + cIntHist;
            fint = fint + cint;
        end

        % Aggregate results for the current contact angle
        lightDirData(angle + 1).histStorage(ii, :) = hist;
        lightDirData(angle + 1).emIntensityStorage{ii} = emInt;
        lightDirData(angle + 1).fIntensityStorage{ii} = fint;
        lightDirData(angle + 1).intHistStorage(ii, :) = intHist;
    end

    % Save all data for this light angle
    save(sprintf('LightAngle_%d_nRays_%d.mat', angle, nRays), 'lightDirData');
    fprintf('Processed Light Angle: %d, LightDir: %s\n', angle, lightdir);
end
toc;
