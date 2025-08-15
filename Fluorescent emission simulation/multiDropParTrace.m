function [x,y,z]=multiDropParTrace(rayNum,choppedRays,theDrop,choppedDirs,logic,environment,planes)
%The purpose of this is for distributed paralellization of this
%function

%Fortran style because im an absolute 
useScreen = logic(1);
useCoverslip = logic(2);
useWaterPlane = logic(3);
forIntensity = logic(4);
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
rRays = choppedRays;
intRays = [];
rayStorage = [];
Rays = [];
dropTrace = [];
pointBin = [];
histBinSize = 360; %1 deg hist
hist = zeros(1,360);
intHist = zeros(1,360);
histsq = zeros(1,360);
intHistsq = zeros(1,360);
intMap = zeros(imgDim(2),imgDim(1));
int = zeros(imgDim(2),imgDim(1));

for j = 1:rayNum
    aRay = ray([rRays(j,1), rRays(j,2)], choppedDirs(j, :), ri);
    %this is the line that actually does the raytracing calculation:
    aRayTraced = traceRay(aRay, theDrop, nm, minimumAmplitude,1000);
    dropTrace = aRayTraced;

    %This one is literally just to look at the rays as they come off
    %the droplet vs the rays after the water plane
    %now histBin is an array containing all of the points where rays hit, now
    %we place it into a histogram for each ray
    rayLength = length(dropTrace);
    for k = 1:rayLength
        rayData = [dropTrace(1,k).location(1) dropTrace(1,k).location(2) dropTrace(1,k).direction(1) dropTrace(1,k).direction(2) dropTrace(1,k).amplitude];
        
        X = real(dropTrace(1,k).direction(1));
        Y = real(dropTrace(1,k).direction(2));
        if isreal(X) && isreal(Y)
            if forIntensity == true
                int = IntensityPlaneAnyCustom(rayData,imgDim,resolution);
            end
            intMap = intMap + int;
            angle = (atan2d(Y,X));
            if angle < 0
                angle = 360-abs(angle);
            end
            point = round(angle);
            if point <= 359
                point = point+1; %an array goes from 1to360, or 0to359 later
            else
                point = 1; %Because the array is actually 0 to 359, we place this at the one position
            end
            intHist(:,point) = intHist(:,point) + abs(dropTrace(1,k).amplitude);
            point = [];
        else
            
        end
    end
    
    %The fifth entry is default 1000, and is the max interactions
    %propagate to screen and store
    if isempty(aRayTraced) == 0
        for l = 1:length(dropTrace)
            aRayTraced = dropTrace(l);
            if useCoverslip
                [loc, norm] = rayPlaneIntersection(aRayTraced, cTop);
                if not(isnan(loc))
                    aRayTraced = propagate(aRayTraced, loc);
                    [aRayTraced, tRef] = refract (aRayTraced, norm, nc) ;
                end
                [loc, norm] = rayPlaneIntersection(aRayTraced, cBottom);
                if not(isnan(loc))
                    aRayTraced = propagate(aRayTraced, loc);
                    [aRayTraced, tRef] = refract (aRayTraced, norm, na) ;
                end
            end
            if useWaterPlane
                [loc, norm] = rayPlaneIntersection(aRayTraced, waterPlane);
                if not(isnan(loc))
                    aRayTraced = propagate(aRayTraced, loc);
                    [aRayTraced, tRef] = refract (aRayTraced, norm, na) ;
                end
            end
            if useScreen
                [loc, norm] = rayPlaneIntersection(aRayTraced, screen);
                if (abs(loc(1))< screenSize)
                    aRayTraced = propagate(aRayTraced, loc);
                end
            else
                aRayTraced = propagate(aRayTraced, propagationDistance);
            end
            %Histogram/intensity generation for water plane
            rayLength = length(aRayTraced);
            for k = 1:rayLength
                X = real(aRayTraced(1,k).direction(1));
                Y = real(aRayTraced(1,k).direction(2));
                if isreal(X) && isreal(Y) && aRayTraced(1,k).currentIndex == 1 
                    angle = (atan2d(Y,X));
                    if angle < 0
                        angle = 360-abs(angle);
                    end
                    point = round(angle);
                    if point <= 359
                        point = point+1; %an array goes from 1to360, or 0to359 later
                    else
                        point = 1; %Because the array is actually 0 to 359, we place this at the one position
                    end
                    hist(:,point) = hist(:,point) + abs(aRayTraced(1,k).amplitude);
                else
                    
                end
            end
        end
    end
    aRayTraced = [];
    dropTrace = [];
end
x = hist;
y = intHist;

z = intMap;
end