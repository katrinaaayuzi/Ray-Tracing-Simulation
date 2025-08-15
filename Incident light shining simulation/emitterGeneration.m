function [x] = emitterGeneration(rayNum,choppedRays,theDrop,choppedDirs,environment)
%This function is used in order to trace rays over a droplet, and get a
%heat map of a particular phase which can be used for the generation of
%emitters in further ray tracing

%Fortran style because im an absolute
nm = environment(1);
na = environment(2);
nc = environment(3);
ri = environment(4);
minimumAmplitude = environment(5);
propagationDistance = environment(6);
resolution = environment(8);
imgDim = [environment(9) environment(10)];
rRays = choppedRays;
intMap = zeros(imgDim(2),imgDim(1));
int = zeros(imgDim(2),imgDim(1));

for j = 1:rayNum
    aRay = ray([rRays(j,1), rRays(j,2)], choppedDirs(j, :), ri);
    %this is the line that actually does the raytracing calculation:
    aRayTraced = traceRayReturnAll(aRay, theDrop, nm, minimumAmplitude,500);
    for ph = 1:length(aRayTraced)
        t = aRayTraced(ph);
        t = propagate(t, propagationDistance);
        aRayTraced(ph) = t;
    end
    
    rayData = [];
    for pi = 1:length(aRayTraced)
        currentRay = [];
        currentRay = aRayTraced(pi);
		rayLength = size(currentRay.trajectory,1);
		if rayLength > 2
			for pj = 1:(rayLength-1)
            %rayData is a n by 6 array where (x loc) (y loc) (x loc2)
            %(y loc2) (amplitude) (index)
            nextNum = pj+1;
            rayData = [rayData;currentRay.trajectory(pj,1),currentRay.trajectory(pj,2),currentRay.trajectory(nextNum,1),currentRay.trajectory(nextNum,2),currentRay.amplitudeHistory(pj),currentRay.indexHistory(pj)];
			end
        end
    end
	
    if isempty(rayData)
    else
        testIndex = max([theDrop.innerIndex theDrop.outerIndex]); %always grabs the larger refractive index because we <3 hydrocarbons
        rDel = find(rayData(:,6) ~= testIndex);
        rayData(rDel,:)=[];
        int = int + IntensityPlaneAnyDrop(unique(rayData,'rows'),imgDim,resolution); %TODO
    end
    
end

x = int;
end