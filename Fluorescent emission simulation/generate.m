function [collectedRi, collectedCA] = generate(CA_start, CA_end, numPoints, vr, Rd)
%Simple and fun little function to generate droplets and clean things up.
%This operates with drop.m as well as dropShapeSolverContactAngle which
%each have a change, so if this file is used, make sure you have the proper
%versions of those files so the droplets will properly generate.
contactAngles = linspace(CA_start, CA_end, numPoints);
collectedRi = zeros(1,numPoints);

for j = 1:numPoints
    contactAngle = cos(deg2rad(contactAngles(j)));
    if contactAngle >= 0 %an attempt to try and get the correct values out without specific hacks
        contactAngle = contactAngle;
        d_mod = 1;
        vrTemp = vr;
    elseif contactAngle < 0
        contactAngle = -contactAngle;
        vrTemp = 1/vr;
        d_mod = -1;
    end
    [solu,fval,exitflag,output]=fsolve(@(x)dropShapeSolverContactAngle(x, Rd, vrTemp, contactAngle), [0.8, 1] );
    if (exitflag<0)
        continue
    end
    d = solu(1);
    Ri(j) = solu(2)*d_mod;
end
collectedRi = Ri;
collectedCA = contactAngles;
end

