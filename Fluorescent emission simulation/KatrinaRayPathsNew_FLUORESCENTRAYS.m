clc;clear all;close all;

% Make sure DropletTracer folder is in your path. 
addpath('DropletTracer')

%% Droplet parameters:

% purposes
options = optimset('Display','off');
tic;
colorO = [0.69, 0.94, 1];
colorI = [1, 0.69, 0.69];
edgeC=[0, 0, 0];

vr=1; %volume ratio.  V_internal/v_external
Rd=1; %droplet radius
dropLocation=[0,0];
dropAngle=0;  % rotation angle (degrees)

angleRadii = Rd;

% Droplet and Phase Info

% Indexes of refraction: 
nm=1.334; % Surrounding Medium
no=1.361; % Outer Phase
ni=1.462;% Inner Phase

%% Raytracing Parameters:

%set number of rays:
rayDensity = 10;
nRays = 100;
% How many droplets
    numPoints = 3; %In this current version, this is the amount of drops between 0 and 90, which gets doubled.
    contactAngles = [15 30 45]; %linspace(1,90,numPoints); %This must ALWAYS run with a 90 degree drop, as the baseline
    interestAngles = []; %[0,10, 27, 40, 60, 90];%[0 12 24 42 60 90];% [0 10 26 45 60 90];% ;%contactAngles; %, 10, 26, 50, 90, 140, 164,170,180]; %These are the angles we gather intensity maps for
    reverseAngles=true; %Make true if you want to look at the other side of the droplet conformations
    
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

Ri = collectedRi(2);
contactAngle=contactAnglesCorrected(2);

%% 
%include coverslip (or any glass substrate)
useCoverslip=false;
if useCoverslip
   cThickness=105 %thickness of coverslip
   nc=1.52 %refractive index of coverslip
   na=1 %refractive index outside of coverslip (air probably)

   % coverslip data (planes that define top and bottom): 
   cTop=[-Rd, 0, 1];
   cBottom=[-Rd-cThickness, 0, 1];
end

%if you use a thicker substrate, check the location of the screen so that it is after the bottom edge of the coverslip
screenDistance=(Rd+20);
%plane to propagate rays to after interaction with the droplet
%if rays don't hit this plane they are thrown away:
%equation of plane is p(1)=p(2)x+p(3)y. 
screen=[+screenDistance, 0, 1]; % Plane that will be the screen.  
screenSize=200

%you can choose to use a screen
useScreen=false; 
propagationDistance=100; %instead propagate all rays after they leave the drop by this amount

%%
%set ray start points: 
useBradGeneration = true;

if useBradGeneration == false
    xStart=linspace(-2.5, 0, 2*Rd*rayDensity); %full width of droplet,
    numRays=length(xStart);
    %numRays = nRays;
    %yStart=(Rd-5)*ones(numRays); %5 units above droplet.
    yStart = linspace(0, 2.5, 2*Rd*rayDensity);
    %set ray start directions:
    rDir=repmat([1, -1],numRays ,1 );
end

%Ray amplitude is between 0 and 1, through away rays that have less than:
minimumAmplitude=0.1;

%% Draw all the drops&rays

for i = 1:length(collectedRi)
    theDrop=drop(dropLocation, Rd, collectedRi(i), no, ni, vr, dropAngle); %invent your droplet
    
    if useBradGeneration == true
        tempFig = figure('visible','off'); %We take the droplet and get the polygon out of it
        [fO, fI]=drawDrop(theDrop);
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
        xStart = uv(:,1);
        yStart = uv(:,2);
    end
    
    figure;

    drawDrop(theDrop,100,edgeC,colorO,colorI)
    
    axis image
    xlabel('x')
    ylabel('y')
    xlim([-5 5])
    ylim([-5 10])
    set(gcf, 'color', 'white')
    
% Store Rays in:

    RaysAll=[];
    
    for r=1:numRays
        rRay=ray([xStart(r), yStart(r)], rDir(r, :), nm);
        rRayTraced=traceRay(rRay, theDrop, nm, minimumAmplitude);
        
        for p=rRayTraced
            
            %coverslip: 
            if useCoverslip
                [loc, norm]=rayPlaneIntersection(p, cTop);
                if not(isnan(loc))
                    p=propagate(p, loc);
                    [p, tRef] = refract (p, norm, nc) ;
                end
                [loc, norm]=rayPlaneIntersection(p, cBottom);
                if not(isnan(loc))
                    p=propagate(p, loc);
                    [p, tRef] = refract (p, norm, na) ;
                end
            end
            
            %screen:
            if useScreen
               [loc, norm]=rayPlaneIntersection(p, screen);
            if (abs(loc(1))< screenSize)
                p=propagate(p, loc);
                RaysAll=[RaysAll, p];
            end
            
            else
            p=propagate(p, propagationDistance);
            RaysAll=[RaysAll, p];
            end
        end
    end
    
   % Draw the rays.  (Modify to draw only selected rays.)
    for Ray=RaysAll
        drawRay(Ray, [0, 0, 0.5]);
        hold on        
    end
    
% Run calculation
 
    for n=1:numDrops
    
    % Make the drop object:
    theDrop=drop(dropLocation, Rd, Ri, no, ni, vr, dropAngle);

    % An Array to store the Rays in: 
    Rays=[];

    % Send Rays through scene
    for ii=1:numRays
        %create ray object: 
        aRay=ray([xStart(ii), yStart(ii)], rDir(ii, :), nm);
        
        %this is the line that actually does the raytracing calculation:
        aRayTraced=traceRay(aRay, theDrop, nm, minimumAmplitude);
        
        for t=aRayTraced
            
            %coverslip: 
            if useCoverslip
                [loc, norm]=rayPlaneIntersection(t, cTop);
                if not(isnan(loc))
                    t=propagate(t, loc);
                    [t, tRef] = refract (t, norm, nc) ;
                end
                [loc, norm]=rayPlaneIntersection(t, cBottom);
                if not(isnan(loc))
                    t=propagate(t, loc);
                    [t, tRef] = refract (t, norm, na) ;
                end
            end
            
            %screen:
            if useScreen
               [loc, norm]=rayPlaneIntersection(t, screen);
            if (abs(loc(1))< screenSize)
               t=propagate(t, loc);
               Rays=[Rays, t];
            end
            
            else
            t=propagate(t, propagationDistance);
            Rays=[Rays, t];
            end       
        end       
    end   

% plot the output: 

% Draw the rays. (Modify to draw only selected rays.)

        for Ray=Rays
            drawRay(Ray, [0, 0, 0.5]);
            hold on 
        end   
    end
    
    if useCoverslip
       plot([-2*d, 2*d], [-d, -d], 'c')
       plot([-2*d, 2*d], [-d-cThickness, -d-cThickness], 'c') 
    end    
end

%%

