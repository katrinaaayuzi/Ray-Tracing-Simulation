clear all; %close all; clc;
% After Response from Sara

    vr = 1;
    Rd = 1;
    n_w = 1.334;
    n_hc = 1.462;
    n_fc = 1.361;
    dropAngle = 0;

    CA_start = 0;
    CA_end = 180;
    numPoints = 31; %In this current version, this is the amount of drops between 0 and 90, which gets doubled.

    [collectedRi, contactAngles] = generate(CA_start, CA_end, numPoints, vr, Rd);

for j = 1:numPoints
    Ri = collectedRi(j);
    contactAngle = deg2rad(contactAngles(j));
    loc = [1 1];
    theDrop = drop(loc, Rd, Ri, n_w, n_hc, vr, dropAngle);
    tempFig = figure('visible','off'); %We take the droplet and get the polygon out of it
    [fO, fI] = drawDrop(theDrop);
    close(tempFig);
    polyGenIU = polyshape({fO(1,:)},{fO(3,:)});
    polyGenII = polyshape({fO(2,:)},{fO(4,:)});
    drops{1} = [polyGenII.Vertices(1,:);polyGenII.Vertices(end/2,:);polyGenII.Vertices(end,:)];
    drops{2} = [polyGenIU.Vertices(1,:);polyGenIU.Vertices(end/2,:);polyGenIU.Vertices(end,:)];
    
    for i = 1:2
        drop_coords = [drops{i}];
        x = drop_coords(:,1);
        y = drop_coords(:,2);
        leftOrthog = (y(2)-y(1))./(x(2)-x(1));
        rightOrthog = (y(3)-y(2))./(x(3)-x(2));
        interfaceFit = (leftOrthog.*rightOrthog.*(y(3)-y(1))+leftOrthog.*(x(2)+x(3))-rightOrthog.*(x(1)+x(2)))./(2*(leftOrthog-rightOrthog));
        interfaceFit(:,2) = -1./leftOrthog.*(interfaceFit-(x(1)+x(2))/2)+((y(1)+y(2))/2);
        interfaceRadii = sqrt((interfaceFit(:,1)-x(1)).^2+(interfaceFit(:,2)-y(1)).^2);
        if isempty(interfaceFit) == 0
            curveCenters(i,:) = interfaceFit(1,:);
            curveRadii(i) = interfaceRadii(1,:);
        end
    end
    
    loc_outer = curveCenters(1,:);
    radii_outer = curveRadii(1);
    loc_inner = curveCenters(2,:);
    radii_inner = curveRadii(2);
    d = norm(loc_outer - loc_inner); %distance between circle centerpoints
    contact_angle = pi - acos(((radii_inner^2)+(radii_outer^2)-(d^2))/(2*radii_inner*radii_outer)); %sanity check
    %contact_angle = acos(((radii_inner^2)+(radii_outer^2)-(d^2))/(2*radii_inner*radii_outer)); %sanity check
    contact_angle_new(j) = rad2deg(contact_angle); %perhaps the error in generating a droplet and translating back gave the 1-3deg CA error in final data
    
    %theta_c is the solid angle, theta_p is the direction of the three-phase
    %interface
    theta_critical_hcw = asind(n_w/n_hc);%critical angle
    theta_critical_hcfc = asind(n_fc/n_hc);
    
    theta_p = asind((radii_outer^2 - radii_inner^2 - d^2)/(2*radii_inner*d)); %Direction the three-phase interface points in global coordinate system
    theta_p_out(j) = theta_p;

    theta_t = 90 - contact_angle_new(j);
    theta_t_prime(j) = asind(sind(theta_t)*(n_hc/n_w));
    theta_n(j) = -(acosd((radii_inner^2 - radii_outer^2 - d^2)/(-2*radii_outer*d))); %'is normal to surface of droplet relative to gravity'
    theta_o(j) = theta_n(j) + theta_t_prime(j);
    

    theta_w = (theta_critical_hcfc - 90);
    %theta_w = (90 - theta_critical_hcfc);
    theta_o_one(j) = theta_o(j) + theta_w;
    theta_o_two(j) = theta_o(j) - theta_w;
    
    %gamma_i = contact_angle_new(j) - 90;
    gamma_i = 90 - contact_angle_new(j);
    gamma_i_a(j) = gamma_i + (90-theta_critical_hcfc);
    gamma_i_b(j) = gamma_i - (90-theta_critical_hcfc);
    %For gamma we don't need to set a normal, as gamma is already relative
    %to the tangent at the triple-phase-interface, gamma_i_one or i_two
    %only needs to fulful the condition, is it larger than pi/2 -
    %theta_critical_hcw
    
    %here you can check the limit programmatically, however it takes care
    %of itself more or less
    if gamma_i_a(j) >= (theta_critical_hcw)
        gamma_i_a(j) = theta_critical_hcw;
    else
        %j
    end
    
    theta_a(j) = theta_n(j) + asind(sind(gamma_i_a(j))*(n_hc/n_w)); %not sure that + or - here is correct
    theta_b(j) = theta_n(j) + asind(sind(gamma_i_b(j))*(n_hc/n_w));
    theta_center(j) =  (theta_a(j) + theta_b(j))/2;
    theta_width(j) = theta_a(j)+90 - (theta_b(j)+90);

end

%
%figure; hold on;  box on; %grid on;
%plot(contact_angle_new,theta_width,'LineWidth',2,'DisplayName','\theta_{center}')

  
figure; hold on;  box on; %grid on;
%plot(contact_angle_new,theta_p_out,'LineWidth',2,'DisplayName','\theta_p')
%plot(contact_angle_new,theta_center,'LineWidth',2,'DisplayName','\theta_{center}')
% plot(contact_angle_new,theta_a,'LineWidth',2,'DisplayName','\theta_{a}')
% plot(contact_angle_new,theta_b,'LineWidth',2,'DisplayName','\theta_{b}')
%plot(contact_angle_new,(theta_n),'LineWidth',2,'DisplayName','\theta_n')

%plot(contact_angle_new,(gamma_i_a),'LineWidth',2,'DisplayName','\gamma_a')
%plot(contact_angle_new,(gamma_i_b),'LineWidth',2,'DisplayName','\gamma_b')
%plot(contact_angle_new,(theta_o_one),'LineWidth',2,'DisplayName','\theta_o + \theta_w')
%plot(theta_p_ref_limOne_out(:,1),(theta_p_ref_limOne_out(:,2)),'LineWidth',2,'\theta_p_r_e_f')
%plot(theta_p_ref_limTwo_out(:,1),(theta_p_ref_limTwo_out(:,2)),'LineWidth',2,'\theta_p_r_e_f')

plot(theta_a,contact_angle_new,'LineWidth',2,'DisplayName','\theta_{a}')
plot(theta_b,contact_angle_new,'LineWidth',2,'DisplayName','\theta_{b}')

%ylim([-90 90]); xlim([0 90]);
% ylabel('angle of emission')
% xlabel('droplet contact angle')
xlim([-90 90]);
ylim([0 90]);
xlabel('Angle of emission')
ylabel('Droplet contact angle')
legend
a = real([contact_angle_new',theta_a',theta_b']);

%% lukas
figure; hold on;  box on; %grid on;
%plot(contact_angle_new,(theta_p_out),'LineWidth',2,'DisplayName','\theta_p')
plot(contact_angle_new,(theta_t_prime),'LineWidth',2,'DisplayName','\theta_t_,_p_r_i_m_e')
plot(contact_angle_new,(theta_o),'LineWidth',2,'DisplayName','\theta_o')
%plot(contact_angle_new,(theta_o_one),'LineWidth',2,'DisplayName','\theta_o + \theta_w')
%plot(contact_angle_new,(theta_o_two),'LineWidth',2,'DisplayName','\theta_o - \theta_w')


%plot(contact_angle_new,rad2deg(test),'LineWidth',2)

%plot(theta_p_ref_limOne_out(:,1),(theta_p_ref_limOne_out(:,2)),'LineWidth',2)
%plot(theta_p_ref_limTwo_out(:,1),(theta_p_ref_limTwo_out(:,2)),'LineWidth',2)
%plot(contact_angle_new,rad2deg(theta_p_ref_crit))
%plot(contact_angle_new,rad2deg(pi/2 -theta_n_out))

ylim([-90 90]); xlim([0 90]);
ylabel('Angle of emission, (^{\circ})')
xlabel('Droplet contact angle (\theta, ^{\circ})')
legend


%%
%plot((theta_p_ref_out(:,2))+90,theta_p_ref_out(:,1),'LineWidth',2)
%plot((theta_p_ref_limOne_out(:,2))+90,theta_p_ref_limOne_out(:,1),'LineWidth',2)
%plot((theta_p_ref_limTwo_out(:,2))+90,theta_p_ref_limTwo_out(:,1),'LineWidth',2)
%plot(theta_o_one+90,contact_angle_new+1,'LineWidth',2,'DisplayName','theta_o_l')
figure; hold on; ylim([0 90]); xlim([0 178]);
plot(theta_center+91,contact_angle_new,'LineWidth',2,'DisplayName','\theta_{o,center,refraction}')
plot(theta_a+91,contact_angle_new,'LineWidth',2,'DisplayName','\theta_{o,withlimit,refraction}')
plot(theta_b+91,contact_angle_new,'LineWidth',2,'DisplayName','\theta_{o,nolimit,ref,refraction}')
legend
%scatter(rad2deg(theta_p_ref_crit)+90,contact_angle_new,'LineWidth',2)


%scatter(rad2deg(theta_o_max)+90,contact_angle_new,'LineWidth',2)
%plot(rad2deg(theta_n_out),contactAngles,'LineWidth',2)



%scatter(rad2deg(theta_p_out)+90-width/2,contact_angle_new,'LineWidth',2)
%plot(rad2deg(theta_o)+90,contact_angle_new,'LineWidth',2)
% plot(rad2deg(refraction_out),contactAngles,'LineWidth',2)
% plot(rad2deg(theta_o),contactAngles,'LineWidth',2)

%plot(rad2deg(new_model),contactAngles,'LineWidth',2)


%%
output = [];
output = [contact_angle_new',rad2deg(theta_p_ref_out'),rad2deg(theta_p_ref_limOne_out'),rad2deg(theta_p_ref_limTwo_out')];
