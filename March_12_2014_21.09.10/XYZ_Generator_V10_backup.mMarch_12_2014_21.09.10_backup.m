function [Section,Re] = XYZ_Generator_V10()
% Here are the inputs of the script:
%     l  - inches . The incremental length of the turbine blade ( l <=
%                   totallength)
%     totalLength - (inches) total length of the blade ~ 4 inches
%     4415 - NACA airfoil code - I decided to use 4415

% Outputs:
%     xyz#.xls  - (x,y,z) points defining a curve. Import these into
%               SolidWorks and use the loft function to create the blade.
%               there are about 20 of these files, you do not have import
%               all of them as long as you import profileLine.sldcrv.
%     profileLine.sldcrv - (x,y,z) points defining a the guideline for the
%               blade. Import it into solidworks. Loft will ask for a
%               guide curve, use this curve to create a smooth blade.

close all;clear all
hub_r = 0.125; %radius of the hub to the beginning of the blade
totalLength=3.5; % inches total length
chord_max=0.75; % inches - the maximum chord width
Cl = 0.8; %optimum coefficient of lift
Om_des = 12000; %Desired rotational speed in RPM

%Blade Planform Profile
%the interpolation technique to derive the planform profile
method = 'pchip';

%Non-uniform wind distribution information
V_c = 40; %Wind speed at center [mph]
V_e = 10; %Wind speed at edge of the flow [mph]
V_r = 4; %radial distance to end wind speed measurement [in]

%Set the number of profiles along the length of the blade
detail = 15;

% set the folder to write the turbine profile data.
dt = datestr(now,'mmmm_dd_yyyy_HH.MM.SS')
mkdir(dt)
oldFolder = cd(dt);

%This will create an empty Section matrix that stores all the values
%row 1 - radial distance
%row 2 - chord length
%row 3 - Optimum chord length by Betz Method
%row 4 - x_move
%row 5 - rotation
%row 6 - apparent wind speed seen by section


Section = zeros(5,detail);
Re = zeros(1,detail);
V_e = 17.6 * V_e; %Convert [mph] to [in/s]
V_c = 17.6 * V_c; %Convert [mph] to [in/s]
V_true = [V_c V_e];
V_r = [0 V_r];

%Convert Om_des to rads/sec
%Multiple Om_des by 2 to go from computational to experimental
Om_des = Om_des * ((2 * pi())/(60));

%This for will write the profiles for the turbine blade
n=0;
i = 1;
for l=0:(totalLength/(detail-1)):(totalLength)
    
    % defining variables for naca4gen.m
    if (l/totalLength) < 0.01
        iaf.designation='4465';
    else
        iaf.designation='4415';
    end
    
    % designation='0008';
    iaf.n=30;
    iaf.HalfCosineSpacing=1;
    iaf.wantFile=1;
    iaf.datFilePath='./'; % Current folder
    iaf.is_finiteTE=0;
    blade = naca4gen(iaf);
    xyz=[blade.xU blade.zU;blade.xL(2:end) blade.zL(2:end)];
    lxyz=length(xyz);
    [IDX, C0]=kmeans(xyz,1); % this calculates the centroid of the airfoil
    
    %These populate the Section matrix
    Section(1,i) = l;
    %This calculates the chord thickness for the current profile
    Section(2,i) = Chord_Thick(l,totalLength,chord_max,method);
    %This calculates how far to move the centroid to line it up
    Section(4,i) = abs((chord_max - Section(2,i))/2);
    
    %This calculates the apparent wind speed on each profile
    %This interpolates the incoming wind speed from the given measurements
    u = interp1(V_r,V_true,l,'linear'); %incoming wind speed [in/s]
    %Use desired Omega, Om_des rather than calculating Omega
    w = (Section(1,i) + hub_r)*Om_des;
    W = sqrt(u^2 + w.^2);
    %apparent wind speed in mph
    Section(6,i) = W/(17.6);
    
    %This calculates how much to twist it using Betz approx.
    Section(5,i) = Twist_Exp(l,hub_r,Om_des,u) * (180/pi());
    
    %Calculate Reynolds number
    kin_vi = 236.16E-4; %kinematic viscosity, in^2/s
    Re(1,i) = (W * Section(2,i))/(kin_vi);
    
    X_calc = (Om_des * (totalLength + hub_r)) / (V_c);
    %Calculate the optimum chord length by Betz
    C_opt = ((2*pi()*(totalLength+hub_r)*8*V_c)/(3*9*Cl*X_calc*W));
    
    Section(3,i) = C_opt;
    
    %Calculations
    tempXYZ=Section(2,i)*xyz;
    
    % this is for aligning the centroids of the airfoils
    %Ct=centroid(tempXYZ); % Download it from Mathworks website

    %This section of code will keep the leading edge straight
    %xShift=[(-1) * Section(3,i) 0];

    %This section of code rotates the profile
    twistAng=Section(5,i)*(pi()/180);
    M=[cos(twistAng) -sin(twistAng); % rotation matrix
        sin(twistAng) cos(twistAng)];

    % rotate by twistAng using rotation matrix
    tempXYZ=M*(tempXYZ');

    %Begin blade profile file generation
    n=n+1;
    %     create a new filename
    nstr=num2str(n); % changing the n to string
    filename=strcat('Profile_',nstr,'.sldcrv');
    
    profile(n,:,:)=[tempXYZ' l*ones(lxyz,1)];
    profileLine(n,:)=[profile(n,1,1) profile(n,1,2) l];
    tempM(:,:)=profile(n,:,:); % this is after realization that dl
    %     dlmwrite does not work with 3D matrice
    dlmwrite(filename,tempM,'delimiter','\t','newline','pc'); % write a file
    %      import these files into Solidworks as curves
    i = i + 1;
end

dlmwrite('profileLine.sldcrv',profileLine,'delimiter','\t','newline','pc');

plot(Section(1,:),Section(2,:),Section(1,:),Section(3,:))
xlabel('Radial Length [in]');
ylabel('Chord Length [in]');
legend('Existing Chord Length', 'Optimum Chord Length')
pause
plot(Section(1,:),Section(5,:))
xlabel('Radial Length [in]');
ylabel('Angle of Twist [deg]');
pause
plot(Section(1,:),Section(6,:))
xlabel('Radial Length [in]');
ylabel('Apparent Wind Speed [mph]');

Re

%This code will write the current version of the code to the
%Same folder as the profile files
Version=dt;
FileNameAndLocation=[mfilename('fullpath')];
New_FileNameAndLocation=[mfilename(dt)];
newbackup=sprintf('%s_backup.m',New_FileNameAndLocation,Version);
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

cd(oldFolder);

end

function [Chord_size] = Chord_Thick(l,totalLength,chord_max,method)

%Parameters denoting the information for the root of the blade
Root_P = 0.0625; %Percent where the root ends
Root_T_S = 0.33; %Percent of max. chord width at beginning of root
Root_T_E = 0.6; %Percent of max. chord width at end of root

%Parameters denoting the information for the expanding of the blade
Expand_P = 0.13; %Percent where the expanding ends
Expand_T_E = 1.0; %Percent of max. chord width at end of root

%Parameters denoting the information for the shrinking of the blade
Shrink_P = 1; %Percent where the shrinking ends
Shrink_T_E = 0.05; %Percent of max. chord width at end of root

x = [0 Root_P Expand_P Shrink_P];
v = [Root_T_S Root_T_E Expand_T_E Shrink_T_E];

xq = (l/totalLength);
Chord_size = chord_max * interp1(x,v,xq,method);

end

function [theta] = Twist_Exp(l,hub_r,Om_des,u)

%---------------------------Input Data-----------------------------
%R = hub_r + totalLength; %Blade Radius, inches
%---------------Calculations for whole blade--------------------
%Omega = X*(V_0/(hub_r + totalLength)); %Rotational speed, s^-1.
r = (l + hub_r); %radius where the angle is calculated, [in]
Alpha = -20 * ((2*pi())/360); %Angle of attack of 10 deg, converted to rads

%--------------Determine Factors at Section-----------------------
alpha = Alpha;

%Use desired Om_des rather than calculating Omega
w = r*Om_des;

phi = atan(w/u);
theta = alpha - phi;

end