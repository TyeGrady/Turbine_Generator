function [] = XYZ_Generator_V3()
% Here are the inputs of the script:
%     l  - inches . The incremental length of the turbine blade ( l <=
%                   totallength)
%     totalLength - (inches) total length of the blade ~ 4 inches
%     tipThick   - (inches) the thickness of the blade tip
%     L0  - inches) the length were tapering starts
%     tipCurveL   - inches)  length were the blade-tip starts to curve
%     tipCurveR   -inches ) the radius of curvature of the blade-tip
%     4415 - NACA airfoil code - I decided to use 4415
%     cl - the percentage length where twisting stops
%     theta(l) - twist angle as a function of blade incremental length (I
%                 suggest you modify this function to create a smooth
%                 twist: fitting a 3rd order polynomial might do it)
%     twistAng - (rad) calculated by theta(l)
% Outputs:
%     xyz#.sldcrv  - (x,y,z) points defining a curve. Import these into
%               SolidWorks and use the loft function to create the blade.
%               there are about 20 of these files, you do not have import
%               all of them as long as you import profileLine.sldcrv.
%     profileLine.sldcrv - (x,y,z) points defining a the guideline for the
%               blade. Import it into solidworks. Loft will ask for a
%               guide curve, use this curve to create a smooth blade.

close all;clear all
hub_r = 0.125; %radius of the hub to the beginning of the blade
totalLength=3.75; % inches total length
chord_max=0.5; % inches - the maximum chord width

%Set the number of profiles
detail = 25;

l=0:0.05:totalLength;
n=0;
% Write the coordinates to a text file

% set the folder to write the turbine profile data.
cd('D:\My Documents\Class Files\SE 120\Propeller Design\MATLAB Code 2\Generated Files');
% Change this folder to a folder on your desktop. The script will write



%This for will write the profile for the turbine blade
i = 1;
for l=0:(totalLength/detail):(totalLength)
    i = i + 1;
    % defining variables for naca4gen.m
    iaf.designation='4415';
    % designation='0008';
    iaf.n=30;
    iaf.HalfCosineSpacing=1;
    iaf.wantFile=1;
    iaf.datFilePath='./'; % Current folder
    iaf.is_finiteTE=0;
    blade = naca4gen(iaf);
    s=5;
    l1=length(blade.xU);
    l2=length(blade.xL(2:end));
    xyz=[blade.xU blade.zU;blade.xL(2:end) blade.zL(2:end)];
    lxyz=length(xyz);
    [IDX, C0]=kmeans(xyz,1); % this calculates the centroid of the airfoil
    
    %New blade setting calculations based on  Wind Energy eqs 2.67-2.70
    %tip speed ratio
    Lambda_D = 5;
    Lambda_R = (Lambda_D * (l + hub_r))/(totalLength+hub_r);
    Phi = (2/3) * atan(1/Lambda_R);
    %Desired angle of attack for the used profile [degrees]
    A_A = 4;
    %Beta is the blade setting Angle, ie- theta we want to rotate the blade
    Beta = Phi - A_A;
    
    
    theta0=20;    % the starting pitch angle at the beginning in deg
    cl = 0.8 ; % default is 0.8 (80% the length)
    theta=@ (l) -(theta0/(cl*totalLength))*l+theta0;
    
    
    %Begin blade generation
    n=n+1;
    %     create a new filename
    nstr=num2str(n); % changing the n to string
    filename=strcat('xyz',nstr,'.xls');
    
    %     calculations
    %This calculates the chord thickness for the current profile
    thickC = Chord_Thick(l,totalLength,chord_max);
    Test_chord(i) = thickC;
    tempXYZ=thickC*xyz;
    z(n,1)=l;
    
    % this is for aligning the centroids of the airfoils
    Ct=centroid(tempXYZ); % Download it from Mathworks website
    
    %This section of code will keep the leading edge straight
    x_move = abs((chord_max - thickC)/2);
    xShift=[(-1) * x_move 0];
    x_test(i) = xShift(1,1);

    
    if l<cl*totalLength % if the l is less than cut-off length
        twist_test = Beta
        twist_exist = theta(l)
        twistAng=theta(l)*pi/180;
        M=[cos(twistAng) -sin(twistAng); % rotation matrix
            sin(twistAng) cos(twistAng)];
    else
        M=eye(2);
    end
    %         shift center of mass
    tempXYZ=tempXYZ+repmat((C0+xShift-Ct),lxyz,1);
    %         rotate by twistAng
    tempXYZ=M*(tempXYZ');
    %         final value  subrated the shift of the center of mass
    
    %Note, l has been made -l to create the mirror of the blade
    profile(n,:,:)=[tempXYZ' -l*ones(lxyz,1)];
    profileLine(n,:)=[profile(n,1,1) profile(n,1,2) -l];
    tempM(:,:)=profile(n,:,:); % this is after realization that dl
    %     dlmwrite does not work with 3D matrice
    dlmwrite(filename,tempM,'delimiter','\t','newline','pc'); % write a file
    %      import these files into Solidworks as curves
end

%These will plot the chord size and the x shift versus the profile number
%plot(Test_chord)
%pause
%plot(x_test)

dlmwrite('profileLine.xls',profileLine,'delimiter','\t','newline','pc');
end

function [Chord_size] = Chord_Thick(l,totalLength,chord_max)

%Parameters denoting the information for the root of the blade
Root_P = 0.0625; %Percent where the root ends
Root_T_S = 0.50; %Percent of max. chord width at beginning of root
Root_T_E = 0.50; %Percent of max. chord width at end of root

%Parameters denoting the information for the expanding of the blade
Expand_P = 0.1975; %Percent where the expanding ends
Expand_T_E = 1; %Percent of max. chord width at end of root

%Parameters denoting the information for the shrinking of the blade
Shrink_P = 1; %Percent where the shrinking ends
Shrink_T_E = 0.50; %Percent of max. chord width at end of root

x = [0 Root_P Expand_P Shrink_P];
v = [Root_T_S Root_T_E Expand_T_E Shrink_T_E];

xq = (l/totalLength);
Chord_size = chord_max * interp1(x,v,xq,'pchip');

end