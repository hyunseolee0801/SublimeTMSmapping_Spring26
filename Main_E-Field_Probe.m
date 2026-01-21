
%% Connect to Scope
    % This should work regardless of the computer you are working on
    scope = scopeconnect; 

%% Connect to Arduino
    % You must go back into the code for the probe object and change the
    % port reference for the computer you are working on
    a = probe;
    
   

    
%% Load Map Points
    % Loads the coordinate points for moving the E-field Probe
    load S1P;
    load S2P;
    A = S1P(1:1000);
    B = S2P(1:1000);
    
%% Record Data
    % Moves probe (a) to each point designated by S1P and S2P and records
    % the number o f samples d esignated by variable samples.
    samples = 20;
    [raw1, raw2] = proberecord(a, A, B, samples);

%% Process Data
    % processes the raw data, does some error checking and filtering and
    % returns the x, y, z components of the E-field at each point
   % [uu, vv, ww] = FieldData(raw1, raw2, 70,0);
    [uu,vv,ww]=postprocraw(raw1, raw2, 70,0);
    %!!!!!switch data values index back from through 10 to through 1000
    
%% Load simulated data

%simvec=xlsread('NHP_283mm_apr13_o2.xlsx',-1);
    
    
%% Rotate

load coordinates

X = x;
Y = y*cos(-pi/2)-z*sin(-pi/2);
Z = y*sin(-pi/2)+z*cos(-pi/2);

max(sqrt(uu.^2+vv.^2+ww.^2))

UU = uu;
VV = vv*cos(-pi/2)-ww*sin(-pi/2);
WW = vv*sin(-pi/2)+ww*cos(-pi/2);

mag = sqrt(uu.^2+vv.^2+ww.^2);
% mag = [0; mag];
% MA = 10*[0:10];
% figure()
% plot(MA,mag)
% hold on
% scatter(MA,mag)
% plot(MA,MA)

% radius = 0.0283;
% X = X.*radius;
% Y = Y.*radius;
% Z = Z.*radius-0.00028952;

% % X1 = -Y;
% % Y1 = X;
% % Z1 = Z;
% % 
% % X2 = Y;
% % Y2 = -X;
% % Z2 = Z;
% % 
% % X3=X2;
% % Y3=Y2;
% % Z3=Z2;

% csvwrite('NHP_283mm_apr13_o1.csv',[X1 Y1 Z1]);
% csvwrite('NHP_283mm_apr13_o2.csv',[X2 Y2 Z2]);
% figure
% scatter3(X,Y,Z,'o');
% xlabel('x');
% ylabel('y');

% 
% uu=simvec(:,1);
% vv=simvec(:,2);
% ww=simvec(:,3);

% % KK = rotmatrix([1,-0.113,0],[2.56,-0.4844,0]/norm([2.56,-0.4844,0])); %remember to change FieldData
% % 
% % X4 = X3*KK(1,1)+Y3*KK(1,2)+Z3*KK(1,3);
% % Y4 = X3*KK(2,1)+Y3*KK(2,2)+Z3*KK(2,3);
% % Z4 = X3*KK(3,1)+Y3*KK(3,2)+Z3*KK(3,3);


%% Plot Data

    % load Michael's Data fron ANSYS
    % current file is LFMS_E-field_ANSYS_7_5_2016
    % Note: the ANSYS data contains x, y, z variables, These are there for
    % reference. Once loaded, you must load the proper coordinates in the
    % proper units from the .mat file below for the visualization to work.
    %load E-field_Coordinates.mat;
    %efieldmap(x, y, z, Ex, Ey, Ez);

    % Make sure you label the graph titles
    % appropriately to not get mixed up between data sets' figures
    
    % Graph the processed e-field components from the probe.
    % Current data file is E-field_LFMS_processed_data_27-Jun-2016
   
    
    %efieldmap(x,y,z,uu,vv,ww,70);  
    efieldmap(X,Y,Z,UU,VV,WW,70);
    % The e-field map function will output the max e-field for the graph. I
    % will usually put a text box on the graph to display that value.
    % Didn't have the time to put that in the code.
    
 %% Other Plots
% i = X.*1000;
% j = Y.*1000;
% k = Z.*1000+37.8;
%  
% Y = y*cos(pi/2)-z*sin(pi/2);
% Z = y*sin(pi/2)+z*cos(pi/2);
% 
% VV = vv*cos(pi/2)-ww*sin(pi/2);
% WW = vv*sin(pi/2)+ww*cos(pi/2);
% diffmap(UU2,VV2,WW2,VX,VY,VZ,i,j,k);
% diffmapavg(UU1,VV1,WW1, UU, VV, WW, X,Y,Z);
absdiffmap(UUsh,VVsh,WWsh, UU, VV, WW, X,Y,Z);
diffmap(UU,VV,WW,zeros(size(UU)),zeros(size(UU)),zeros(size(UU)),X4,Y4,Z4,28.4)
% vecdiff(UU,VV,WW,zeros(size(UU)),zeros(size(UU)),zeros(size(UU)),X4,Y4,Z4,28.4)

