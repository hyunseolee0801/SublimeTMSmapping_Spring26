

%% Connect to Scope
    % This should work regardless of the computer you are working on
freq=30*10^3;
    scope = scopeconnect; 
     fprintf(scope,':TRIG:PULS:SWE SING');
     T=5*10-6;          
          fprintf(scope,strcat(':CHAN2:OFFS 0'));
          fprintf(scope,strcat(':CHAN2:PROB 1'));
          fprintf(scope,strcat(':CHAN2:SCAL?'));
          fprintf(scope,strcat(':CHAN1:OFFS 0'));
          fprintf(scope,strcat(':CHAN1:PROB 1'));
          fprintf(scope,strcat(':CHAN1:SCAL?'));
          scal=fscanf(scope)


          fprintf(scope,':TIM:SCAL?')
     
          fprintf(scope,strcat(':TIM:SCAL 0.000005'));
          fprintf(scope,':TIM:SCAL?')
            tcale = fscanf(scope)

          fprintf(scope,strcat(':TRIG:PULS:SOUR EXT'));
          fprintf(scope,strcat(':ACQ:TYPE AVER'));
          fprintf(scope,strcat(':ACQ:AVER 16'));
          fprintf(scope,strcat(':ACQ:AVER?'));
          Navg = fscanf(scope)
          fprintf(scope,strcat(':ACQ:SAMP?'));
          Srate = str2double(fscanf(scope));
          Tsamp=Srate/freq
          
          fprintf(scope,':RUN');
     fprintf(scope,':TRIGger:STATus?');
        STOP = fscanf(scope);

%% Connect to Arduino
    % You must go back into the code for the probe object and change the
    % port reference for the computer you are working on
    a = probe;

     
   

    
%% Load Map Points
    % Loads the coordinate points for moving the E-field Probe
    load S1P;
    load S2P;
    %pids=[288 424 54 60 59 39 40 41 1 1000 972 977 976  979 604 603  751  750 740];
    A = S1P(1:1000);
    B = S2P(1:1000);
%     load  C:\Users\ljg24\Desktop\postproc\generatecoil\cartcoopt.mat cartco order bottom top;
%     A = bottom(order);
%     B = top(order);
%     
    pids=1:numel(A);

%% Record Data
    % Moves probe (a) to each point designated by S1P and S2P and records
    % the number o f samples d esignated by variable samples.
    samples = 1;
    nsamp=numel(A);
    rawb=zeros([samples       16384/2        nsamp]);
    for i=1:1
    [raw1, raw2] = proberecord(a, A, B, samples);
    rawa(:,:,:)=reshape(raw2,[samples       16384/2        nsamp]);
    rawb(:,:,:)=reshape(raw1,[samples       16384/2        nsamp]);

    save data7.mat rawa rawb Tsamp
    end
%% Process Data
    % processes the raw data, does some error checking and filtering and
    % returns the x, y, z components of the E-field at each point
   % [uu, vv, ww] = FieldData(raw1, raw2, 70,0);
   load data7.mat
%  load data6C90corr.mat rawb
%  rawa=sum(rawb,1);
%  load data6phi90corr.mat rawb
%  rawb=sum(rawb,1);
 pids=1:1000;
[dataa,datab]=cleandata(rawa,rawb,Tsamp);
save probe7data.mat dataa datab;
%!!!!!switch data values index back from through 10 to tlhshrough 1000
%%

%% [uu,vv,ww]=cleandatashortcut(dataa,datab);
load  probe7data.mat dataa datab;
pids=1:1000;
[cartco,probevec1,probevec2]=arduino2cart(A,B);
[uu,vv,ww]=cleandataluis(dataa,datab,pids,cartco,probevec1,probevec2);
X = cartco(:,1);
Y = cartco(:,2);
Z = cartco(:,3);

TR=delaunay(X,Y);

subplot(2,2,1),
trisurf(TR,X,Y,Z,abs(dataa(:))/max(dataa),'facecolor','interp','edgealpha',0)
view(2)
subplot(2,2,2),
trisurf(TR,X,Y,Z,abs(datab(:))/max(dataa),'facecolor','interp','edgealpha',0)
view(2)
%%
axis equal

max(sqrt(uu.^2+vv.^2+ww.^2))

 UU = uu;
 VV = vv;
 WW = ww;

mag = sqrt(uu.^2+vv.^2+ww.^2);
figure(1)
subplot(2,2,1),
trisurf(TR,X,Y,Z,abs(UU),'facecolor','interp','edgealpha',0)
axis equal
view(2)
subplot(2,2,2),
trisurf(TR,X,Y,Z,abs(VV),'facecolor','interp','edgealpha',0)
axis equal
view(2)
subplot(2,2,3),
trisurf(TR,X,Y,Z,abs(WW),'facecolor','interp','edgealpha',0)
axis equal
view(2)
subplot(2,2,4),
trisurf(TR,X,Y,Z,abs(mag(:)),'facecolor','interp','edgealpha',0)
xlabel('x')
ylabel('y')
axis equal
view(2)
figure(2)
subplot(2,2,1),
quiver3(X,Y,Z,uu,zeros(size(VV)),zeros(size(VV)));
view(2)
subplot(2,2,2),
quiver3(X,Y,Z,zeros(size(VV)),vv,zeros(size(VV)));
view(2)
subplot(2,2,3),
quiver3(X,Y,Z,zeros(size(VV)),zeros(size(VV)),ww);
view(2)
subplot(2,2,4),
quiver3(X,Y,Z,UU,VV,WW);
view(2)


figure(3)
fm=scatteredInterpolant(X,Y,mag);
[Xc,Yc]=ndgrid(-1:.001:1,-1:.001:1);
[C,h]=contour(.07*Xc,.07*Yc,fm(Xc,Yc)/max(mag),[.4 .5 .6 .7 .8 .9 .95]);
hold on
scatter(0,0)
clabel(C,h)
%%


trisurf(TR,0.07*X,0.07*Y,Z,(mag/max(mag)),'facecolor','interp','edgealpha',0)
colormap(parula)
axis equal
view([0 0 1])
axis off
ax = gca;
mymap = colormap(ax);
newmap=zeros([128 3]);
newmap(:,1)=interp1(sqrt(1:64)',mymap(:,1),1:7/127:8);
newmap(:,2)=interp1(sqrt(1:64)',mymap(:,2),1:7/127:8);
newmap(:,3)=interp1(sqrt(1:64)',mymap(:,3),1:7/127:8);

colormap(newmap);
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

    
    %    efieldmap(x_alt,y_alt,z_alt,UU_STIM_A,VV_STIM_A,WW_STIM_A,70);
%mag=sqrt(UU_STIM_A.^2,VV_STIM_A.^2,WW_STIM_A.^2);
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

