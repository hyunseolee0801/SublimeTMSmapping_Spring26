
%% Connect to Scope (This as scope commands)
    % This should work regardless of the computer you are working on
freq=30*10^3;
    scope = scopeconnect; 
     fprintf(scope,':TRIG:PULS:SWE SING');
     T=5*10-6;          
          fprintf(scope,strcat(':CHAN2:DISP OFF')); %removes chan 2 from display
          fprintf(scope,strcat(':CHAN1:OFFS 0')); %0 offset
          fprintf(scope,strcat(':CHAN1:PROB 1'));%amplification factor
          fprintf(scope,strcat(':CHAN1:SCAL?'));%what is the scale
          scal=fscanf(scope)%read scale


          fprintf(scope,strcat(':TIM:SCAL 0.000002'));%set new time scale
          fprintf(scope,':TIM:SCAL?')%write time scale
            tcale = fscanf(scope)

          fprintf(scope,strcat(':TRIG:PULS:SOUR EXT'));%set trigger
          fprintf(scope,strcat(':ACQ:TYPE AVER'));%set acquision type
          fprintf(scope,strcat(':ACQ:AVER 8'));%set number of averages
          fprintf(scope,strcat(':ACQ:AVER?'));%read out number of averages
          Navg = fscanf(scope)
          fprintf(scope,strcat(':ACQ:SAMP?'));%Figure out the sampling rate
          Srate = str2double(fscanf(scope));
          Tsamp=Srate/freq
          
          fprintf(scope,':RUN');
     fprintf(scope,':TRIGger:STATus?');%get the trigger status
        STOP = fscanf(scope);

%% Connect to Arduino
    % You must go back into the code for the probe object and change the
    % port reference for the computer you are working on
    a = probe;

     
   

    
%% Load Map Points
    % Loads the coordinate points for moving the E-field Probe
  % load C:/Users/ljg24/Desktop/postproc/generatecoil/cartcoopt.mat;
   
    %pids=[288 424 54 60 59 39 40 41 1 1000 972 977 976  979 604 603  751  750 740];
    load S1P
    load S2P
    A = S1P;
    B = S2P;
     load  C:\Users\ljg24\Desktop\postproc\generatecoil\cartcoopt.mat cartco order bottom top;
     A = bottom(order);
     B = top(order);
%     
    pids=1:numel(A);

%% Record Data
    % Moves probe (a) to each point designated by S1P and S2P and records
    % the number o f samples d esignated by variable samples.
    samples = 1;
    nsamp=numel(A);
    rawb=zeros([samples       16384        nsamp]);
    for i=1:1
    [raw1, raw2] = proberecord(a, A, B, samples);
    rawb(:,:,:)=reshape(raw1,[samples       16384        nsamp]);
    save datasetflip7C30.mat rawb Tsamp
    end
%% Process Data
    % processes the raw data, does some error checking and filtering and
    % returns the x, y, z components of the E-field at each point
   % [uu, vv, ww] = FieldData(raw1, raw2, 70,0);
 load datasetflip7phi30.mat
 rawa=sum(rawb,1);
 load datasetflip7C30.mat
 rawb=sum(rawb,1);
 pids=1:numel(A);
[dataa,datab]=cleandata(rawa,rawb,Tsamp);
save probeflip7data.mat dataa datab;
%!!!!!switch data values index back from through 10 to through 1000
%%

%% [uu,vv,ww]=cleandatashortcut(dataa,datab);
     load  C:\Users\ljg24\Desktop\postproc\generatecoil\cartcoopt.mat cartco order bottom top;

     A = bottom(order);
     B = top(order);
addpath('C:/Users/ljg24/Desktop/postproc/generatecoil')
load  probeflip7data.mat dataa datab;
pids=1:numel(A);
[cartco,probevec1,probevec2]=arduino2cart(A,B);
datab=datab-mean(datab);
dataa=dataa-mean(datab);
[uu,vv,ww]=cleandataluis(dataa,datab,pids,cartco,probevec1,probevec2);
X = cartco(:,1);
Y = cartco(:,2);
Z = cartco(:,3);


%%

datab=(datab)/max(abs(dataa));
dataa=(dataa)/max(abs(dataa));
f=scatteredInterpolant(X,Y,dataa);
plot(ye,f(ye,xe))
hold on
plot(ye,f(xe,ye))
f2=scatteredInterpolant(X,Y,datab);
plot(ye,f2(xe,ye))
plot(ye,f2(ye,xe))

%%
TR=delaunay(X,Y);
load simres.mat res1 res2;
datab=abs(datab)/max(abs(datab));
dataa=abs(dataa)/max(abs(datab));
res1=abs(res1)/max(abs(res1));
res2=abs(res2)/max(abs(res2));
subplot(3,2,1),
trisurf(TR,X,Y,Z,(dataa(:)),'facecolor','interp','edgealpha',0)
view(2)
subplot(3,2,2),
trisurf(TR,X,Y,Z,(datab(:)),'facecolor','interp','edgealpha',0)
view(2)
subplot(3,2,3),
trisurf(TR,X,Y,Z,(res1(order)),'facecolor','interp','edgealpha',0)
view(2)
subplot(3,2,4),
trisurf(TR,X,Y,Z,(res2(order)),'facecolor','interp','edgealpha',0)
view(2)

subplot(3,2,5),
trisurf(TR,X,Y,Z,(res1(order))-dataa(:),'facecolor','interp','edgealpha',0)
view(2)
subplot(3,2,6),
trisurf(TR,X,Y,Z,(res2(order))-(datab(:)),'facecolor','interp','edgealpha',0)
view(2)
for i=1:6
subplot(3,2,i),
axis equal
colorbar
end
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
[Xc,Yc]=ndgrid(-1:.0005:1,-1:.0005:1);
[C,h]=contour(.07*Xc,.07*Yc,fm(Xc,Yc)/max(mag),[.5 .6 .8 .9]);
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
%light
%lighting gouraud
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

