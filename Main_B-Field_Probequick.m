
%% Connect to Scope
    % This should work regardless of the computer you are working on
freq=30*10^3;
    scope = scopeconnect; 
     fprintf(scope,':TRIG:PULS:SWE SING');
     T=5*10-6;          
          fprintf(scope,strcat(':CHAN2:DISP OFF'));
          fprintf(scope,strcat(':CHAN1:OFFS 0'));
          fprintf(scope,strcat(':CHAN1:PROB 1'));
          fprintf(scope,strcat(':CHAN1:SCAL?'));
          scal=fscanf(scope)


          fprintf(scope,':TIM:SCAL?')
     
          fprintf(scope,strcat(':TIM:SCAL 0.000002'));
          fprintf(scope,':TIM:SCAL?')
            tcale = fscanf(scope)

          fprintf(scope,strcat(':TRIG:PULS:SOUR EXT'));
          fprintf(scope,strcat(':ACQ:TYPE AVER'));
          fprintf(scope,strcat(':ACQ:AVER 8'));
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
  % load C:/Users/ljg24/Desktop/postproc/generatecoil/cartcoopt.mat;
   
    %pids=[288 424 54 60 59 39 40 41 1 1000 972 977 976  979 604 603  751  750 740];
    load S1P
    load S2P
    A = S1P;
    B = S2P;
%      load  C:\Users\ljg24\Desktop\postproc\generatecoil\cartcoopt.mat cartco order bottom top;
%      A = bottom(order);
%      B = top(order);
% %     
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
    save Bnormal85.mat rawb Tsamp
    end
%% Process Data
    % processes the raw data, does some error checking and filtering and
    % returns the x, y, z components of the E-field at each point
   % [uu, vv, ww] = FieldData(raw1, raw2, 70,0);
 load Bnormal85.mat
 rawa=sum(rawb,1);
 load Bnormal85.mat
 rawb=sum(rawb,1);
 pids=1:numel(A);
[dataa,datab]=cleandata(rawa,rawb,Tsamp);


save Bnormal85res.mat dataa A B;
%%
load Bnormal85res.mat dataa A B;
addpath('C:/Users/ljg24/Desktop/postproc/generatecoil')
[cartco,probevec1,probevec2]=arduino2cart(A,B);

uu=dataa(:).*cartco(:,1);
vv=dataa(:).*cartco(:,1);
ww=dataa(:).*cartco(:,1);

X = cartco(:,1);
Y = cartco(:,2);
Z = cartco(:,3);

f=scatteredInterpolant(X,Y,dataa);
TR=delaunay(X,Y);
subplot(1,2,1),
trisurf(TR,X,Y,Z,(dataa(:)),'facecolor','interp','edgealpha',0)
axis equal
axis off
%%light
%lighting
[XX,YY]=ndgrid(min(X(:)):.01:max(X(:)),min(Y(:)):.01:max(Y(:)));
view(2)
subplot(1,2,2),contour(XX,YY,abs(f(XX,YY))+abs(f(XX,YY)))
%%

trisurf(TR,0.07*X,0.07*Y,Z,(dataa)/max(dataa),'facecolor','interp','edgealpha',0)
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