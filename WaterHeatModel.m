% Matlab/Octave code for temperature modelling video. 
%
% Samuli Siltanen September 2019

%% Preliminary computations and definitions

% Parameters for controlling plot appearance
msize = 20;
lwidth = 1.5;
predline = 1.5;
thinline = 1;
fsize = 16;
tempfsize = 24;
GrayLine = [.6 .6 .6];
LightGray = [.85 .85 .85];

% Here are the data points. Time in minutes, temperature in C.
timevec_orig = [1,2,3,4,5,6,7,8,9,10,11,12,13,...
    14,15,16,17,18,19,20,21,22,...
    23,24,25,26,27,28,29,30,31,32,33,34,...
    35,36,37,38,39,40,41,42,43,44,...
    45,46,47,48,49,50,51,52,53,54,55,...
    56,57,58,59,60,61,62,63,64,65,...
    66,67,68,69,70,71,72,73,74,75,...
    76,77,78,79,80,81,82,83,84,85,...
    86,87,88,89,90,91,92,93,94,95,...
    96,97,98,99,100,101,102,103,104,105,...
    106,107,108,109,110,111,112,113,114,115,116,...
    117,118,119,120,121,122,123,124,125,...
    126,127,128,129,130,131,132,133,134,135,136,...
    137,138,139,140,141,142,143,144,145,146];
tempvec0_orig = [88,84,81,78,75,73,70,68,66,64,62,61,59,...
    58,57,55,54,53,52,51,50,49,...
    48,48,47,46,45,45,44,44,43,42,42,41,...
    41,41,40,40,39,39,38,38,38,37,...
    37,37,36,36,36,35,35,35,35,34,34,...
    34,34,33,33,33,33,33,32,32,32,...
    32,32,31,31,31,31,31,31,31,30,...
    30,30,30,30,30,30,30,30,29,29,...
    29,29,29,29,29,29,29,29,29,28,...
    28,28,28,28,28,28,28,28,28,28,...
    28,28,28,28,27,27,27,27,27,27,27,...
    27,27,27,27,27,27,27,27,27,...
    27,27,27,27,27,26,26,26,26,26,26,...
    26,26,26,26,26,26,26,26,26,26];

% Let's pick out a representative part of the measurements.
timevec = timevec_orig(5:10:round(.55*end));
tempvec0 = tempvec0_orig(5:10:round(.55*end));

% Record an early time-temperature pair
earlymeastime = timevec_orig(15);
earlymeastemp = tempvec0_orig(15);

% Record prediction time and temperature
latemeastime = timevec_orig(90);
latemeastemp = tempvec0_orig(90);

% Background temperature (best guess; the room had different temperatures
% away from the boiling kettle...)
roomtemp = 26;

% Ideal time coordinate vector, extended clearly outside the measurement
% time vector. That way we can study predictions of the model for earlier
% and later times than those measured.
tMAX = max(timevec);
idtimevec = linspace(0,2*tMAX,1000);
idtimevecneg = linspace(-tMAX,2*tMAX,1000);


%% Build an ideal model by fitting a polynomial model to log-temperatures

% Subtract room temperature
tempvec = tempvec0-roomtemp;

% Linear fit for exponential heat behaviour
logtempvec = log(tempvec);
p = polyfit(timevec,logtempvec,1);

% Measured values
figure(1)
clf
plot(timevec,tempvec0,'r.','markersize',msize)
hold on
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
set(gca,'xtick',0:10:80 ,'fontsize',fsize)
ylim([20 100])
xlim([0 85])
pbaspect([2 1 1])
box off

% Plot log of measured points
figure(2)
clf
plot(timevec,logtempvec,'r.','markersize',msize)
set(gca,'xtick',0:10:80 ,'fontsize',fsize)
xlim([0 85])
pbaspect([2 1 1])
box off

% Plot linear fit to log-points
figure(3)
clf
plot(timevec,polyval(p,timevec),'k','linewidth',lwidth)
hold on
plot(timevec,logtempvec,'r.','markersize',msize)
set(gca,'xtick',0:10:80 ,'fontsize',fsize)
xlim([0 85])
pbaspect([2 1 1])
box off

% Illustrate the exponential model fit
figure(4)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
%plot(timevec,roomtemp+exp(polyval(p,timevec)),'k','linewidth',lwidth)
plot(timevec,roomtemp+exp(p(2))*exp(p(1)*timevec),'k','linewidth',lwidth)
plot(timevec,tempvec0,'r.','markersize',msize)
set(gca,'xtick',0:10:80 ,'fontsize',fsize)
ylim([20 100])
xlim([0 85])
pbaspect([2 1 1])
box off




%% Study sensitivity of forward prediction.

% Choose time instant to which the prediction goes
predtime = latemeastime;
indpred = min(find(abs(idtimevec-predtime)==min(abs(idtimevec-predtime))));

% Given the time-temperature pair at some early measurement time,
% we determine the exponential function passing through that point
% and decaying to room temperature with the same rate than the above model.
timeA = earlymeastime;
tempA = earlymeastemp;
logtempA = log(tempA-roomtemp);
pA = p;
pA(2) = logtempA-pA(1)*timeA; % Find constant in the linear model for logtemp
indA = min(find(abs(idtimevec-timeA)==min(abs(idtimevec-timeA))));

timeB = earlymeastime-2;
tempB = earlymeastemp-2;
logtempB = log(tempB-roomtemp);
pB = p;
pB(2) = logtempB-pB(1)*timeB; % Find constant in the linear model for logtemp
indB = min(find(abs(idtimevec-timeB)==min(abs(idtimevec-timeB))));

timeC = earlymeastime+2;
tempC = earlymeastemp+2;
logtempC = log(tempC-roomtemp);
pC = p;
pC(2) = logtempC-pC(1)*timeC; % Find constant in the linear model for logtemp
indC = min(find(abs(idtimevec-timeC)==min(abs(idtimevec-timeC))));

% Evaluate predictions at last measured time
predictionA = roomtemp+exp(pA(2)+pA(1)*predtime);
expfunB = roomtemp+exp(pB(2)+pB(1)*idtimevec);
predictionB = roomtemp+exp(pB(2)+pB(1)*predtime);
expfunC = roomtemp+exp(pC(2)+pC(1)*idtimevec);
predictionC = roomtemp+exp(pC(2)+pC(1)*predtime);
truepred = tempvec0(end);


figure(5)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
% Plot exponential fit
plot(idtimevec(indA:indpred),roomtemp+exp(pA(2)+pA(1)*idtimevec(indA:indpred)),'k','linewidth',predline)
plot(timeA,tempA,'r.','markersize',msize)
plot(predtime,predictionA,'k.','markersize',msize)
text(predtime-7,predictionA+10,[num2str(predictionA,'%0.1f')],'fontsize',tempfsize)
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off


figure(6)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
% Plot exponential fit
pp = plot(idtimevec(indA:indpred),roomtemp+exp(pA(2)+pA(1)*idtimevec(indA:indpred)),'k','linewidth',predline);
set(pp,'color',GrayLine)
plot(idtimevec(indB:indpred),expfunB(indB:indpred),'k','linewidth',predline)
plot(timeA,tempA,'r.','markersize',msize)
plot(timeB,tempB,'b.','markersize',msize)
pp = plot(predtime,predictionA,'k.','markersize',msize);
set(pp,'color',GrayLine)
plot(predtime,predictionB,'k.','markersize',msize)
text(predtime-7,predictionA+10,[num2str(predictionB,'%0.1f')],'fontsize',tempfsize)
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off


figure(7)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(indA:indpred),roomtemp+exp(pA(2)+pA(1)*idtimevec(indA:indpred)),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevec(indB:indpred),expfunB(indB:indpred),'k','linewidth',predline);
set(pp,'color',GrayLine)
plot(idtimevec(indC:indpred),expfunC(indC:indpred),'k','linewidth',predline)
pp = plot(timeA,tempA,'r.','markersize',msize);
pp = plot(timeB,tempB,'b.','markersize',msize);
set(pp,'color',GrayLine)
plot(timeC,tempC,'b.','markersize',msize)
pp = plot(predtime,predictionA,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(predtime,predictionB,'k.','markersize',msize);
set(pp,'color',GrayLine)
plot(predtime,predictionC,'k.','markersize',msize)
text(predtime-7,predictionA+10,[num2str(predictionC,'%0.1f')],'fontsize',tempfsize)
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off


figure(8)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(indA:indpred),roomtemp+exp(pA(2)+pA(1)*idtimevec(indA:indpred)),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevec(indB:indpred),expfunB(indB:indpred),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevec(indC:indpred),expfunC(indC:indpred),'k','linewidth',predline)
set(pp,'color',GrayLine)
pp = plot(timeA,tempA,'r.','markersize',msize);
pp = plot(timeB,tempB,'b.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(timeC,tempC,'b.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(predtime,predictionA,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(predtime,predictionB,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(predtime,predictionC,'k.','markersize',msize);
set(pp,'color',GrayLine)
plot(latemeastime,latemeastemp,'r.','markersize',msize)
tt = text(predtime-7,predictionA+10,[num2str(latemeastemp,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',[1 0 0])
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off





%% Study sensitivity of backward prediction over a long time interval.

% Given the time-temperature pair at the late measurement time,
% we determine the exponential function passing through that point
% and decaying to room temperature with the same rate than the above model.
timeG = latemeastime;
tempG = latemeastemp;
logtempG = log(tempG-roomtemp);
pG = p;
pG(2) = logtempG-pG(1)*timeG; % Find constant in the linear model for logtemp
indG = min(find(abs(idtimevec-timeG)==min(abs(idtimevec-timeG))));

timeH = latemeastime-1;
tempH = latemeastemp-1;
logtempH = log(tempH-roomtemp);
pH = p;
pH(2) = logtempH-pH(1)*timeH; % Find constant in the linear model for logtemp
indH = min(find(abs(idtimevec-timeH)==min(abs(idtimevec-timeH))));

timeI = latemeastime+1;
tempI = latemeastemp+1;
logtempI = log(tempI-roomtemp);
pI = p;
pI(2) = logtempI-pI(1)*timeI; % Find constant in the linear model for logtemp
indI = min(find(abs(idtimevec-timeI)==min(abs(idtimevec-timeI))));

% Evaluate predictions at time zero
predictionG = roomtemp+exp(pG(2));
expfunG = roomtemp+exp(pG(2)+pG(1)*idtimevec);
predictionH = roomtemp+exp(pH(2));
expfunH = roomtemp+exp(pH(2)+pH(1)*idtimevec);
predictionI = roomtemp+exp(pI(2));
expfunI = roomtemp+exp(pI(2)+pI(1)*idtimevec);


figure(19)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(timeG,tempG,'r.','markersize',msize);
tt = text(timeG,tempG+8,[num2str(tempG,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',[1 0 0])
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off


figure(20)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(1:indG),roomtemp+exp(pG(2)+pG(1)*idtimevec(1:indG)),'k','linewidth',predline);
tt = text(2,predictionG+3,[num2str(predictionG,'%0.0f')],'fontsize',tempfsize);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(0,predictionG,'k.','markersize',msize);
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off



figure(21)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(1:indG),roomtemp+exp(pG(2)+pG(1)*idtimevec(1:indG)),'k','linewidth',predline);
set(pp,'color',GrayLine)
tt = text(2,predictionG+3,[num2str(predictionG,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
pp = plot(idtimevec(1:indH),expfunH(1:indH),'k','linewidth',predline);
tt = text(2,predictionH+3,[num2str(predictionH,'%0.0f')],'fontsize',tempfsize);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeH,tempH,'b.','markersize',msize);
pp = plot(0,predictionG,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(0,predictionH,'k.','markersize',msize);
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off



figure(22)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(1:indG),roomtemp+exp(pG(2)+pG(1)*idtimevec(1:indG)),'k','linewidth',predline);
set(pp,'color',GrayLine)
tt = text(2,predictionG+3,[num2str(predictionG,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
pp = plot(idtimevec(1:indH),expfunH(1:indH),'k','linewidth',predline);
set(pp,'color',GrayLine)
tt = text(2,predictionH+3,[num2str(predictionH,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
pp = plot(idtimevec(1:indI),expfunI(1:indI),'k','linewidth',predline)
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeH,tempH,'b.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(timeI,tempI,'b.','markersize',msize);
pp = plot(0,predictionG,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(0,predictionH,'k.','markersize',msize);
set(pp,'color',GrayLine)
tt = text(20,80,[num2str(predictionI,'%0.0f')],'fontsize',tempfsize);
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off





figure(23)
clf
plot(idtimevec,roomtemp*ones(size(idtimevec)),'b','linewidth',thinline)
hold on
pp = plot(idtimevec(1:indG),roomtemp+exp(pG(2)+pG(1)*idtimevec(1:indG)),'k','linewidth',predline);
set(pp,'color',GrayLine)
tt = text(2,predictionG+3,[num2str(predictionG,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
pp = plot(idtimevec(1:indH),expfunH(1:indH),'k','linewidth',predline);
set(pp,'color',GrayLine)
tt = text(2,predictionH+3,[num2str(predictionH,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
pp = plot(idtimevec(1:indI),expfunI(1:indI),'k','linewidth',predline)
set(pp,'color',GrayLine)
pp = plot(timeH,tempH,'b.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(timeI,tempI,'b.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(0,predictionG,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(0,predictionH,'k.','markersize',msize);
set(pp,'color',GrayLine)
% pp = plot(predtime,predictionI,'k.','markersize',msize);
% set(pp,'color',GrayLine)
tt = text(20,80,[num2str(predictionI,'%0.0f')],'fontsize',tempfsize);
set(tt,'color',GrayLine)
plot(0,100,'r.','markersize',msize)
set(gca,'xtick',0:10:130,'fontsize',fsize)
ylim([20 100])
xlim([0 120])
pbaspect([2 1 1])
box off






%% Evaluate estimates for the boiling time counted backward from 100 minutes

expfunG = roomtemp+exp(pG(2)+pG(1)*idtimevecneg);
tmpG = abs(roomtemp+exp(pG(2)+pG(1)*idtimevecneg)-100);
indG100 = find(tmpG == min(tmpG));
timeG100 = idtimevecneg(indG100);

expfunH = roomtemp+exp(pH(2)+pH(1)*idtimevecneg);
tmpH = abs(roomtemp+exp(pH(2)+pH(1)*idtimevecneg)-100);
indH100 = find(tmpH == min(tmpH));
timeH100 = idtimevecneg(indH100);

expfunI = roomtemp+exp(pI(2)+pI(1)*idtimevecneg);
tmpI = abs(roomtemp+exp(pI(2)+pI(1)*idtimevecneg)-100);
indI100 = find(tmpI == min(tmpI));
timeI100 = idtimevecneg(indI100);

indGneg = min(find(abs(idtimevecneg-timeG)==min(abs(idtimevecneg-timeG))));
indHneg = min(find(abs(idtimevecneg-timeH)==min(abs(idtimevecneg-timeH))));
indIneg = min(find(abs(idtimevecneg-timeI)==min(abs(idtimevecneg-timeI))));



figure(31)
clf
plot(idtimevecneg,roomtemp*ones(size(idtimevecneg)),'b','linewidth',thinline)
hold on
pp = plot(idtimevecneg(1:indGneg),roomtemp+exp(pG(2)+pG(1)*idtimevecneg(1:indGneg)),'k','linewidth',predline);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeG100,100,'k.','markersize',msize);
pp = plot([timeG100,timeG100],[0,100],'k--','linewidth',thinline);
plot(0,100,'r.','markersize',msize)
pp = plot([0,0],[0,100],'r--','linewidth',thinline);
set(gca,'xtick',-30:10:130,'fontsize',fsize)
ylim([20 100])
xlim([-30 110])
pbaspect([2 1 1])
box off



figure(32)
clf
plot(idtimevecneg,roomtemp*ones(size(idtimevecneg)),'b','linewidth',thinline)
hold on
pp = plot(idtimevecneg(1:indGneg),roomtemp+exp(pG(2)+pG(1)*idtimevecneg(1:indGneg)),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevecneg(1:indHneg),expfunH(1:indHneg),'k','linewidth',predline);
pp = plot(timeH,tempH,'b.','markersize',msize);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeG100,100,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot([timeG100,timeG100],[0,100],'k--','linewidth',thinline);
set(pp,'color',GrayLine)
pp = plot([timeH100,timeH100],[0,100],'k--','linewidth',thinline);
pp = plot(timeH100,100,'k.','markersize',msize);
plot(0,100,'r.','markersize',msize)
pp = plot([0,0],[0,100],'r--','linewidth',thinline);
set(gca,'xtick',-30:10:130,'fontsize',fsize)
ylim([20 100])
xlim([-30 110])
pbaspect([2 1 1])
box off



figure(33)
clf
plot(idtimevecneg,roomtemp*ones(size(idtimevecneg)),'b','linewidth',thinline)
hold on
pp = plot(idtimevecneg(1:indGneg),roomtemp+exp(pG(2)+pG(1)*idtimevecneg(1:indGneg)),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevecneg(1:indHneg),expfunH(1:indHneg),'k','linewidth',predline);
set(pp,'color',GrayLine)
pp = plot(idtimevecneg(1:indIneg),expfunI(1:indIneg),'k','linewidth',predline);
pp = plot(timeH,tempH,'b.','markersize',msize);
pp = plot(timeI,tempI,'b.','markersize',msize);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeG100,100,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot([timeG100,timeG100],[0,100],'k--','linewidth',thinline);
set(pp,'color',GrayLine)
pp = plot([timeH100,timeH100],[0,100],'k--','linewidth',thinline);
set(pp,'color',GrayLine)
pp = plot([timeI100,timeI100],[0,100],'k--','linewidth',thinline);
pp = plot(timeH100,100,'k.','markersize',msize);
set(pp,'color',GrayLine)
pp = plot(timeI100,100,'k.','markersize',msize);
plot(0,100,'r.','markersize',msize)
pp = plot([0,0],[0,100],'r--','linewidth',thinline);
set(gca,'xtick',-30:10:130,'fontsize',fsize)
ylim([20 100])
xlim([-30 110])
pbaspect([2 1 1])
box off



figure(34)
clf
pp = patch([timeH100,timeH100,timeI100,timeI100],[0,100,100,0],LightGray);
set(pp,'EdgeColor',LightGray)
hold on
plot(idtimevecneg,roomtemp*ones(size(idtimevecneg)),'b','linewidth',thinline)
pp = plot(idtimevecneg(1:indHneg),expfunH(1:indHneg),'k','linewidth',predline);
pp = plot(idtimevecneg(1:indIneg),expfunI(1:indIneg),'k','linewidth',predline);
pp = plot(timeH,tempH,'b.','markersize',msize);
pp = plot(timeI,tempI,'b.','markersize',msize);
pp = plot(timeG,tempG,'r.','markersize',msize);
pp = plot(timeH100,100,'k.','markersize',msize);
pp = plot(timeI100,100,'k.','markersize',msize);
plot(0,100,'r.','markersize',msize)
pp = plot([0,0],[0,100],'r--','linewidth',thinline);
set(gca,'xtick',-30:10:130,'fontsize',fsize)
ylim([20 100])
xlim([-30 110])
pbaspect([2 1 1])
box off


%% Display results
disp(' ')
disp(['Temperature at future time ',num2str(latemeastime),' is between ',num2str(predictionB), ' and ',num2str(predictionC)])
disp(['Long-time backward prediction: Temperature at time zero is between ',num2str(predictionH), ' and ',num2str(predictionI)])
disp(['Long-time backward prediction: Boiling time was between ',num2str(round(timeH100)),' and ',num2str(round(timeI100))])



