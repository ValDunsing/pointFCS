clear all
close all
% Script for analysis of 1channel pointFCS data

% global movingaveragewindowsize blocksize spatialfilter depletioncorrection...
%     intensityfilter intensitythreshold intensityfilterwindow

% loadandpool=0; % 0 load new raw file, 1 load previously computed ACFs
% calculateACFonly=0; % Calculate ACF only and return afterwards
 
%!!!======================================%!!!
mac=0; %set 1 when working on mac, 0 for windows computer!
%!!!======================================%!!!

ACFselection=1; % 1 divide timetrace in ACFsegmentnumb segments and select or discard segment ACFs
ACFsegmentnumb=5; % for SMOG, VSVG etc.
%ACFsegmentnumb=9; %for Fog 40k
%ACFsegmentnumb=14; %for Fog 60k
depletioncorrection=0;
savebackgroundaverage=0; % 1: calculate average background correction and save to file
backgroundcorrection=1; % 1: load background average and correct by subtraction
cut=1; % set to 0 when calculating bg, removes the first segment of the time trace (e.g. the time trace is divided in ACFsegmentnumb+1 segmets and the first removed)
twocomponents=1; % referes to two diffusing components
binningwindow_GUI=1000;
%binningwindow=10000;
binningwindow=500;
triplett=1;
flickering=0; % refers to an additional photophysics term
timetotal=90; % Total Measurement time in seconds
S=6; %put to 1000 for receptor diffusion
G0=1;
momentsonly=0;
scrsz=   get(0,'ScreenSize');

pixeltime=1.23e-06;
%pixeltime=7.56e-06;
binningtime=pixeltime;
% intensityselection=1;
% global timedelay
%  global goon1 mintimedelay maxtimedelay
% if ACFselection==1
%     global curveincl correlationcurves corfit ItraceII goon
% end
% if intensityselection==1
%     global segmentsincl Ifull timelinebinned numberofsegments segmentlength correctionfitseg ff
% end



% Import: exported .mat file from function overnight

% % % % [nameout, pathout]=uigetfile('*.mat', 'Select mat file to open');
% % % % fprintf(1,'\nImporting out file...');
% % % % Imp=load([pathout nameout], '-ascii');
%Imp=load([pathout nameout])

%timedelay=Imp.DTime;
%timedelay=Imp(:,2);

% Selection of counts to include based on lifetime histogram
% goon1=0;
% backgroundcorrection
% 
% while goon1==0
%     pause(5)
% end

% Read out of raw data

%rohdaten2=Imp.TrueTimens;
% % % % rohdaten2=Imp(:,1);

%rohdaten2(Imp.DTime<mintimedelay | Imp.DTime>maxtimedelay)=NaN;
%rohdaten2(Imp(:,2)<mintimedelay | Imp(:,2)>maxtimedelay)=NaN;
%rohdaten2(isnan(rohdaten2))=[];
%rohdaten=(rohdaten2-min(rohdaten2));%Zeit f???ngt immer von Null an**
%rohdaten(rohdaten==0)=[];
% % % % rohdaten=rohdaten2;

% Binning interval of raw data
% % % % binningtime=5e-6; %in seconds
% % % % pixeltime=binningtime;

% Construction of time series:
% first column: Time, second column: Counts
%Itrace=zeros(1000*uint32(((max(rohdaten)/1e9/binningtime)+1)/1000),2);

%Itrace=zeros((max(rohdaten)/1e9/binningtime)+1,2);
% % % % Itrace=zeros(uint32(timetotal/binningtime),2);
% % % % maxindex=max(rohdaten)/1e9;
% % % % %Itrace(:,1)=0:binningtime:maxindex;
% % % % Itrace(:,1)=0:binningtime:timetotal-binningtime;
% % % % for i=1:size(rohdaten,1)
% % % % 
% % % % 
% % % % Itrace(ceil(rohdaten(i)/1e9/binningtime),2)=Itrace(ceil(rohdaten(i)/1e9/binningtime),2)+1;
% % % % end
%return
path= uigetdir; 
path2=uigetdir('D:\Shared\Valentin\LSM Zeiss Data');
files=dir([path '/*.tif']);

if backgroundcorrection==1
    [backgroundavgfile,filepath_bg]=uigetfile;
end


for ibatch=1:size(files,1)
    if ACFselection==1
    global curveincl correlationcurves corfit ItraceII goon
    end
    
    [namedata,remain]=strtok(files(ibatch).name,'.');
    inputfilename=files(ibatch).name;
%[inputfilename, path]=uigetfile('D:\Shared\Valentin\*.tif');
% tic


%Loading of inputimage (.tif)
inputimage=double(imread([path '/' inputfilename]));
rows=size(inputimage,1);
coloumns=size(inputimage,2);


%Creating vector from inputimage that shows intensity
vectorofintensity=zeros(1,rows*coloumns);
for i=1:rows
    vectorofintensity((i-1)*coloumns+1:i*coloumns)=inputimage(i,:);
end



%Creating vector for the time
vectoroftime=zeros(1,(rows*coloumns));
for i=1:(rows*coloumns);
    vectoroftime(1,i)=i*(1.53/(10^6));
end




% % % % fluorescenceseries=Itrace(:,2);
fluorescenceseries=vectorofintensity';
if cut==1
    index_cut=length(fluorescenceseries)/(ACFsegmentnumb+1)+1;
    fluorescenceseries=fluorescenceseries(index_cut:end);
end

if savebackgroundaverage==1
if ibatch==1
    backgroundintensities=zeros(length(fluorescenceseries),size(files,1));
end
backgroundintensities(:,ibatch)=fluorescenceseries;
end

if backgroundcorrection==1
    backgroundavgcorrection=load([filepath_bg '\' backgroundavgfile]);
    bg_correction=zeros(size(fluorescenceseries));
    if cut==0
        bg_correction=backgroundavgcorrection.correctionfit_bg;
    else
        bg_correction=backgroundavgcorrection.correctionfit_bg_cut;
    end
    fluorescenceseries=fluorescenceseries-bg_correction;
    if mean(fluorescenceseries)<mean(bg_correction)
        sprintf('Background correction causes negative average value for file:')
        inputfilename
    end
        
        
    %if mean(fluorescenceseries(length(fluorescenceseries)/2+1:end))<0 || mean(fluorescenceseries(length(fluorescenceseries)/2+1:end))>mean(fluorescenceseries(1:length(fluorescenceseries)/2))
    %    fluorescenceseries=fluorescenceseries+bg_correction;
    %    sprintf('Background correction causes negative values and is not applied for file:')
    %    inputfilename
    %end
        
        
end

timeline=1:1:length(fluorescenceseries);
timeline=timeline'*pixeltime;
clear Itrace

bleachingfraction=1-nanmean(fluorescenceseries(end-1000+1:end))/nanmean(fluorescenceseries(1:1000));

% Figure: Binned fluorescence time series and Brightness
lastbinlength=mod(length(fluorescenceseries),binningwindow)+binningwindow;
fluorescenceseriesbinned=mean(reshape(fluorescenceseries(1:end-lastbinlength),binningwindow,length(fluorescenceseries(1:end-lastbinlength))/binningwindow),1);
fluorescenceseriesbinned=[fluorescenceseriesbinned mean(fluorescenceseries(end-lastbinlength+1:end))];

timeline_binned=1:1:length(fluorescenceseriesbinned);
timeline_binned=(timeline_binned-1)*binningwindow*pixeltime;
%multiexpfit=@(a,b,c,ta,tb,tc,x)a*exp(ta*x)+b*exp(tb*x)+c*exp(tc*x);
f=fit(timeline_binned',fluorescenceseriesbinned','exp2');
%f=fit(timeline_binned',fluorescenceseriesbinned',multiexpfit,'StartPoint',[0.8 -0.5 0.1 -0.1 0.1 -0.05]);
correctionfit=f.a.*exp(f.b*timeline)+f.c.*exp(f.d*timeline);
%correctionfit=f.a.*exp(f.ta*timeline)+f.b.*exp(f.tb*timeline)+f.c.*exp(f.tc*timeline);

correctionfitbinned=mean(reshape(correctionfit(1:end-lastbinlength),binningwindow,length(correctionfit(1:end-lastbinlength))/binningwindow),1);
correctionfitbinned=[correctionfitbinned mean(correctionfit(end-lastbinlength+1:end))];
fluorescencebinnedresiduals=fluorescenceseriesbinned-correctionfitbinned;
    
%     linebrightnessbinned=var(reshape(linefluorescenceseries(1:end-lastbinlength),binningwindow,length(linefluorescenceseries(1:end-lastbinlength))/binningwindow),1);
%     linebrightnessbinned=[linebrightnessbinned var(linefluorescenceseries(end-lastbinlength+1:end))];
%     meanbrightness=mean(linebrightnessbinned);
%     meanbrightnessfit=meanbrightness*ones(1,length(linebrightnessbinned));
%     linebrightnessresiduals=linebrightnessbinned-meanbrightness;
hI=figure('OuterPosition',[scrsz(1) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Binned Membrane Fluorescence time series');
plot(timeline_binned,fluorescenceseriesbinned)
hold on
plot(f)
xlabel('time in s')
ylabel('fluorescence Ch1')

fluorescenceseries_uncorr=fluorescenceseries;

if depletioncorrection==1
    %linefluorescenceseries_corrected=depletioncorrectionfunc(timeline,linefluorescenceseries,fsin);
    fluorescenceseries_corrected=depletioncorrectionfunc(timeline,fluorescenceseries,f); %Funktion umschreiben, sodass sie auf Fit zur?ckgreift!!
    lastbinlength=mod(length(fluorescenceseries),binningwindow)+binningwindow;
    fluorescenceseries_corrected_binned=mean(reshape(fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    fluorescenceseries_corrected_binned=[fluorescenceseries_corrected_binned mean(fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('OuterPosition',[scrsz(1) scrsz(2) 2*scrsz(3)/3 scrsz(4)/2],'Name','Corrected Membrane Fluorescence time series')
    plot(timeline_binned,fluorescenceseries_corrected_binned);
    xlabel('time in s')
    ylabel('fluorescence Ch1')
    
    fluorescenceseries=fluorescenceseries_corrected;
end

% Calculate quantiles to check for negative intensity offset
n=length(fluorescenceseries)/binningwindow;
Q=quantile(fluorescenceseriesbinned,n);
min(Q)
hQ=figure('OuterPosition',[scrsz(1) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Intensity Quantiles');
plot([1:1:n],Q)
xlabel('Quantile number')
ylabel('I') 
saveas(hQ,[path2 '\' inputfilename(1:end-4) '_I_Quantiles.fig'])
 
% if Q(1)<0
%     %fluorescenceseries=fluorescenceseries-Q(1);
%     %fluorescenceseries=fluorescenceseries-max(Q(Q<0));
%     fluorescenceseries=fluorescenceseries;
%     sprintf('Background correction causes negative values:')
%     inputfilename
%     %fluorescenceseries_uncorr=fluorescenceseries-Q(1);
% end
    
save([path2 '\' inputfilename(1:end-4) '_fluorescenceseries_bleachcorr' '.mat'], 'fluorescenceseries');
save([path2 '\' inputfilename(1:end-4) '_fluorescenceseries' '.mat'], 'fluorescenceseries_uncorr');




fluorescenceseries_moments=fluorescenceseries(1/10*length(fluorescenceseries)+1:end);
%fluorescenceseries_moments=fluorescenceseries;
brightness=nanvar(fluorescenceseries_moments)/mean(fluorescenceseries_moments)-1;
brightnessS=nanvar(fluorescenceseries_moments)/mean(fluorescenceseries_moments)-1.06;
centralmoments=zeros(1,6);
centralmoments(1)=mean(fluorescenceseries_moments);
for i=1:5
    centralmoments(i+1)=cemo(fluorescenceseries_moments,i);
end

if momentsonly==1
    nameout=inputfilename;
    fid3=fopen([path2 '\' nameout(1:end-4) '_final_Sfactor.txt'],'a'); % adjust path if necessary!
    fid4=fopen([path2 '\' nameout(1:end-4) '_final_centralmoments.txt'],'a'); % adjust path if necessary!
    
    outputSfactor=[0 0];
    outputcentralmoments=centralmoments;
         
    fprintf(fid3,'%e\n',outputSfactor');
    fprintf(fid4,'%e\n',outputcentralmoments');
    
    saveas(hI,[path2 '\' nameout(1:end-4) ' I_t.fig'])
else

fprintf('Computing ACF...\n');

[tcorr,fcorr,sigmas]=autocorrpFCSmultipletau(fluorescenceseries);
tcorr=tcorr*pixeltime;

scrsz=   get(0,'ScreenSize');
figure('OuterPosition',[2*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Autocorrelation function')
    % semilogx(tcorr,fcorr,'b');
    % xlabel('time')
    % ylabel('autocorrelation')
    % hold on
semilogx(tcorr,fcorr,'r.')
figure('OuterPosition',[2*scrsz(3)/3 scrsz(2) scrsz(3)/3 scrsz(4)/2],'Name','Weights')
weights=fcorr./sigmas;
semilogx(tcorr,weights,'b.')
    
    
% % Fit of correlation function with diffusion model
% fprintf('Fitting...\n');
% lb=pixeltime;
% %ub=1;
% ub=max(tcorr);
% tfit=tcorr;
% %S=5;
% %tfit=lb:linetime:ub;
% %tfit=tcorr2(1:40);
% 
% % Weights
% weights=abs(weights);
% %weights=abs(fcorr2./sigmas);
% % weights=zeros(size(tfit));
% % for i=1:length(tfit)
% %     weights(i)=1/i;
% % end    
% tcorrfit=tcorr(1:length(tfit));
% fcorrfit=fcorr(1:length(tfit));
% weightsfit=weights(1:length(tfit));
% 
% if triplett==1
%     lsfitfunc=@(x,t)((1-x(4)+x(4).*exp(-t./x(5)))./(1-x(4))).*(1./x(1)).*((1+t./x(2)).^-1).*(1+t./(x(2).*x(3).^2)).^-0.5;
%     x0=[500,0.001, S, 0.2, 3*10^-6];
%     fixed=[false false true false false];
%     [N,taud,Sfit,Tfraction,Ttau,CI,fitcurve,residuals] = autocorrfit3Ddifftrip(tcorrfit,fcorrfit,lb,ub,weightsfit,lsfitfunc,x0,fixed);
% else
%     lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
%     x0=[50,0.001 S];
%     fixed=[false false true];
%     [N,taud,Sfit,CI,fitcurve,residuals] = autocorrfit3Ddiff(tcorrfit,fcorrfit,lb,ub,weightsfit,lsfitfunc,x0,fixed);
% end
% 
% % Figure: Plot of data and fit
% %N=10;
% %taud=0.002;
% figure('OuterPosition',[scrsz(3) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Fit curve')
% subplot(2,1,1),semilogx(tcorrfit,fitcurve,'-r')
% hold on
% subplot(2,1,1),semilogx(tcorrfit,fcorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')
% subplot(2,1,2),semilogx(tcorrfit,residuals,'bx')
% xlabel('time')
% ylabel('residuals')
% hold on
% subplot(2,1,2),semilogx(tcorrfit,zeros(size(tcorrfit)),'-k');

% Again: Binning of corrected intensity trace for GUI
lastbinlength=mod(length(fluorescenceseries),binningwindow_GUI)+binningwindow_GUI;
fluorescenceseriesbinned_GUI=mean(reshape(fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
fluorescenceseriesbinned_GUI=[fluorescenceseriesbinned_GUI mean(fluorescenceseries(end-lastbinlength+1:end))];

pixeltime=binningtime;
    %segmentlength=1000*uint32(length(fluorescenceseries)/(ACFsegmentnumb*1000));
    segmentlength=(length(fluorescenceseries)/(ACFsegmentnumb));
    if G0==1
        G0s=zeros(1,ACFsegmentnumb);
    end
    %Isegmentlength_GUI=uint32(segmentlength/binningwindow_GUI);
     Isegmentlength_GUI=segmentlength/binningwindow_GUI;
    Itrace=zeros(segmentlength,ACFsegmentnumb);
    Itrace_GUI=zeros(Isegmentlength_GUI,ACFsegmentnumb);
    %correlationcurves=zeros(1,ACFsegmentnumb);
    %sigmascurves=correlationcurves;
    %corfit=zeros(size(correlationcurves,1),size(correlationcurves,2))
    %x0=[N,taud,Sfit];
    for i=1:size(Itrace,2)
        Itrace(:,i)=fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace_GUI(:,i)=fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        if G0==1
            [tcorriG0,fcorriG0,sigmasiG0]=autocorrpFCSmultipletauG0(Itrace(:,i));
            tcorri=tcorriG0(2:end);
            fcorri=fcorriG0(2:end);
            sigmasi=sigmasiG0(2:end);
            G0s(i)=fcorriG0(1);
        else
            [tcorri,fcorri,sigmasi]=autocorrpFCSmultipletau(Itrace(:,i));
        end
        tcorri=tcorri*pixeltime;
        correlationcurves(:,i)=fcorri';
        sigmascurves(:,i)=sigmasi';
        lbi=min(tcorri);
        ubi=max(tcorri);
%         [Ni,taudi,Sfiti,CIi,fitcurvei,residualsi] = autocorrfit3Ddiff(tcorri*linetime,fcorri,lbi,ubi,1./sigmasi,lsfitfunc,x0);
%         x0=[Ni,taudi,Sfiti];
%         corfit(:,i)=fitcurvei';
        
    end
    
    fprintf(' Done Computing CFS (segments)!\n');
    figure('OuterPosition',[scrsz(1) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],'Name','ACFs segments')
    for i=1:size(correlationcurves,2)
        fcorri=correlationcurves(:,i)';
        semilogx(tcorri,fcorri,'.b')
%         hold on
%         semilogx(tcorri,fitcurvei,'-r') 
        pause
        hold on
    end
        
    
    corfit=zeros(length(tcorri),ACFsegmentnumb);
    
%     % Binning for better display
%     Binfactor=100;
%     ItraceII=zeros(size(Itrace,1)/Binfactor,size(Itrace,2));
%     for i=1:size(Itrace,2)
%         ItraceII(:,i)=nanmean(reshape(Itrace(:,i),Binfactor,length(Itrace(:,1))/Binfactor),1);
%     end


%    ItraceII=Itrace;
    ItraceII=Itrace_GUI;
%    segmenttime=linetime:linetime:linetime*segmentlength;
    segmenttime_GUI=pixeltime:pixeltime*binningwindow_GUI:pixeltime*segmentlength;
%    ItraceII=[segmenttime' ItraceII];
    ItraceII=[segmenttime_GUI' ItraceII];
%     segmenttime_binned=linetime:Binfactor*linetime:linetime*segmentlength;
%     ItraceII=[segmenttime_binned' ItraceII];
    correlationcurves=[tcorri' correlationcurves];
    corfit=[tcorri' corfit];
    curveincl=ones(size(corfit,2)-1,1);

    % GUI for selection of segment ACFs
    Corrselection
    goon=0;
    while goon==0
        pause(5)
    end
    fprintf('Done!')
    fluorescenceseries2=fluorescenceseries;
    for i=1:length(curveincl)
        if curveincl(i)==0
            fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
        end
    end
    
    % Calculation of final ACF
    % ------> whole time trace
    if ACFselection==0
    [tcorrfinal,fcorrfinal,sigmasfinal]=autocorrpFCSmultipletau(fluorescenceseries);
    tcorrfinal=tcorrfinal*pixeltime;
    weightsfinal=abs(fcorrfinal./sigmasfinal);
    end
    
    % ------> average of segments
    fcorrsegfinal=zeros(length(tcorri),1);
    sigmassegfinal=fcorrsegfinal;
    weightseg=fcorrsegfinal;
    if G0==1
        G0final=0;
    end
    for i=1:length(curveincl);
        if curveincl(i)==1
             fcorrsegfinal=fcorrsegfinal+correlationcurves(:,1+i);
            %fcorrsegfinal=fcorrsegfinal+correlationcurves(:,i+1)./sigmascurves(:,i);
             sigmassegfinal=sigmassegfinal+sigmascurves(:,i);
             if G0==1
                 G0final=G0final+G0s(i);
             end
            %sigmassegfinal=sigmassegfinal+sigmascurves(:,i);
            %weightsseg=weightseg+1./sigmascurves(:,i);
        end
    end
    fcorrsegfinal=fcorrsegfinal/sum(curveincl);
    %fcorrsegfinal=fcorrsegfinal./weightsseg;
    sigmassegfinal=sigmassegfinal/sum(curveincl);
    %fcorrsegfinal=nansum(correlationcurves(:,2:end)./sigmascurves,2);
    %fcorrsegfinal=fcorrsegfinal./nansum(sigmascurves,2);
    %sigmassegfinal=nanmean(sigmascurves,2);
    weightssegfinal=abs(fcorrsegfinal./sigmassegfinal);
    tcorrsegfinal=tcorri;
    
    if G0==1
        G0final=G0final/sum(curveincl);
    end
    
    
    % Fit of final correlation function with diffusion model
    fprintf('Final Fitting...\n');
    lbseg=pixeltime;
    %ub=1;
    ubseg=max(tcorrsegfinal);
    
    %tfit=lb:linetime:ub;
    %tfit=tcorr2(1:40);

    % Weights
    % weights=zeros(size(tfit));
    % for i=1:length(tfit)
    %     weights(i)=1/i;
    % end    
    
    if ACFselection==0
        ub=max(tcorrfinal);
        lb=pixeltime;
        tfit=tcorrfinal;
        tcorrfit=tcorrfinal(1:length(tfit));
        fcorrfit=fcorrfinal(1:length(tfit));
        weightsfit=weightsfinal(1:length(tfit));
        weightsfit(weightsfit==0)=10^-4;
    end
    
    fcorrsegfinal=fcorrsegfinal(5:end);
    weightssegfinal=weightssegfinal(5:end);
    tcorrsegfinal=tcorrsegfinal(5:end);
    
   
    if triplett==1
        if twocomponents==0
            if flickering==0
                alpha=1;
                lsfitfunc=@(x,t)((1-x(4)+x(4).*exp(-t./x(5)))./(1-x(4))).*(1./x(1)).*((1+(t./x(2)).^x(6)).^-1).*(1+(t./(x(2))).^x(6).*x(3).^-2).^-0.5;
                fitsuccess=0;
                failed=0;
                while fitsuccess<1
                if failed==0
                    N0=1/(mean(fcorrsegfinal(1:3)));
                    %x0=[N0,10^-2, S, 0.1, 20*10^-6, alpha];
                    %lbfit=[0, 10^-5, S, 0, 1*10^-6, alpha];
                    %lbfit=[0, 10^-5, S, 0, 0.1*10^-6, 0.5*alpha];
                    %ubfit=[100, 60*10^-6, S, 0.9, 10*10^-6, alpha];
                    %ubfit=[5000, 1, S, 0.9, 50*10^-6, 1.25*alpha];

                    %x0=[N0,10^-3, S, 0.1, 50*10^-6, alpha];
                    x0=[N0,1*10^-4, S, 0.1, 20*10^-6, alpha];
                    %lbfit=[0, 10^-5, S, 0, 10*10^-6, 0.5*alpha];
                    lbfit=[0, 5*10^-5, S, 0, 0.1*10^-6, alpha];
                    ubfit=[5000, 100, S, 0.9, 100*10^-6, alpha];
                    %ubfit=[5000, 100, S, 0.9, 200*10^-6, alpha];
                end
                fitshow=figure('Name','Find initial conditions...');
                fitfunctry=lsfitfunc(x0,tcorrsegfinal);
                semilogx(tcorrsegfinal,fitfunctry,'--g')
                hold on
                semilogx(tcorrsegfinal,fcorrsegfinal,'.g')
                fixed=[false false false false false false];

                if ACFselection==0
                [Nfinal,taudfinal,Sfitfinal,Tfractionfinal,Ttaufinal,CIfinal,fitcurvefinal,residualsfinal] = autocorrfit3Ddifftrip(tcorrfit,fcorrfit,lb,ub,weightsfit,lsfitfunc,x0,fixed);
                else
                    try    
                        [Nsegfinal,taudsegfinal,Sfitsegfinal,Tfractionsegfinal,Ttausegfinal,alphasegfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit3Ddifftripbounds(tcorrsegfinal,fcorrsegfinal',lbfit,ubfit,weightssegfinal',lsfitfunc,x0,fixed);
                        fitsuccess=1;
                        sprintf('Fitting successful')
                    catch
                    failed=1;
                    sprintf('Adjust initial parameters!')
                    if fitsuccess==0
                        x01=input('Initial parameters [N1,taud1,S,T1,tauT1,alpha]?');
                    end
                    close(fitshow)
                    end
                end
                end
            else
                alpha=1;
                lsfitfuncflick=@(x,t)((1-x(4)+x(4).*exp(-t./x(5)))./(1-x(4))).*(1./x(1)).*((1+(t./x(2)).^x(6)).^-1).*(1+(t./(x(2))).^x(6).*x(3).^-2).^-0.5;
              
            end
        else
            alpha=1;
            ffree=0.2;
            
            
%             tcorrsegfinal=tcorrsegfinal(5:end);
%             fcorrsegfinal=fcorrsegfinal(5:end);
%             weightssegfinal=weightssegfinal(5:end);
            
            lsfitfunctwocomp=@(x,t)((1-x(6)+x(6).*exp(-t./x(7)))./(1-x(6))).*(1./x(1)).*(x(2).*((1+(t./x(3)).^x(8)).^-1).*(1+(t./(x(3))).^x(8).*x(5).^-2).^-0.5+(1-x(2)).*((1+(t./x(4)).^x(9)).^-1).*(1+(t./(x(4))).^x(9).*x(5).^-2).^-0.5);
            lsfitfunctwocomp2D=@(x,t)((1-x(6)+x(6).*exp(-t./x(7)))./(1-x(6))).*(1./x(1)).*(x(2).*((1+(t./x(3)).^x(8)).^-1).*(1+(t./(x(3))).^x(8).*x(5).^-2).^-0.5+(1-x(2)).*((1+(t./x(4)).^x(9)).^-1));
            N0=1/(mean(fcorrsegfinal(1:3)));
            x0=[N0,ffree,3*10^-4,5*10^-3, S, 0.1, 30*10^-6, alpha,alpha];
            
            %x0=[N0,ffree,3*10^-4,1*10^-3, S, 0.1, 30*10^-6, alpha,alpha];   % for receptor analysis
      
            lbfit=[0,0,5*10^-5,10^-3, S, 0, 1*10^-7, alpha,alpha];
    
            ubfit=[5000,1,5*10^-3,10^-1, S, 0.9, 100*10^-6, alpha, alpha];
            
            fixed=[false false false false false false false false false];
            
% %             x0=[N0,ffree,85*10^-6,10^-3, S, 0.1, 5*10^-6, alpha, alpha];
% %             
% %             lbfit=[0,0,85*10^-6,4*10^-4, S, 0, 1*10^-7, alpha, alpha];
% %     
% %             ubfit=[5000,1,85*10^-6,10^-1, S, 0.9, 10*10^-6, alpha, alpha];
% %             
% %             fixed=[false false false false false false false false false];
            %%
            % Two components diffusing in 3D:
            %[Nsegfinal,ffreefit,taud1segfinal,taud2segfinal,Sfitsegfinal,Tfractionsegfinal,Ttausegfinal,alpha1segfinal,alpha2segfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit3Dtwocompdifftripbounds(tcorrsegfinal,fcorrsegfinal',lbfit,ubfit,weightssegfinal',lsfitfunctwocomp,x0,fixed);
            %%
            % Two components where one is 3D diffusion, the other 2D
            % diffusion on apical membrane
            [Nsegfinal,ffreefit,taud1segfinal,taud2segfinal,Sfitsegfinal,Tfractionsegfinal,Ttausegfinal,alpha1segfinal,alpha2segfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit2Dtwocompdifftripbounds(tcorrsegfinal,fcorrsegfinal',lbfit,ubfit,weightssegfinal',lsfitfunctwocomp2D,x0,fixed);
        end
    else
    lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    x0=[5,0.01 S];
    fixed=[false false true];
    % ------> Fit whole time trace
    if ACFselection==0
    [Nfinal,taudfinal,Sfitfinal,CIfinal,fitcurvefinal,residualsfinal] = autocorrfit3Ddiff(tcorrfit,fcorrfit,lb,ub,weightsfit,lsfitfunc,x0,fixed);
    else
    % ------> Fit of segment average
    [Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit3Ddiff(tcorrsegfinal,fcorrsegfinal',lbseg,ubseg,weightssegfinal',lsfitfunc,x0,fixed);
    end
    end
    
    if twocomponents==0
    usseg=0.5*(CIsegfinal(1:3,2)-CIsegfinal(1:3,1));
    
    if ACFselection==0
        us=0.5*(CIfinal(1:3,2)-CIfinal(1:3,1));

    % Figure: Plot of data and fit
    %N=10;
    %taud=0.002;
    %-----> whole time trace
    h=figure('OuterPosition',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2],'Name','Final Fit curve');
    subplot(2,1,1),semilogx(tcorrfit,fitcurvefinal,'-r')
    hold on
    subplot(2,1,1),semilogx(tcorrfit,fcorrfit,'bx')
    xlabel('time')
    ylabel('autocorrelation')
    textN=['N = %.1f ' char(177) ' %.1f'];
    strtextN=sprintf(textN,Nfinal,us(1));
    texttaud=[' = %.3f ' char(177) ' %.3f'];
    strtexttaud=['\tau' sprintf(texttaud,taudfinal,us(2)) ' s'] ;
%     textw0=['w0 = %.2f ' char(177) ' %.2f'];
%     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
    dim = [0.55 0.6 0.3 0.3];
%     str = {strtextC,strtextD,strtextw0};
    str = {strtextN,strtexttaud};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    subplot(2,1,2),semilogx(tcorrfit,residualsfinal,'bx')
    xlabel('time')
    ylabel('residuals')
    hold on
    subplot(2,1,2),semilogx(tcorrfit,zeros(size(tcorrfit)),'-k');
    end
    
    %-----> segment average
    if ACFselection==1
    hh=figure('OuterPosition',[scrsz(1) scrsz(4)/4 2*scrsz(3)/3 3*scrsz(4)/4],'Name','Final (segment averaged) Fit curve');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    set(groot,'defaultLineLineWidth',1.5)
    subplot('Position',positionvector1),semilogx(tcorrsegfinal,fitcurvesegfinal,'-r')
    hold on
    subplot('Position',positionvector1),semilogx(tcorrsegfinal,fcorrsegfinal,'bx')
    set(gca,'FontSize',18,'FontWeight','bold','LineWidth',1.0,'XGrid','on','YGrid','on');
    xlabel('time','FontSize',24,'FontWeight','bold')
    ylabel('Autocorrelation','FontSize',24,'FontWeight','bold')
    
    subplot('Position',positionvector2),semilogx(tcorrsegfinal,residualssegfinal,'bx')
    set(gca,'FontSize',18,'FontWeight','bold','LineWidth',1.0,'XGrid','on','YGrid','on');
    xlabel('time','FontSize',24,'FontWeight','bold')
    ylabel('residuals','FontSize',24,'FontWeight','bold')
    hold on
    subplot('Position',positionvector2),semilogx(tcorrsegfinal,zeros(size(tcorrsegfinal)),'-k');
    textNseg=['N = %.1f ' char(177) ' %.1f'];
    strtextNseg=sprintf(textNseg,Nsegfinal,usseg(1));
    texttaudseg=[' = %.2f ' char(177) ' %.2f'];
    strtexttaudseg=['\tau' sprintf(texttaudseg,taudsegfinal*10^3,usseg(2)*10^3) ' ms'] ;
    textSseg=['S = %.1f '];
    strtextSseg=sprintf(textSseg,Sfitsegfinal);
    textalpha1seg=['alpha = %.2f '];
    strtextalpha1seg=sprintf(textalpha1seg,alphasegfinal);
%     textw0=['w0 = %.2f ' char(177) ' %.2f'];
%     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
    dim = [0.55 0.55 0.3 0.3];
%     str = {strtextC,strtextD,strtextw0};
    strseg = {strtextNseg,strtexttaudseg,strtextSseg,strtextalpha1seg};
    t=annotation('textbox',dim,'String',strseg,'FitBoxToText','on');
    set(t,'FontSize',28)
    end
    else
        usseg=0.5*(CIsegfinal(1:5,2)-CIsegfinal(1:5,1));
        hh=figure('OuterPosition',[scrsz(1) scrsz(4)/4 2*scrsz(3)/3 3*scrsz(4)/4],'Name','Final (segment averaged) Fit curve two components');
        positionvector1=[0.1 0.35 0.8 0.55];
        positionvector2=[0.1 0.1 0.8 0.15];
        set(groot,'defaultLineLineWidth',1.5)
        subplot('Position',positionvector1),semilogx(tcorrsegfinal,fitcurvesegfinal,'-r')
        hold on
        subplot('Position',positionvector1),semilogx(tcorrsegfinal,fcorrsegfinal,'bx')
        set(gca,'FontSize',18,'FontWeight','bold','LineWidth',1.0,'XGrid','on','YGrid','on');
        xlabel('time','FontSize',24,'FontWeight','bold')
        ylabel('Autocorrelation','FontSize',24,'FontWeight','bold')

        subplot('Position',positionvector2),semilogx(tcorrsegfinal,residualssegfinal,'bx')
        set(gca,'FontSize',18,'FontWeight','bold','LineWidth',1.0,'XGrid','on','YGrid','on');
        xlabel('time','FontSize',24,'FontWeight','bold')
        ylabel('residuals','FontSize',24,'FontWeight','bold')
        hold on
        subplot('Position',positionvector2),semilogx(tcorrsegfinal,zeros(size(tcorrsegfinal)),'-k');
        textNseg=['N = %.1f ' char(177) ' %.1f'];
        strtextNseg=sprintf(textNseg,Nsegfinal,usseg(1));
        textffreeseg=[' = %.2f ' char(177) ' %.2f'];
        strtextffreeseg=['f free' sprintf(textffreeseg,ffreefit,usseg(2))] ;
        
        texttaud1seg=[' = %.2f ' char(177) ' %.2f'];
        strtexttaud1seg=['\tau_1' sprintf(texttaud1seg,taud1segfinal*10^3,usseg(3)*10^3) ' ms'] ;
        texttaud2seg=[' = %.2f ' char(177) ' %.2f'];
        strtexttaud2seg=['\tau_2' sprintf(texttaud2seg,taud2segfinal*10^3,usseg(4)*10^3) ' ms'] ;
        textSseg=['S = %.1f '];
        strtextSseg=sprintf(textSseg,Sfitsegfinal);
        textalpha1seg=['alpha_1 = %.2f '];
        strtextalpha1seg=sprintf(textalpha1seg,alpha1segfinal);
        textalpha2seg=['alpha_2 = %.2f '];
        strtextalpha2seg=sprintf(textalpha2seg,alpha2segfinal);
    %     textw0=['w0 = %.2f ' char(177) ' %.2f'];
    %     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
        dim = [0.55 0.55 0.3 0.3];
    %     str = {strtextC,strtextD,strtextw0};
        strseg = {strtextNseg,strtextffreeseg,strtexttaud1seg,strtexttaud2seg,strtextSseg,strtextalpha1seg,strtextalpha2seg};
        t=annotation('textbox',dim,'String',strseg,'FitBoxToText','on');
        set(t,'FontSize',28)
    end
    
    
    % Save results
    %path2=uigetdir('D:\Shared\Valentin\LSM Zeiss Data');
    nameout=inputfilename;
    if mac==1
        fid1=fopen([path2 '/' nameout(1:end-4) '_final_ACF.txt'],'a'); % adjust path if necessary!
        fid2=fopen([path2 '/' nameout(1:end-4) '_final_fitparameters.txt'],'a'); % adjust path if necessary!
    
        saveas(hh,[path2 '/' nameout(1:end-4) ' ACF.fig'])
    else
    fid1=fopen([path2 '\' nameout(1:end-4) '_final_ACF.txt'],'a'); % adjust path if necessary!
    fid2=fopen([path2 '\' nameout(1:end-4) '_final_fitparameters.txt'],'a'); % adjust path if necessary!
    fid3=fopen([path2 '\' nameout(1:end-4) '_final_Sfactor.txt'],'a'); % adjust path if necessary!
    fid4=fopen([path2 '\' nameout(1:end-4) '_final_centralmoments.txt'],'a'); % adjust path if necessary!
    fid5=fopen([path2 '\' nameout(1:end-4) '_final_curveincl.txt'],'a'); % adjust path if necessary!
    fid6=fopen([path2 '\' nameout(1:end-4) '_bg_offset.txt'],'a'); % adjust path if necessary!

    saveas(hh,[path2 '\' nameout(1:end-4) ' ACF.fig'])
    end
    
    %I0=nanmean(fluorescenceseries);
    meanintensity_corrected=nanmean(fluorescenceseries);
    
    Sfactorfit=(G0final-1/Nsegfinal)*meanintensity_corrected;
    Sfactoramp=(G0final-fcorrsegfinal(1))*meanintensity_corrected;
    
    output=zeros(length(tcorrsegfinal),3);
    output(:,1)=tcorrsegfinal;
    output(:,2)=fcorrsegfinal;
    output(:,3)=weightssegfinal;
    if twocomponents==0
    if triplett==1
        outputparameters=zeros(9,3);
        outputparameters(6,1)=alphasegfinal;
        outputparameters(7,1)=Tfractionsegfinal;
        outputparameters(8,1)=Ttausegfinal;
        outputparameters(9,1)=fluorescenceseriesbinned(1);
    else
        outputparameters=zeros(5,3);
    end
    outputparameters(1:5,1)=[Nsegfinal;taudsegfinal;Sfitsegfinal;bleachingfraction;meanintensity_corrected];
    outputparameters(1:3,2:end)=CIsegfinal(1:3,:);
  
         
    
    else
        outputparameters=zeros(12,3);
        outputparameters(6,1)=alpha1segfinal;
        outputparameters(7,1)=alpha2segfinal;
        outputparameters(8,1)=Tfractionsegfinal;
        outputparameters(9,1)=Ttausegfinal;
        outputparameters(10,1)=bleachingfraction;
        outputparameters(11,1)=meanintensity_corrected;
        outputparameters(12,1)=fluorescenceseriesbinned(1);
        outputparameters(1:5,1)=[Nsegfinal;ffreefit;taud1segfinal;taud2segfinal;Sfitsegfinal];
        outputparameters(1:5,2:end)=CIsegfinal(1:5,:);
    end
    outputSfactor=[Sfactoramp Sfactorfit];
    outputcentralmoments=centralmoments;
    fprintf(fid1,'%e\t %e\t %e\n',output');
    fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
    fprintf(fid3,'%e\n',outputSfactor');
    fprintf(fid4,'%e\n',outputcentralmoments');
    fprintf(fid5,'%e\n',curveincl');
    fprintf(fid6,'%e\n',Q(1));
end
close all
fclose all;
clearvars -global curveincl correlationcurves corfit ItraceII goon
clearvars sigmascurves
end

if savebackgroundaverage==1
    binningwindow_bg=100;
    backgroundintensity_avg=mean(backgroundintensities,2);
    lastbinlength=mod(length(backgroundintensity_avg),binningwindow_bg)+binningwindow_bg;
    backgroundintensity_avgbinned=mean(reshape(backgroundintensity_avg(1:end-lastbinlength),binningwindow_bg,length(backgroundintensity_avg(1:end-lastbinlength))/binningwindow_bg),1);
    backgroundintensity_avgbinned=[backgroundintensity_avgbinned mean(backgroundintensity_avg(end-lastbinlength+1:end))];

    timeline_binned_bg=1:1:length(backgroundintensity_avgbinned);
    timeline_binned_bg=(timeline_binned_bg-1)*binningwindow_bg*pixeltime;

    figure('Name','Background Average')
    plot(timeline_binned_bg,backgroundintensity_avgbinned)
    hold on
    fbg=fit(timeline_binned_bg',backgroundintensity_avgbinned','exp2');
    plot(fbg)
    correctionfit_bg=fbg.a.*exp(fbg.b*timeline)+fbg.c.*exp(fbg.d*timeline);
    save([path2 '\background'],'correctionfit_bg');
end
