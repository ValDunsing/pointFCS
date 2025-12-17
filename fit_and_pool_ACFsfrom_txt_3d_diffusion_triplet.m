% Fit programm for .txt ACF files exported from ZEISS ZEN FCS window

S=6;
triplett=1;
Zen=0;
scrsz=   get(0,'ScreenSize');

path= uigetdir; 
files=dir([path '/*.txt']);

[namedata,remain]=strtok(files(1).name,'.');
namefile=files(1).name;
ACFdata1=load([path '/' namefile]);
    %tcorrsegfinal=ACFdata(:,1);
    %fcorrsegfinal=ACFdata(:,2)-1;
tcorrsegfinal1=ACFdata1(1:end,1);
if Zen==1
    fcorrsegfinal1=ACFdata1(1:end,2)-1;
    weightssegfinal1=(1:1:length(tcorrsegfinal1))/length(tcorrsegfinal1)';
else
    fcorrsegfinal1=ACFdata1(1:end,2);
    weightssegfinal1=ACFdata1(1:end,3);
end

ACFs=zeros(size(ACFdata1,1),size(files,1));
weights=zeros(size(ACFdata1,1),size(files,1));
ACFs(:,1)=fcorrsegfinal1;
weights(:,1)=weightssegfinal1;

%fitparameterlist=zeros(1,6);
hh=figure('OuterPosition',[scrsz(1) scrsz(4)/4 2*scrsz(3)/3 3*scrsz(4)/4],'Name','Final Fit curve');
for i=2:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    ACFdata=load([path '/' namefile]);
    %tcorrsegfinal=ACFdata(:,1);
    %fcorrsegfinal=ACFdata(:,2)-1;
    tcorrsegfinal=ACFdata(1:end,1);
    if Zen==1
        ACFs(:,i)=ACFdata(1:end,2)-1;
        weights(:,i)=(1:1:length(tcorrsegfinal))/length(tcorrsegfinal)';
    else
        ACFs(:,i)=ACFdata(1:end,2);
        weights(:,i)=ACFdata(1:end,3);
    end
    %weightssegfinal=ones(size(tcorrsegfinal));
    
end


fcorrsegfinal=mean(ACFs,2);
weightssegfinal=mean(weights,2);
%weightssegfinal=ones(length(tcorrsegfinal),1);
%weightssegfinal=1./(1:1:length(tcorrsegfinal))';

% fcorrsegfinal=fcorrsegfinal(2:end);
% weightssegfinal=weightssegfinal(2:end);
% tcorrsegfinal=tcorrsegfinal(2:end);

fcorrsegfinal=fcorrsegfinal(5:end);
weightssegfinal=weightssegfinal(5:end);
tcorrsegfinal=tcorrsegfinal(5:end);

if triplett==1
    alpha=1;
    lsfitfunc=@(x,t)((1-x(4)+x(4).*exp(-t./x(5)))./(1-x(4))).*(1./x(1)).*((1+(t./x(2)).^x(6)).^-1).*(1+(t./(x(2))).^x(6).*x(3).^-2).^-0.5;
    fitsuccess=0;
    failed=0;
    while fitsuccess<1
    if failed==0
        N0=1/(mean(fcorrsegfinal(1:3)));
        %x0=[N0, 10^-4, S, 0.2, 8.41*10^-6, alpha];
        %lbfit=[0, 20^-5, S, 0, 8.41*10^-6, alpha];
        %ubfit=[5000, 0.01, S, 0.9, 8.41*10^-6, alpha];
        
        x0=[N0,3*10^-4, S, 0.2,50*10^-6, alpha];
        lbfit=[0,5*10^-5, S, 0, 0.1*10^-6, alpha];
        ubfit=[5000, 100, S, 0.9, 100*10^-6, alpha];
        %ubfit=[5000, 100, S, 0.9, 6*10^-6, alpha];
    end
%             fitshow=figure('Name','Find initial conditions...');
    fitfunctry=lsfitfunc(x0,tcorrsegfinal);
%             semilogx(tcorrsegfinal,fitfunctry,'--g')
%             hold on
%             semilogx(tcorrsegfinal,fcorrsegfinal,'.g')
    fixed=[false false false false false false];
    try    
        [Nsegfinal,taudsegfinal,Sfitsegfinal,Tfractionsegfinal,Ttausegfinal,alphasegfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit3Ddifftripbounds(tcorrsegfinal',fcorrsegfinal',lbfit,ubfit,weightssegfinal',lsfitfunc,x0,fixed);
        %[Nsegfinal,taudsegfinal,Sfitsegfinal,Tfractionsegfinal,Ttausegfinal,alphasegfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit3Ddifftrip(tcorrsegfinal',fcorrsegfinal',lbfit,ubfit,weightssegfinal',lsfitfunc,x0,fixed);

        fitsuccess=1
        sprintf('Fitting successful')
    catch
    failed=1;
    sprintf('Adjust initial parameters!')
    if fitsuccess==0
        x01=input('Initial parameters [N1,taud1,S,T1,tauT1,alpha]?');
    end
%             close(fitshow)
    end
    end
fitparameterlist=[Nsegfinal taudsegfinal Sfitsegfinal Tfractionsegfinal Ttausegfinal alphasegfinal];
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
    fitparameterlis=[Nsegfinal taudsegfinal Sfitsegfinal 0 0 0];
end

usseg=0.5*(CIsegfinal(1:3,2)-CIsegfinal(1:3,1));

%-----> segment average

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
textalphaseg=['alpha = %.2f '];
strtextalphaseg=sprintf(textalphaseg,alphasegfinal);
%     textw0=['w0 = %.2f ' char(177) ' %.2f'];
%     strtextw0=[sprintf(textw0,w0final,us(3)) ' \mu m'];
dim = [0.55 0.55 0.3 0.3];
%     str = {strtextC,strtextD,strtextw0};
strseg = {strtextNseg,strtexttaudseg,strtextSseg,strtextalphaseg};
t=annotation('textbox',dim,'String',strseg,'FitBoxToText','on');
set(t,'FontSize',28)


 output=zeros(length(tcorrsegfinal),3);
 output(:,1)=tcorrsegfinal;
 output(:,2)=fcorrsegfinal;
 output(:,3)=weightssegfinal;


path2=uigetdir;
filename=sprintf('fog_lp0p25_2.txt'); %Does that also work for windows, or do we need \ there?
fid1=fopen([path2 '\' filename],'a'); % adjust path if necessary!
fid2=fopen([path2 '\' filename(1:end-4) '_avg_ACF.txt'],'a'); % adjust path if necessary!

saveas(hh,[path2 '\' filename '_ACF.fig'])


fprintf(fid1,'sample\t N\t tau\t S\t T\t tauT\t alpha\n');
    %fprintf(fid1,files(i).name(1:end-24));
fprintf(fid1,files(1).name);
fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\n',fitparameterlist');
fprintf(fid2,'%e\t %e\t %e\n',output');
fclose all;