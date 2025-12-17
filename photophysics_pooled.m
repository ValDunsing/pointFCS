path= uigetdir; 
files=dir([path '/*final_fitparameters.txt']);
Ttaus=zeros(1,size(files,1));
Tfraction=Ttaus;
alphas=Ttaus;
Sfits=Ttaus;
tauds=Ttaus;
%tauerrors=zeros(1,size(files,1));
% namefile1=files(1).name;
% datafile1=load([path '/' namefile1]);
% ACFlength=size(datafile1,1);
% ACFvalues=zeros(ACFlength,size(files,1));
% timepoints=zeros(ACFlength,size(files,1));
% weightsdata=zeros(ACFlength,size(files,1));
meantau=0;
for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    tauds(i)=parametersi(2,1);
    Sfits(i)=parametersi(3,1);
    alphas(i)=parametersi(6,1);
    Tfraction(i)=parametersi(7,1);
    Ttaus(i)=parametersi(8,1);
    %tauerrors(i)=0.5*abs(parametersi(2,3)-parametersi(2,2));
    %meantau=meantau+taus(i)/tauerrors(i);
end
output=[Sfits' tauds'*10^6 Ttaus'.*10^6 alphas']
%meantau=meantau/sum(1./tauerrors);
% meantaufinal=taus;
% meantau=mean(meantaufinal);
% %meantaufinal(8)=[];
% meantaustd=std(meantaufinal);
% %meantauerror=sqrt(sum(tauerrors.^2));
% figure('Name','Transit times')
% errorbar(1:1:length(taus),taus,tauerrors,'xb','MarkerSize',8,'LineWidth',1.5)
% hold on
% errorbar(length(taus)+1,meantau,meantaustd,'sr','MarkerSize',8,'LineWidth',1.5)
% xlabel('Cell','FontSize',18)
% ylabel('\tau [s]','FontSize',18)
% 
% tauAlexa=1/3*(24.54+24.96+24.96)*10^-06;
% DAlexa=435;
% w0=sqrt(4*tauAlexa*DAlexa);
% 
% Ds=w0^2./(4*meantaufinal);
% Dmeanweighted=w0^2./(4*meantau);
% Dmean=mean(Ds);
% Dstd=std(Ds);