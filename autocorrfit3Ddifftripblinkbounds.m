function [N,taud,Sfit,T1fraction,T1tau,T2fraction,T2tau,alpha,CI,fitcurve,residuals] = autocorrfit3Ddifftripblinkbounds(tcorr,fcorr,lb,ub,weights,fitfunc,x0,fixed)
%[fcorrfitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,tcorr,fcorr,fitfunc,x0,'Weights',weights);
weightsvector=sqrt(abs(weights));
lsfitfunc=@(x) weightsvector.*(fitfunc(x,tcorr)-fcorr);
[fcorrfitparameters,fcorrfitresnorm,residuals,fcorrfitexitflag,fcorrfitoutput,lambda,jacobian] = lsqnonlin(lsfitfunc,x0,lb,ub);
MSE=(residuals*residuals')/(length(fcorr)-length(fcorrfitparameters));
% COVB=inv(jacobian'*jacobian)*MSE;
% size(COVB)
N=real(fcorrfitparameters(1));
taud=real(fcorrfitparameters(2));
Sfit=real(fcorrfitparameters(3));
T1fraction=fcorrfitparameters(4);
T1tau=fcorrfitparameters(5);
T2fraction=fcorrfitparameters(6);
T2tau=fcorrfitparameters(7);
alpha=fcorrfitparameters(8);

% Compute confidence intervals of coefficients
%CI=0;
fitparametersnotfixed=fcorrfitparameters(fixed==false);
%try 
    %CInotfixed = nlparci(fitparametersnotfixed,residuals,'covar',COVB);
    %CInotfixed = nlparci(fitparametersnotfixed,residuals,'jacobian',jacobian);
    CInotfixed = nlparci(fcorrfitparameters,residuals,'jacobian',jacobian);
%catch
%    CInotfixed=zeros(length(fcorrfitparameters),2);
%end
% CInotfixed = nlparci(fitparametersnotfixed,residuals,'covar',COVB);

CI=zeros(length(fcorrfitparameters),2);
k=1;
for i=1:length(fcorrfitparameters)
    %if fixed(i)==false
        CI(i,:)=CInotfixed(k,:);
        k=k+1;
    %else
     %   CI(i,:)=[0 0];
    %end
end    
fitcurve=((1-T1fraction+T1fraction.*exp(-tcorr./T1tau))./(1-T1fraction)).*((1-T2fraction+T2fraction.*exp(-tcorr./T2tau))./(1-T2fraction)).*1./N.*((1+(tcorr./taud).^alpha).^-1).*(1+(tcorr./taud).^alpha.*Sfit.^-2).^-0.5;
%[fcorrfitparameters,fcorrfitresnorm,fcorrfitresiduals,fcorrfitexitflag,fcorrfitoutput]