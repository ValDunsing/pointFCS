function [N,taud,Sfit,Tfraction,Ttau,alpha,CI,fitcurve,residuals] = autocorrfit3Ddifftrip(tcorr,fcorr,lb,ub,weights,fitfunc,x0,fixed)
[fcorrfitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,tcorr,fcorr,fitfunc,x0,'Weights',weights);
%[fcorrfitparameters,fcorrfitresiduals] = lsqcurvefit(lsfitfunc,x0,tcorrfit,fcorrfit,lb,ub,'Weights',weights);
N=real(fcorrfitparameters(1));
taud=real(fcorrfitparameters(2));
Sfit=real(fcorrfitparameters(3));
Tfraction=fcorrfitparameters(4);
Ttau=fcorrfitparameters(5);
alpha=fcorrfitparameters(6);

% Compute confidence intervals of coefficients
%CI=0;
fitparametersnotfixed=fcorrfitparameters(fixed==false);
try 
    CInotfixed = nlparci(fitparametersnotfixed,residuals,'covar',COVB);
catch
    CInotfixed=zeros(length(fcorrfitparameters),2);
end
% CInotfixed = nlparci(fitparametersnotfixed,residuals,'covar',COVB);

CI=zeros(length(fcorrfitparameters),2);
k=1;
for i=1:length(fcorrfitparameters)
    if fixed(i)==false
        CI(i,:)=CInotfixed(k,:);
        k=k+1;
    else
        CI(i,:)=[0 0];
    end
end    
fitcurve=((1-Tfraction+Tfraction.*exp(-tcorr./Ttau))./(1-Tfraction)).*1./N.*((1+(tcorr./taud).^alpha).^-1).*(1+(tcorr./taud).^alpha.*Sfit.^-2).^-0.5;

