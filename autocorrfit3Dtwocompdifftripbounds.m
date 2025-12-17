function [N,f1,taud1,taud2,Sfit,Tfraction,Ttau,alpha1,alpha2,CI,fitcurve,residuals] = autocorrfit3Dtwocompdifftripbounds(tcorr,fcorr,lb,ub,weights,fitfunc,x0,fixed)
%[fcorrfitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,tcorr,fcorr,fitfunc,x0,'Weights',weights);
weightsvector=sqrt(abs(weights));
lsfitfunc=@(x) weightsvector.*(fitfunc(x,tcorr)-fcorr);
[fcorrfitparameters,fcorrfitresnorm,residuals,fcorrfitexitflag,fcorrfitoutput,lambda,jacobian] = lsqnonlin(lsfitfunc,x0,lb,ub);
MSE=(residuals*residuals')/(length(fcorr)-length(fcorrfitparameters));
COVB=inv(jacobian'*jacobian)*MSE;
size(COVB)
N=real(fcorrfitparameters(1));
f1=real(fcorrfitparameters(2));
taud1=real(fcorrfitparameters(3));
taud2=real(fcorrfitparameters(4));
Sfit=real(fcorrfitparameters(5));
Tfraction=fcorrfitparameters(6);
Ttau=fcorrfitparameters(7);
alpha1=fcorrfitparameters(8);
alpha2=fcorrfitparameters(9);

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
fitcurve=((1-Tfraction+Tfraction.*exp(-tcorr./Ttau))./(1-Tfraction)).*1./N.*(f1.*((1+(tcorr./taud1).^alpha1).^-1).*(1+(tcorr./taud1).^alpha1.*Sfit.^-2).^-0.5+(1-f1).*((1+(tcorr./taud2).^alpha2).^-1).*(1+(tcorr./taud2).^alpha2.*Sfit.^-2).^-0.5);
%[fcorrfitparameters,fcorrfitresnorm,fcorrfitresiduals,fcorrfitexitflag,fcorrfitoutput]