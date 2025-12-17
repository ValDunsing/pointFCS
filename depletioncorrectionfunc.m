%global movingaveragewindowsize blocksize spatialfilter depletioncorrection...
%    intensityfilter intensitythreshold
function [linefluorescenceseries_corrected] =depletioncorrectionfunc(timeline,linefluorescenceseries,f)
% Depletion correction of time trace;
% Exponential
correctionfit=f.a.*exp(f.b*timeline)+f.c.*exp(f.d*timeline);
%correctionfit=f.a.*exp(f.ta*timeline)+f.b.*exp(f.tb*timeline)+f.c.*exp(f.tc*timeline);
linefluorescenceseries_corrected=linefluorescenceseries./sqrt(correctionfit./correctionfit(1))+correctionfit(1)*(1-sqrt(correctionfit./correctionfit(1)));
% Sinusoidal
%correctionfit=f.a1*sin(f.b1*timeline+f.c1)+f.a2*sin(f.b2*timeline+f.c2)+f.a3*sin(f.b3*timeline+f.c3)+f.a4*sin(f.b4*timeline+f.c4)+f.a5*sin(f.b5*timeline+f.c5)+f.a6*sin(f.b6*timeline+f.c6)+f.a7*sin(f.b7*timeline+f.c7)+f.a8*sin(f.b8*timeline+f.c8)+f.d1;
%linefluorescenceseries_corrected=linefluorescenceseries./correctionfit;

