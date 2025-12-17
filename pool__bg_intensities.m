pixeltime=1.23e-06;
path= uigetdir; 
files=dir([path '/*fluorescenceseries.mat']);
namefile=files(1).name;
datai=load([path '/' namefile]);

binningwindow_bg=100;

%backgroundintensities=zeros(size(files,1),length(datai.fluorescenceseries));
backgroundintensities=zeros(size(files,1),length(datai.fluorescenceseries_uncorr));

for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    datai=load([path '/' namefile]);
    %backgroundintensities(i,:)=datai.fluorescenceseries;
    backgroundintensities(i,:)=datai.fluorescenceseries_uncorr;
end

%backgroundintensities_binned=zeros(size(files,1),length(datai.fluorescenceseries)/binningwindow_bg);
backgroundintensities_binned=zeros(size(files,1),length(datai.fluorescenceseries_uncorr)/binningwindow_bg);
for i=1:size(files,1)
    backgroundintensity_i=backgroundintensities(i,:);
%    lastbinlength=mod(length(datai.fluorescenceseries),binningwindow_bg)+binningwindow_bg;
    lastbinlength=mod(length(datai.fluorescenceseries_uncorr),binningwindow_bg)+binningwindow_bg;
    backgroundintensity_binned=mean(reshape(backgroundintensity_i(1:end-lastbinlength),binningwindow_bg,length(backgroundintensity_i(1:end-lastbinlength))/binningwindow_bg),1);
    backgroundintensities_binned(i,:)=[backgroundintensity_binned mean(backgroundintensity_i(end-lastbinlength+1:end))];
end
backgroundintensity_avgbinned=mean(backgroundintensities_binned,1)-std(backgroundintensities_binned,1); %adjust and remove std  if no negative values occur
%backgroundintensity_stdbinned=std(backgroundintensities_binned,1);


% backgroundintensity_avg=mean(backgroundintensities,1);
% lastbinlength=mod(length(backgroundintensity_avg),binningwindow_bg)+binningwindow_bg;
% backgroundintensity_avgbinned=mean(reshape(backgroundintensity_avg(1:end-lastbinlength),binningwindow_bg,length(backgroundintensity_avg(1:end-lastbinlength))/binningwindow_bg),1);
% backgroundintensity_avgbinned=[backgroundintensity_avgbinned mean(backgroundintensity_avg(end-lastbinlength+1:end))];
timeline=1:1:size(backgroundintensities,2);
timeline=timeline'*pixeltime;
timeline_binned_bg=1:1:length(backgroundintensity_avgbinned);
timeline_binned_bg=(timeline_binned_bg-1)*binningwindow_bg*pixeltime;

h=figure('Name','Background Average');
plot(timeline_binned_bg,backgroundintensity_avgbinned)
hold on
%f = @(b,x) b(1).*exp(b(2).*x)+b(3).*exp(b(4).*x) + b(5);
%ft = fittype('a.*exp(b.*x)+c.*exp(d.*x)+e');
fbg=fit(timeline_binned_bg',backgroundintensity_avgbinned','exp2');
%fbg=fit(timeline_binned_bg',backgroundintensity_avgbinned',ft,'StartPoint',[0.2346,-1.1020,0.1662,-0.0186,0.05]);
plot(fbg)
saveas(h,[path '\bg_fit.fig'])
%correctionfit_bg=fbg.a.*exp(fbg.b*timeline)+fbg.c.*exp(fbg.d*timeline)+fbg.e;
correctionfit_bg=fbg.a.*exp(fbg.b*timeline)+fbg.c.*exp(fbg.d*timeline);
af=fbg.a;
bf=fbg.b;
cf=fbg.c;
df=fbg.d;
%ef=fbg.e;
%correctionfit_bg=fbg.b(1).*exp(fbg.b(2)*timeline)+fbg.b(3).*exp(fbg.b(4)*timeline)+b(5);

hh=figure('Name','Background Average cut');
index_bin_cut=length(timeline_binned_bg)/6+1;
index_cut=length(timeline)/6+1;
plot(timeline_binned_bg(index_bin_cut:end),backgroundintensity_avgbinned(index_bin_cut:end))
hold on
fbg=fit(timeline_binned_bg(index_bin_cut:end)',backgroundintensity_avgbinned(index_bin_cut:end)','exp2');
%fbg=fit(timeline_binned_bg(index_bin_cut:end)',backgroundintensity_avgbinned(index_bin_cut:end)',ft,'StartPoint',[af,bg,cf,df,0.1]);
plot(fbg)
saveas(hh,[path '\bg_fit_cut.fig'])
correctionfit_bg_cut=fbg.a.*exp(fbg.b*timeline(index_cut:end))+fbg.c.*exp(fbg.d*timeline(index_cut:end));%+fbg.e;
%correctionfit_bg_cut=fbg.a.*exp(fbg.b*timeline(index_cut:end))+fbg.c.*exp(fbg.d*timeline(index_cut:end))+fbg.e;
%fbg.e
save([path '\background'],'correctionfit_bg');
save([path '\background_cut'],'correctionfit_bg_cut');