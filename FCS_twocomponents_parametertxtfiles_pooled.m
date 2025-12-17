clear all
close all
%filename N Int tau bleaching brightness
path= uigetdir; 
files=dir([path '/*final_fitparameters.txt']);
Ints=zeros(size(files,1),1);
Ns=Ints;
taus1=Ints;
taus2=Ints;
bleaching_fractions=Ints;
alphas1=Ints;
alphas2=Ints;
Ss=Ints;
ffrees=Ints;
I0=Ints;
T=Ints;
Ttau=Ints;

for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    Ns(i)=parametersi(1,1);
    Ints(i)=parametersi(5,1);
    ffrees(i)=parametersi(2,1);
    taus1(i)=parametersi(3,1);
    taus2(i)=parametersi(4,1);
    bleaching_fractions(i)=parametersi(10,1);
    alphas1(i)=parametersi(6,1);
    alphas2(i)=parametersi(7,1);
    T(i)=parametersi(8,1);
    Ttau(i)=parametersi(9,1);
    Ints(i)=parametersi(11,1);
    I0(i)=parametersi(12,1);
    Ss(i)=parametersi(5,1);   
end
Bs=Ints./Ns;
output=[Ns ffrees taus1 taus2 alphas1 alphas2 T Ttau Ints I0 Bs bleaching_fractions Ss];

path2= uigetdir;

%filename=sprintf('/2019-09-11_rpl3_2_comp.txt'); %Does that also work for windows, or do we need \ there?
filename=sprintf('/2025-07-28_SMOG_scab_RNAi_bg_std.txt'); %Does that also work for windows, or do we need \ there?
%filename=sprintf('/2024-02-29_yw_dextran647_e4.txt'); %Does that also work for windows, or do we need \ there?

fid1=fopen([path2  filename],'a'); % adjust path if necessary!
fprintf(fid1,'sample\t N\t f free\t tau free\t tau bound\t alpha free\t alpha bound\t T\t tauT\t Imean\t I0\t B\t bleaching\t S\n');
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;