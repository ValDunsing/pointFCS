clear all
close all
%filename N Int tau bleaching brightness
path= uigetdir; 
files=dir([path '/*final_fitparameters.txt']);
Ints=zeros(size(files,1),1);
Ns=Ints;
taus=Ints;
bleaching_fractions=Ints;
alphas=Ints;
Ss=Ints;
I0=Ints;
T1=Ints;
T2=Ints;
Ttau1=Ints;
Ttau2=Ints;

for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    Ns(i)=parametersi(1,1);
    Ints(i)=parametersi(10,1);
    taus(i)=parametersi(2,1);
    Ttau1(i)=parametersi(5,1);
    Ttau2(i)=parametersi(7,1);
    bleaching_fractions(i)=parametersi(9,1);
    alphas(i)=parametersi(8,1);
    T1(i)=parametersi(4,1);
    T2(i)=parametersi(6,1);
    I0(i)=parametersi(11,1);
    Ss(i)=parametersi(3,1);   
end
Bs=Ints./Ns;
output=[Ns Ints taus bleaching_fractions Bs alphas T1 Ttau1 T2 Ttau2 I0 Ss];

path2= uigetdir;

filename=sprintf('/2019-06-18_Rpl3-GFP_blink_50us.txt'); %Does that also work for windows, or do we need \ there?

fid1=fopen([path2  filename],'a'); % adjust path if necessary!
fprintf(fid1,'sample\t N\t I\t tau\t bleaching\t B\t alpha\t T1\t tauT1\t T2\t tauT2\t I0\t S\n');
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;