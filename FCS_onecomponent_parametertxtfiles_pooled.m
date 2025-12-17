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
Tfrac=Ints;
Ttau=Ints;
I0=Ints;

for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    Ns(i)=parametersi(1,1);
    Ints(i)=parametersi(5,1);
    taus(i)=parametersi(2,1);
    bleaching_fractions(i)=parametersi(4,1);
    alphas(i)=parametersi(6,1);
    Tfrac(i)=parametersi(7,1);
    Ttau(i)=parametersi(8,1);
    I0(i)=parametersi(9,1);
    Ss(i)=parametersi(3,1);
end
Bs=Ints./Ns;
output=[Ns Ints taus bleaching_fractions Bs alphas Tfrac Ttau I0];

path2= uigetdir;

filename=sprintf('/SMOG_GFP_lp1.txt');
%filename=sprintf('/2019-09-11_Rpl3.txt'); %Does that also work for windows, or do we need \ there?

fid1=fopen([path2  filename],'a'); % adjust path if necessary!
fprintf(fid1,'sample\t N\t I\t tau\t bleaching\t B\t alpha \t T\t tauT\t I0 \n');
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;