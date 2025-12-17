% FCS-calibrated imaging
% calculate mean intensities of 1xEGFP images in an 7x7 array (approx. width
% osberavtion volume in 72nm pixel-size images)
clear all;
close all;
fclose all; 

% get directory
% store original FCS and image files in two seperate folders
%path_FCS=uigetdir;
path_FCS='/Volumes/PortableSSD/Valentin/Dokumente/Postdoc/Data/Data_Claudio/Smog-GFP ScabRNAi 18-7-25/E3/FCS';
path_image='/Volumes/PortableSSD/Valentin/Dokumente/Postdoc/Data/Data_Claudio/Smog-GFP ScabRNAi 18-7-25/E3';
path_results='/Volumes/PortableSSD/Valentin/Dokumente/Postdoc/Data/Data_Claudio/Smog-GFP ScabRNAi 18-7-25/dist';
%path_image=uigetdir;
%path_results= uigetdir;

%results file, sample + pixel positions + intensity in kHz
filenameresults=sprintf('Invagination_pos_dist_6.txt'); %Does that also work for windows, or do we need \ there?
%fid1=fopen([path_results '\' filenameresults],'a');
fid1=fopen([path_results '/' filenameresults],'a');
fprintf(fid1,'sample\t x_FCS\t y_FCS\t x_invag\t y_invag\t dist\n'); % \n Zeileumbruch

files_FCS=dir([path_FCS '/Image*.czi']); % added "Image" by Valentin 2025-09-03
files_image=dir([path_image '/Image*.czi']); % added "Image" by Valentin 2025-09-03


% load image and image metadata
for i=1:size(files_FCS,1) % von 1 bis zum letzten Eintrag in Dimension 1 der Matrix
    i
    % load FCS metadata
    namefile_FCS=files_FCS(i).name
    %data2=bfopen_justinfo([path_FCS '\' namefile_FCS]); % test for windows 11
    data2=bfopen_metadata([path_FCS '/' namefile_FCS]);% old version 
    imagemetadata_FCS=data2{1,2};
    position_FCS_X=str2double(imagemetadata_FCS.get('Global Information|Image|S|Scene|Position|X #1'));
    position_FCS_Y=str2double(imagemetadata_FCS.get('Global Information|Image|S|Scene|Position|Y #1'));
    
    % load czi image metadata
    namefile_image=[namefile_FCS(1:6) num2str(str2double(namefile_FCS(7:end-4))-1) '.czi']
    %namefile_image=files_image(i).name
    %data=bfopen([path_image '\' namefile_image]);
    data=bfopen([path_image '/' namefile_image]);
    imagemetadata=data{1,2};
    pixelnumb=str2double(imagemetadata.get('Global Information|Image|SizeX #1'));
    pixelsize=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'))*1e6;  %in microM
    imagesize=pixelsize*pixelnumb; % in microM
    integrationtime=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|CameraIntegrationTime #1')); %in S
    position_image_X=str2double(imagemetadata.get('Global Information|Image|S|Scene|Position|X #1'));
    position_image_Y=str2double(imagemetadata.get('Global Information|Image|S|Scene|Position|Y #1'));     
    
    % Calculate x,y spot position in image 
    X_pxl=round(((position_FCS_X-position_image_X)/imagesize*pixelnumb)+(pixelnumb/2))
    Y_pxl=round((pixelnumb/2)-(position_FCS_Y-position_image_Y)/imagesize*pixelnumb)
    
    % Load tif image and set position of invagination
    
        %namefile_image_control=[namefile_FCS(1:6) num2str(str2double(namefile_FCS(7:end-4))-2) '.czi'] % adjustment for control images
        %data=bfopen([path_image '/' namefile_image_control]);
        %imageframe=data{1,1}{1,1};
    
    
    %imageframe=double(imread([path_image '\' namefile_image(1:end-4) '.tif' ])); %commented for control images
    imageframe=double(imread([path_image '/' namefile_image(1:end-4) '.tif' ])); %commented for control images

    
    %imageframe=imageframe'; %xy position according to Fiji
    
    % Mean filter:
    kernel = ones(3, 3) / 9; % 3x3 mean kernel
    imageframe_filt = conv2(imageframe, kernel, 'same'); % Convolve keeping size of I
    
    figure('Name','Image')
    h=imagesc(imageframe_filt)
    hold on
    plot(X_pxl,Y_pxl, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
    hold on
    [X_inv,Y_inv] = getpts;
    dist=(X_pxl-X_inv)*pixelsize; %lateral distance in um, we only care about the lateral distance
    hold on
    plot(X_inv(1),Y_inv(1), 'm+', 'MarkerSize', 30, 'LineWidth', 2);
    %saveas(h,[path_results '\' namefile_image(1:end-4) '_spots.fig'])
    saveas(h,[path_results '/' namefile_image(1:end-4) '_spots.fig'])
    
    
    fprintf(fid1,namefile_image(1:end-4));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\n',[X_pxl Y_pxl X_inv(1) Y_inv(1) dist]'); % %e Platzhalter, e wegen scientific
end



    