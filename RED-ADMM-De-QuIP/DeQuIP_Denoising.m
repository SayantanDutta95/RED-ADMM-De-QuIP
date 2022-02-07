clc
clear 
close all
addpath testsets/
addpath De-QuIP/


folderTest  = fullfile('testsets','Kodak24_gray'); %%% test dataset
%folderTest  = 'testsets\BSD68'; Set12  Urban100_gray  BSD68  LIVE1_gray  Kodak24_gray  new_test_all
output_Dir = 'D:\De-QuIP_denoisr\result_DeQuIP\Kodak24_gray\';  % save denoised image
if(exist(output_Dir, 'dir') == 0)
    mkdir(output_Dir);
end

%%% read images
ext         =  {'*.jpg','*.png','*.bmp'};
filePaths   =  [];
for i = 1 : length(ext)
    filePaths = cat(1,filePaths, dir(fullfile(folderTest,ext{i})));
end

noiseSigma  = 100;  %%% image noise level
Gaussian    = 1;
showResult  = 1;
pauseTime   = 1;
%%% PSNR and SSIM
PSNRs_in = zeros(1,length(filePaths));
SSIMs_in = zeros(1,length(filePaths));

PSNRs = zeros(1,length(filePaths));
SSIMs = zeros(1,length(filePaths));

for i = 1:length(filePaths)
    
        %%% read images
    label = imread(fullfile(folderTest,filePaths(i).name));
    imageName = filePaths(i).name;    % read image name
    [~,nameCur,extCur] = fileparts(filePaths(i).name);
    label = im2double(label);
    label = 255* label;

  %  randn('seed',0);
    input = label + noiseSigma*randn(size(label));  % noise Gaussian
    
    
    tstart = tic;
    output = DeQuIP_Denoiser_function(input, noiseSigma, Gaussian);
    telapsed = toc(tstart);


    % normalized
    label = label/255;
    input = input/255;
    output = output/255;
    %%% calculate PSNR and SSIM
    [PSNRin, SSIMin] = Cal_PSNRSSIM(im2uint8(label),im2uint8(input),0,0);

    [PSNRCur, SSIMCur] = Cal_PSNRSSIM(im2uint8(label),im2uint8(output),0,0);
    if showResult
        imshow(cat(2,im2uint8(label),im2uint8(input),im2uint8(output)));
        title([filePaths(i).name ,'    ',num2str(PSNRCur,'%2.2f'),'dB','    ',num2str(SSIMCur,'%2.4f')])
        drawnow;
        pause(pauseTime)
    end
    PSNRs_in(i) = PSNRin;
    SSIMs_in(i) = SSIMin;
    
    PSNRs(i) = PSNRCur;
    SSIMs(i) = SSIMCur;
    
    % save gray image
    images_Name = strcat(output_Dir, imageName(1 : end - 4), '_Gaus_denoised_dequip_',...
        num2str(noiseSigma),'.png');
    imwrite(output, images_Name)
end


disp([mean(PSNRs_in),mean(SSIMs_in)]);
disp([mean(PSNRs),mean(SSIMs)]);

% save data
PSNRs_in_mean = mean(PSNRs_in);
SSIMs_in_mean = mean(SSIMs_in);

PSNRs_mean = mean(PSNRs);
SSIMs_mean = mean(SSIMs);
data_file_Name = strcat(output_Dir, 'results_sigma_', num2str(noiseSigma),'.mat');
save(data_file_Name, 'PSNRs_mean', 'SSIMs_mean', 'PSNRs', 'SSIMs','PSNRs_in_mean','SSIMs_in_mean')

