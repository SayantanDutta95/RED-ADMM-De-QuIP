% Sample code of the papers:
% 
% Sayantan Dutta, Nwigbo Kenule Tuador, Jérome Michetti, Bertrand Georgeot,
% Duong Hung Pham, Adrian Basarab, and Denis Kouamé, "Quantum denoising-based
% super-resolution algorithm applied to dental tomography images", 2022 IEEE 
% International Symposium on Biomedical Imaging (ISBI), 2022.
%
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé,
% "A Novel Image Denoising Algorithm Using Concepts of Quantum Many-Body Theory,"
% arXiv preprint arXiv(2021).
%
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé,
% "Image Denoising Inspired by Quantum Many-Body physics,"
% 2021 IEEE International Conference on Image Processing (ICIP), 2021,
% pp. 1619-1623, doi: 10.1109/ICIP42928.2021.9506794.
%
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé,
% "Plug-and-Play Quantum Adaptive Denoiser for Deconvolving Poisson Noisy Images",
% IEEE Access, 2021, vol. 9, pp. 139771-139791, doi: 10.1109/ACCESS.2021.3118608.
%
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé,
% "Quantum mechanics-based signal and image representation: Application
% to denoising," IEEE Open Journal of Signal Processing, vol. 2, pp. 190–206, 2021.
%
% Zhao, N., Wei, Q., Basarab, A., Dobigeon, N., Kouamé, D., & Tourneret, J. Y. (2016),
% "Fast Single Image Super-Resolution Using a New Analytical Solution for l_2-l_2 Problems",
% IEEE Transactions on Image Processing, 25(8), 3683-3697.
%
% One should cite all these papers for using the code.
%---------------------------------------------------------------------------------------------------
%
% MATLAB code prepard by Sayantan Dutta
% E-mail: sayantan.dutta@irit.fr and sayantan.dutta110@gmail.com
% Please contact Sayantan Dutta for the latest MATLAB code and any further communications.
%
%---------------------------------------------------------------------------------------------------
%
% The following script shows an example of a single image super-resolution
% algorithm (RED+ADMM+De-QuIP) using the quantum adaptive denoiser Denoising
% by Quantum Interactive Patches (De-QuIP) in a Regularization by Denoising (RED)
% approache with the ADMM framework
%
%---------------------------------------------------------------------------------------------------
%
% Download the Quantum adaptive denoiser "De-QuIP" from the following link:
% https://github.com/SayantanDutta95/De-QuIP-Denoising
% and save it in a file name 'De-QuIP' before running the code.
%
%---------------------------------------------------------------------------------------------------
%%
clc;
clear;
close all;

% configure the path
% quantum adaptive denoising functions De-QuIP
addpath(genpath('./De-QuIP/'));
% SD, FP, and ADMM methods
addpath(genpath('./minimizers/'));
% contains the default params
addpath(genpath('./parameters/'));
% contains basic functions
addpath(genpath('./helper_functions/'));
% test images for the debluring and super resolution problems, 
addpath(genpath('./test_images/'));

light_mode = false;
if light_mode
    fprintf('Running in light mode. ');
    fprintf('Tune the parameter properly for ADMM and De-QuIP and run the code in a sub optimal but faster mode.\n');
else
    fprintf('Light mode option is off. ');
    fprintf('Tune the parameter properly for the ADMM and for the quantum denoiser De-QuIP.\n');
end

%% read the original image

% choose the image and adjust the parameters for the ADMM and for the quantum denoiser De-QuIP
file_name = '11b_196.png'; 

fprintf('Reading %s image...', file_name);
orig_im = imread(['./test_images/' file_name]);
orig_im = double(orig_im);

% make the dimension even
orig_im = orig_im(1 : end - 1,1 : end - 1);


orig_im = 255.* orig_im./max(orig_im(:));

figure; imagesc(orig_im); colormap gray; axis on; colorbar
fprintf(' Done.\n');
pause(0.1);

%% define the degradation model

% choose the secenrio: 'GaussianBlur', or 'Downscale'
degradation_model = 'Downscale';

fprintf('Test case: %s degradation model.\n', degradation_model);

switch degradation_model
    
    case 'GaussianBlur'
        % filter size
        psf_sz = 9;
        % std of the Gaussian filter
        gaussian_std = 3;
        % create gaussian filter
        psf = fspecial('gaussian', psf_sz, gaussian_std);
        % use fft to solve a system of linear equations in closed form
        use_fft = true;
        % scaling factor
        scale = 2;
        % create a function handle to blur the image
        ForwardFunc = ...
            @(in_im) imfilter(in_im,psf,'conv','same','circular');
        % the psf is symmetric, i.e., the ForwardFunc and BackwardFunc
        % are the same
        BackwardFunc = ForwardFunc;
        % special initialization (e.g. the output of other method)
        % set to identity mapping
        InitEstFunc = @(in_im) in_im;
        
    case 'Downscale'
        % filter size
        psf_sz = 9;
        % std of the Gaussian filter
        gaussian_std = 3;
        % create gaussian filter
        psf = fspecial('gaussian', psf_sz, gaussian_std);
        % scaling factor
        scale = 2;
        
        % compute the size of the low-res image
        lr_im_sz = [ceil(size(orig_im,1)/scale),...
                    ceil(size(orig_im,2)/scale)];        
        % create the degradation operator
        H = CreateBlurAndDecimationOperator(scale,lr_im_sz,psf);
        % downscale
        ForwardFunc = @(in_im) reshape(H*in_im(:),lr_im_sz);        
        % upscale
        BackwardFunc = @(in_im) reshape(H'*in_im(:),scale*lr_im_sz);
        % special initialization (e.g. the output of other method)
        % use bicubic upscaler
        InitEstFunc = @(in_im) imresize(in_im,scale,'bicubic');
        
    otherwise
        error('Degradation model is not defined');
end


%% degrade the original image

% blur the clean image
[FB,FBC,F2B,input_im] = HXconv(orig_im,psf);
input_im = input_im(1:scale:end,1:scale:end,:);
            
randn('state', 0);

% add noise
fprintf(' Adding noise...');
BSNRdb =30; 
N = numel(input_im);
input_sigma = norm(input_im-mean(mean(input_im)),'fro')/sqrt(N*10^(BSNRdb/10));
input_im = input_im + input_sigma*randn(size(input_im));

% convert to YCbCr color space if needed
input_luma_im = PrepareImage(input_im);
orig_luma_im = PrepareImage(orig_im);

if strcmp(degradation_model,'Downscale')
    % upscale using bicubic
    input_im_bicubic = imresize(input_im,scale,'bicubic');
    input_im_bicubic = input_im_bicubic(1:size(orig_im,1), 1:size(orig_im,2), :); 
end
fprintf(' Done.\n');
psnr_input = ComputePSNR(orig_im, input_im_bicubic);


%% minimize the Laplacian regularization functional via ADMM

fprintf('Restoring using RED: ADMM method\n');

switch degradation_model
    case 'GaussianBlur'
        params_admm = GetGaussianDeblurADMMParams(light_mode, psf, use_fft);
    case 'Downscale'
        assert(exist('use_fft','var') == 0);
        params_admm = GetSuperResADMMParams(light_mode, psf);        
    otherwise
        error('Degradation model is not defined');
end
tic
[est_admm_im, psnr_admm, ssim_out, data] = RunADMM(input_luma_im,...
                                   ForwardFunc,...
                                   BackwardFunc,...
                                   InitEstFunc,...
                                   input_sigma,...
                                   params_admm,...
                                   orig_luma_im, input_im_bicubic);
out_admm_im = MergeChannels(input_im,est_admm_im);
toc

fprintf('Done.\n');

%% display final results

fprintf('Input PSNR = %f \n', psnr_input);
fprintf('RED: ADMM PSNR = %f \n', psnr_admm);

%% write images

output_Dir = './results/RED-De-QuIP';  % save HR image
if(exist(output_Dir, 'dir') == 0)
    mkdir(output_Dir);
end

figure;
subplot(131);imagesc(orig_im); colormap gray; axis on;
subplot(132);imagesc(input_im); colormap gray; axis on; title(sprintf('Noisy: \n  PSNR %.2f', ComputePSNR(orig_im, input_im_bicubic)));
subplot(133);imagesc(est_admm_im);colormap gray; axis on; title(sprintf('Noisy: \n  PSNR %.2f', ComputePSNR(orig_im, est_admm_im)));

% save image
images_Name = strcat(output_Dir, file_name(1 : end - 4), '_red_dequip_.png');
imwrite(est_admm_im, images_Name) 
% save data
save(output, 'orig_im', 'input_im', 'est_admm_im', 'data');

font_size = 18;
figure;
plot(data.psnr, 'k-', 'LineWidth', 2); hold on
ylabel('PSNR (dB)','Interpreter','latex','Fontsize',font_size);
xlabel('Iterations','Interpreter','latex','Fontsize',font_size);
legend({'RED-De-QuIP'},'Interpreter','latex', 'Location','southeast','Fontsize',font_size,'box','off');
axis([0 100 20 42])
set(gca,'FontSize',font_size);
