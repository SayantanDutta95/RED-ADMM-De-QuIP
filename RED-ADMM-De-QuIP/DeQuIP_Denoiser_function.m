% Sample code of the papers:
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
% "Quantum mechanics-based signal and image representation: Application
% to denoising," IEEE Open Journal of Signal Processing, vol. 2, pp. 190–206, 2021.
% 
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé,
% "Despeckling Ultrasound Images Using Quantum Many-Body Physics,"
% 2021 IEEE International Ultrasonics Symposium (IUS), 2021, pp. 1-4,
% doi: 10.1109/IUS52206.2021.9593778.
%
% One should cite all these papers for using the code.
%---------------------------------------------------------------------------------------------------
% MATLAB code prepard by Sayantan Dutta
% E-mail: sayantan.dutta@irit.fr and sayantan.dutta110@gmail.com
% 
% This script shows an example of our image denoising algorithm 
% Denoising by Quantum Interactive Patches (De-QuIP)
%---------------------------------------------------------------------------------------------------


function [f_x_est] = DeQuIP_Denoiser_function(ima_nse, effective_sigma)


% ima_nse         = noisy image
% orig_ima        = clean image
% effective_sigma = noise level


% Choose parameters
WP=7;              % patch size
hW=10;              % half window size
factor_thr= 2.5;    % thresholding factor

del = effective_sigma;           
threshold = factor_thr * del;
delta = hW;         %< 2*hW+WP; % half window size for the searching zone

% Choose approximate parameter value

p = 0.1;    % proportionality constant
d =  10;    % reduced dimention
fact = 1.7; % Planck as a factor of height



% Start denoising process

% Divide image into small patches
ima_patchs = spatial_patchization(ima_nse, WP); % Divide noisy image

% Denoising process
ima_patchs_fil = DeQuIP_denoising(ima_patchs, hW, delta, p, d, fact, threshold, del);

% Reproject the small patches to construct the denoised image
ima_fil_DeQuIP = reprojection_UWA(ima_patchs_fil);

f_x_est = ima_fil_DeQuIP;

end
