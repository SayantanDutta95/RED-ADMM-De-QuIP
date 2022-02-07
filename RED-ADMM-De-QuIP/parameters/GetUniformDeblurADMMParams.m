% Copyright 2017 Google Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     https://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Returns set of parameters
% set light_mode = true to run the code in a sub optimal but faster mode
% set light_mode = false to obtain the results reported in the RED paper
% psf - the Point Spread Function, used by the solver if use_fft == true
%       this feature is supported only for deblurring

function params = GetUniformDeblurADMMParams(light_mode, psf, use_fft)

params = GetUniformDeblurFPParams(light_mode, psf, use_fft);

% admm parameter
params.beta = 1e-3;

% number of denoising applications
params.inner_denoiser_iters = 1;

% relaxation parameter of ADMM
params.alpha = 2;


return


