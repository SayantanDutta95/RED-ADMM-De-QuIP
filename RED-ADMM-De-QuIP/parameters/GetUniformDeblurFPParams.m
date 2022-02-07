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

function params = GetUniformDeblurFPParams(light_mode, psf, use_fft)

params = GetUniformDeblurSDParams(light_mode);

% number of outer iterations
if light_mode
    params.outer_iters = 100;
else
    params.outer_iters = 200;
end

% use fft for solving a linear system in a closed form
params.use_fft = use_fft;
params.psf = psf;

% number of inner iterations, set to 'nan' if use_fft == true
if use_fft
    params.inner_iters = nan;
else
    params.inner_iters = 200;
end

return





