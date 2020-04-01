function N = WhiteNoise(b, level, seed)
%
%      N = WhiteNoise(b, level, seed);
%
%  This function generates Gaussian white noise for the 
%  data b. 
%
%  Input:  b - array containing data
%      level - scalar in [0, 1] specifiying level (percentage) of
%              noise.  For example, level = 0.01 implies
%                  norm(N)/norm(b) = 0.01, or 1% noise
%              Default is level = 0.01.
%       seed - Used to set the random number generator.
%              Default is seed = 0.
%
%  Output: N - array same dimension as b, containing pseudo-random
%              values drawn from a normal distribution with mean zero
%              and standard deviation one, and scaled as described above.
%

% Check inputs and set default values.
if nargin == 1, level = [];, seed = [];, end
if nargin == 2, seed = [];, end
if isempty(level), level = 0.01;, end
if isempty(seed), seed = 0;, end

% Generate noise.
randn('seed', seed);
N = randn(size(b));
N = N / norm(N(:));
N = level*norm(b(:))*N;
