function y = smooth_selfdefine(x,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   smooth variable y 
%%%   convolute the vaiable y with a causal Gaussian filter 
%%%  For example:
%%%  y = smooth_selfdefine(x)
%%%  y = smooth_selfdefine(x 'conv_step', 10, 'len_fix', true) % the width
%%%  of the Gaussian filter is 10, and length of y is same as x
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% default setting
conv_step = 10; 
len_fix = true;

if ~isempty(varargin)
    conv_step_index = find(strcmpi(varargin,'conv_step'));
    if ~isempty(conv_step_index)
        conv_step = varargin{conv_step_index+1};
    end
    len_fix_index = find(strcmpi(varargin,'len_fix'));
    if ~isempty(len_fix_index)
        len_fix = varargin{len_fix_index+1};
    end
end

% normalized causal Gaussian kernal 
kernal = gaussmf(0:1:conv_step,[conv_step 0]);
norm_kernal = sum(kernal);

% convolution
yy = conv(x,kernal/norm_kernal);

% keep the length as the same or not
if len_fix
    y = yy(1,1:size(x,2));
else
    y = yy;
end

