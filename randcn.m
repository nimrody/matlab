function n = randcn(varargin)

% n = randcn(r, c, d, ...)
% 
%   Generates an array r x c x d x ..  of complex random numbers (iid)
%   with unit power. If only one argument is supplied, generates
%   a random matrix r x r.


n = 1/sqrt(2)*(randn(varargin{:}) + 1j*randn(varargin{:}));


