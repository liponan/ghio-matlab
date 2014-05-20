% 2-D HIO written by Po-Nan Li @ Academia Sinica 2012
function [R, Sup] = shrinkwrap(Fabs, S, n, gen, n2, varargin) % Fabs, S, n, unknown, alpha

% pre-allocated spaces
R   = zeros( size(Fabs, 1), size(Fabs, 2), gen+1);
Sup = false( size(Fabs, 1), size(Fabs, 2), gen+1);
Sup(:,:,1) = S;

% OSS
if empty(varargin)
    alpha = [];
else
    alpha = varargin{1};
end

% first run with given support
R(:,:,1) = hio2d(Fabs, S, n, alpha);

% make Gaussian kernel
fwhm = 3;
sig = fwhm / 2.355;
x = 1:size(S,2);
y = 1:size(S,1);
[X, Y] = meshgrid(x./(size(S,2)/2), y./(size(S,1)/2));
R = sqrt(X.^2 + Y.^2);


% shrink-wrap
for g = 1:gen
    G = exp(-(R./sqrt(2)./sig).^2);
    Sup(:,:,g+1) = ifft2( fft2(R(:,:,g)) .* fft2(G), 'symmetric');
    Sup(:,:,g+1) = ( Sup(:,:,g+1) >= 0.04*max(max(R(:,:,g))) );
    R(:,:,g+1) = hio2d(fft2(R(:,:,g)), Sup(:,:,g+1), n2, alpha);
    sig = sig * 0.99;
end
    