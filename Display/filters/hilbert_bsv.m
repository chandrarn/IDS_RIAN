function x = hilbert_bsv(h_sig)

% this calculates the Hilbert transform
% adapted from a code written by SDK
% this only works if the input is a row vector

xfft = fft(h_sig);

n = length(xfft);
h = zeros(1, n);
n2 = n/2;
n2r = round(n/2);
nc = n2/n2r;

if nc == 1
    h(1) = 1;
    h(n2 + 1) = 1;
    h(2: n2) = 2;
elseif nc ~= 1
    h(1) = 1;
    h(2: n2r - 1) = 2;
end

xprod = xfft.*h;

x = -imag(ifft(xprod));







