function [p_aut,text_aut,a,b] = im2auto_corr(blur,auto_size)
%% FFT

% blurred_im_f = padarray(blurred_im,1*[h,w],'symmetric');
% 
% ff=fft2(blurred_im_f);      % Fourier transform of the image
% ff = ff(h1+1:2*h1,w1+1:2*w1);

ff=fft2(blur);
fmag=abs(ff);    % Fourier Magnitude
fphase = ff ./ (fmag);   % Fourier Phase  + 0.00001

fmag(1, 1) = 0;

% Scale to display

dfmag = mat2gray(fmag);
dfmag = fftshift(dfmag);
mag = 10;
dfmag = log(1+mag*dfmag);  %  enhance visulization


% Images reconstructed from the magnitude and phase
m_im=fftshift(ifftn(fmag));
p_im=real(ifftn(fphase));
p_im = abs(p_im);
p_im_new = p_im(3:end-2,3:end-2);
% p_im = mat2gray(p_im);
% p_im_new = padarray(p_im_new,[2,2]);
% p_im (p_im<0)= 0;
dp_im = p_im_new;             %  enhance visulization
% dp_im(dp_im>0.01) = 0;%0.01``````````````
dp_im = mat2gray(dp_im);

%% autocorrelation
f_paut = fft2(p_im_new);
f_paut2 = real(f_paut .* conj(f_paut));
p_aut = real(ifft2(f_paut2));
p_aut = fftshift(stretch(p_aut));
% p_aut = xcorr2(p_im_new,p_im_new);

% Take the autocorrelation of the origina (blurred) image
% f_baut = fft2(blurred_im);
% f_baut2 = real(f_baut .* conj(f_baut));
% baut = real(ifft2(f_baut2));
% baut(1, 1) = 0;

%--------------------------------------------------------------
%% detect direction and length
% add a filter before it to avoide noise
Sobfh = [-1 0 1;-2 0 2;-1 0 1];
Sobfv = Sobfh';
% Sobfh = [-1 -2 0 2 1;-2 -4 0 4 2;-1 -2 0 2 1]./32;
% Sobfv = Sobfh';
% p_im(p_im>0.01)=0;
p_eh = imfilter(p_im_new,Sobfh);
p_ev = imfilter(p_im_new,Sobfv);
p_edge = sqrt(p_eh.*p_eh+p_ev.*p_ev);
dp_edge = p_edge;             %  enhance visulization
dp_edge(dp_edge>0.02) = 0;
dp_edge = mat2gray(dp_edge);
% dp_edge = mat2gray(p_edge);


e_paut = fft2(p_edge);
e_paut2 = real(e_paut .* conj(e_paut));
e_aut = real(ifft2(e_paut2));
e_aut = fftshift(stretch(e_aut));

%% give more general visulize result
text_aut = p_aut;%p_aut;m_in
[counts,binLocations] = imhist(text_aut);
idx = find((counts==0)|(binLocations==0));
counts(idx) = [];
binLocations(idx) = [];
[mval,mloc] = max(counts);
counts(1:mloc) =[];
binLocations(1:mloc)  = [];
textv = [counts,binLocations];
% for choose scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%auto_s = min (floor(length(textv)./2),auto_size);
auto_s = min (floor(length(textv)-1),auto_size);
idxtv = find(counts<=50);
if size(idxtv,1)~=1
        text_aut((text_aut <= textv(idxtv(1)-1,2))) = 0;
end
text_aut((text_aut < textv(end-auto_s,2))) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b] = find(text_aut==max(text_aut(:)));

aa = (text_aut==max(text_aut(:)));
text_aut(aa)=0;
text_aut = mat2gray(text_aut);
text_aut(aa)=1;
end
