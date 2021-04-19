function [k, lambda_dark, lambda_grad, S] = blind_deconv_main(blur_B, k, ...
                                    lambda_dark, lambda_grad, threshold, opts)
% Do single-scale blind deconvolution using the input initializations
% 
% I and k. The cost function being minimized is: min_{I,k}
%  |B - I*k|^2  + \gamma*|k|_2 + lambda_dark*|I|_0 + lambda_grad*|\nabla I|_0
%
%% Input:
% @blur_B: input blurred image 
% @k: blur kernel
% @lambda_dark: the weight for the L0 regularization on intensity
% @lambda_grad: the weight for the L0 regularization on gradient
%
% Ouput:
% @k: estimated blur kernel 
% @S: intermediate latent image
%
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
%        Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 
%   [2] Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang,
%        Deblurring Text Images via L0-Regularized Intensity and Gradient
%        Prior, CVPR, 2014. 
%
%   Author: Jinshan Pan (sdluran@gmail.com)
%   Date  : 03/22/2016

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2013-08-11
H = size(blur_B,1);    W = size(blur_B,2);
blur_B_w = wrap_boundary_liu(blur_B, opt_fft_size([H W]+size(k)-1));
blur_B_tmp = blur_B_w(1:H,1:W,:);
Bx = conv2(blur_B_tmp, dx, 'valid');
By = conv2(blur_B_tmp, dy, 'valid');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
for iter = 1:opts.xk_iter
    %%
   if lambda_dark~=0
       S = L0Deblur_dark_chanel(blur_B_w, k, lambda_dark, lambda_grad, 2.0);
       S = S(1:H,1:W,:);
   else
       %% L0 deblurring
       S = L0Restoration(blur_B, k, lambda_grad, 2.0);
   end
   %% Necessary for refining gradient ???
   %
  [latent_x, latent_y, threshold]= threshold_pxpy_v1(S,max(size(k)),threshold); 
  %% 
%   latent_x = conv2(S, dx, 'valid');
%   latent_y = conv2(S, dy, 'valid');
  k_prev = k;
  %
  k = estimate_psf(Bx, By, latent_x, latent_y, 2, size(k_prev));
  %%
%   fprintf('pruning isolated noise in kernel...\n');
  CC = bwconncomp(k,8);
  for ii=1:CC.NumObjects
      currsum=sum(k(CC.PixelIdxList{ii}));
      if currsum<.1 
          k(CC.PixelIdxList{ii}) = 0;
      end
  end
  k(k<0) = 0;
  k(k(:)<3e-4) = 0;
  k=k/sum(k(:));
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Parameter updating
  if lambda_dark~=0;
      lambda_dark = max(lambda_dark/1.1, 1e-4);
  else
      lambda_dark = 0;
  end
  %lambda_dark = lambda_dark/1.1;  %% for natural images
  if lambda_grad~=0;
      lambda_grad = max(lambda_grad/1.1, 1e-4);
  else
      lambda_grad = 0;
  end
  %
%   figure(1); 
%   S(S<0) = 0;
%   S(S>1) = 1;
%   subplot(1,3,1); imshow(blur_B,[]); title('Blurred image');
%   subplot(1,3,2); imshow(S,[]);title('Interim latent image');
%   subplot(1,3,3); imshow(k,[]);title('Estimated kernel');
%   drawnow;
  %imwrite(S,'tmp.png')
%   kw = k - min(k(:));
%   kw = kw./max(kw(:));
%   imwrite(kw,'tmp_kernel.png')
%   mat_outname=sprintf('test3_blur_55_interim_kernel_new/interim_kernel_%d.mat',iter);
%   save(mat_outname,'k');
end;
k(k<0) = 0;  
k = k ./ sum(k(:));
end

function S = L0Deblur_dark_chanel(Im, kernel, lambda, wei_grad, kappa)
%%
% Image restoration with L0 regularized intensity and gradient prior
% The objective function:
% S = argmin ||I*k - B||^2 + \lambda |D(I)|_0 + wei_grad |\nabla I|_0
%% Input:
% @Im: Blurred image
% @kernel: blur kernel
% @lambda: weight for the L0 intensity prior
% @wei_grad: weight for the L0 gradient prior
% @kappa: Update ratio in the ADM
%% Output:
% @S: Latent image
%
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
%        Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 

if ~exist('kappa','var')
    kappa = 2.0;
end
%% pad image
% H = size(Im,1);    W = size(Im,2);
% Im = wrap_boundary_liu(Im, opt_fft_size([H W]+size(kernel)-1));
%%
S = Im;
betamax = 1e5;
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(Im);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
%%
KER = psf2otf(kernel,sizeI2D);
Den_KER = abs(KER).^2;
%%
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
if D>1
    Denormin2 = repmat(Denormin2,[1,1,D]);
    KER = repmat(KER,[1,1,D]);
    Den_KER = repmat(Den_KER,[1,1,D]);
end
Normin1 = conj(KER).*fft2(S);
%% pixel sub-problem
%%
dark_r = 35; %% Fixed size!
%mybeta_pixel = 2*lambda;
%[J, J_idx] = dark_channel(S, dark_r);
mybeta_pixel = lambda/(graythresh((S).^2));
maxbeta_pixel = 2^3;
while mybeta_pixel< maxbeta_pixel
    %% 
    [J, J_idx] = dark_channel(S, dark_r);
    u = J;
    if D==1
        t = u.^2<lambda/mybeta_pixel;
    else
        t = sum(u.^2,3)<lambda/mybeta_pixel;
        t = repmat(t,[1,1,D]);
    end
    u(t) = 0;
    %
    clear t;
    u = assign_dark_channel_to_pixel(S, u, J_idx, dark_r);
    %% Gradient sub-problem
    beta = 2*wei_grad;
    %beta = 0.01;
    while beta < betamax
        Denormin   = Den_KER + beta*Denormin2 + mybeta_pixel;
        %
        h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
        v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
        if D==1
            t = (h.^2+v.^2)<wei_grad/beta;
        else
            t = sum((h.^2+v.^2),3)<wei_grad/beta;
            t = repmat(t,[1,1,D]);
        end
        h(t)=0; v(t)=0;
        clear t;
        %
        Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
        Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
        %
        FS = (Normin1 + beta*fft2(Normin2) + mybeta_pixel*fft2(u))./Denormin;
        S = real(ifft2(FS));
        %%
        beta = beta*kappa;
        if wei_grad==0
            break;
        end
    end
    mybeta_pixel = mybeta_pixel*kappa;
end
%
end

function [J, J_index] = dark_channel(I, patch_size)
% function J = dark_channel(I, patch_size);

% Computes the "Dark Channel" of corresponding RGB image.
% -Finds from the input image the minimum value among all 
%  pixels within the patch centered around the location of the 
%  target pixel in the newly created dark channel image 'J'
%  J is a 2-D image (grayscale).

% Example: J = dark_channel(I, 15); % computes using 15x15 patch

% Check to see that the input is a color image
% if ndims(I) == 3
%     [M N C] = size(I);
%     J = zeros(M, N); % Create empty matrix for J
%     J_index = zeros(M, N); % Create empty index matrix
% else
%     error('Sorry, dark_channel supports only RGB images');
% end
%% for grayscale image
%
[M, N, C] = size(I);
J = zeros(M, N); % Create empty matrix for J
J_index = zeros(M, N); % Create empty index matrix

% Test if patch size has odd number
if ~mod(numel(patch_size),2) % if even number
    error('Invalid Patch Size: Only odd number sized patch supported.');
end

% pad original image
%I = padarray(I, [floor(patch_size./2) floor(patch_size./2)], 'symmetric');
I = padarray(I, [floor(patch_size./2) floor(patch_size./2)], 'replicate');

% Compute the dark channel 
for m = 1:M
        for n = 1:N
            patch = I(m:(m+patch_size-1), n:(n+patch_size-1),:);
            tmp = min(patch, [], 3);
            [tmp_val, tmp_idx] = min(tmp(:));
            J(m,n) = tmp_val;
            J_index(m,n) = tmp_idx;
        end
end

end

function [outImg] = assign_dark_channel_to_pixel(S, dark_channel_refine, dark_channel_index, patch_size)
%% assign dark channel value to image pixel
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
%        Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 
% 
[M N C] = size(S);
%outImge = zeros(M, N); % 

% pad original image
padsize = floor(patch_size./2);
S_padd = padarray(S, [padsize padsize], 'replicate');

% assign dark channel to pixel
for m = 1:M
        for n = 1:N
            patch = S_padd(m:(m+patch_size-1), n:(n+patch_size-1),:);
            if ~isequal(min(patch(:)), dark_channel_refine(m,n))
                patch(dark_channel_index(m,n)) = dark_channel_refine(m,n);
            end
            for cc = 1:C
                S_padd(m:(m+patch_size-1), n:(n+patch_size-1),cc) = patch(:,:,cc);
            end
        end
end


outImg = S_padd(padsize + 1: end - padsize, padsize + 1: end - padsize,:);
%% boundary processing
outImg(1:padsize,:,:) = S(1:padsize,:,:);  outImg(end-padsize+1:end,:,:) = S(end-padsize+1:end,:,:);
outImg(:,1:padsize,:) = S(:,1:padsize,:);  outImg(:,end-padsize+1:end,:) = S(:,end-padsize+1:end,:);

%figure(2); imshow([S, outImg],[]);
end




