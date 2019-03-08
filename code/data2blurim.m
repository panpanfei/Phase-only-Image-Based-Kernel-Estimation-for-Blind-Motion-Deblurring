function [blur, blurc]= data2blurim(sharp_imagec, len, theta, motionb)

if size(sharp_imagec,3)~=1
    sharp_image=rgb2gray(sharp_imagec);
else
    sharp_image=(sharp_imagec);
end

if motionb == 1
k_in = fspecial('motion', len, theta);
 
blur = imfilter(sharp_image,k_in,'conv','symmetric');%conv2(sharp_image,k_in,'valid');%
blurred_im = blur;
%
% kernel_in = kernel2Matrix1(h,w,k_in);
% 
% tmp_Kt_B = kernel_in* reshape(sharp_image', [h*w 1]);
% blur = reshape(tmp_Kt_B, [w h])';
% 
% blurred_im = blur;

% [kh,kw] = size(k_in);
% halfkh = (kh-1)/2;
% halfkw = (kw-1)/2;
% 
% blurred_im = blur(halfkh+1:end-halfkh,halfkw+1:end-halfkw);
%
% sharp_image_pad = wrap_boundary_liu(sharp_image, opt_fft_size([h w]+size(k_in)-1));
% sizey = size(sharp_image_pad);
% otfk  = psf2otf(k_in, sizey); 
% Fyout = (otfk).*fft2(sharp_image_pad);%conj
% blur = real(ifft2(Fyout));
% blur = blur(1:h, 1:w, :);
% blurred_im = blur;
%

if size(sharp_imagec,3)~=1
    blurc =  blur;
    blurc(:,:,1)= imfilter(sharp_imagec(:,:,1),k_in,'conv','symmetric');% conv2(sharp_imagec(:,:,1),k_in,'valid');%
    blurc(:,:,2)= imfilter(sharp_imagec(:,:,2),k_in,'conv','symmetric');% conv2(sharp_imagec(:,:,1),k_in,'valid');%
    blurc(:,:,3)= imfilter(sharp_imagec(:,:,3),k_in,'conv','symmetric');% conv2(sharp_imagec(:,:,1),k_in,'valid');%
else 
    blur = blur;
    blurc = blur;
end
else
% No blur
blur = sharp_image;
blurc = sharp_imagec;
% im_name = sprintf('./data//kitti/image/%06d_10.png',2);
% blurred_im=im2double(rgb2gray(imread(im_name)));% ring line circle
end
