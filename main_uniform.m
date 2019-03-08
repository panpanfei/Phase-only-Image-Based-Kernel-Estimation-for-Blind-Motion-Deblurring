%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Works better for uniform linear mothion blur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbstop if error
addpath('./code/');
load('para.mat')
%%
% linear motion for test (synthetic), motionb=1
% for real blur image, motionb=0
motionb = 1;
% the length of the motion kernel (pixel)
len     = 30;
% the dirction of motion kernel (cloclwise, degree)
theta   = 10;
% kernel size
kernel_size = 35;
% half-kernel size
sn = floor(kernel_size/2);
% auto-corralation for show
auto_size = 30;
%% read image and initial
tic
blur_imagec = im2double(imread('./data/Lenna.png'));

[blur, blurc]= data2blurim(blur_imagec, len, theta, motionb);

%% Auto-correlation
[p_aut,text_aut,centrh,centrw ]= im2auto_corr(blur,auto_size);
% text_aut is the scaled cross-correlation map
%% find the bright peak point with is direction and length


%% coarse kernel
if motionb==1
%%linear motion 
    [blurlen, bluranle] = auto2motion(text_aut);
    cPSF = fspecial('motion', blurlen, bluranle);
else
%%complex motion
    cPSF = text_aut(centrh-sn:centrh+sn,centrw-sn:centrw+sn);
end
%% latant image and its kernel
[PSF, deblurring] = kernelwithLatent(blurc, cPSF, para);
toc
%% color result
if motionb~=1
    % smooth filter for noise, can be removed
    deblurring = bilateral_filter(deblurring, 1, 0.05);
end
figure,imshow([blurc,deblurring]);
imwrite(PSF,'./result/ker.png')
imwrite(deblurring,'./result/im.png')


