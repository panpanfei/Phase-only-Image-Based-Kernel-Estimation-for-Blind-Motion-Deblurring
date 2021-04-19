%% citation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Code is created based on the method described in the following paper 
% %   [1] Pan, Liyuan, , Richard Hartley, Miaomiao Liu, Yuchao Dai. Phase-only image based kernel 
% %       estimation for single image blind deblurring. CVPR. 2019.
% %   [2] S. Cho, J. Wang, and S. Lee, Handling outliers in non-blind image deconvolution. 
% %       ICCV 2011.
% %   [3] Li Xu, Cewu Lu, Yi Xu, and Jiaya Jia. Image smoothing via l0 gradient minimization.
% %       ACM Trans. Graph., 30(6):174, 2011
% %   [4] Krishnan, Dilip, and Rob Fergus. Fast image deconvolution using hyper-Laplacian priors.
% %       Advances in Neural Information Processing Systems. 2009.
% %   [5] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
% %       Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 
% %   [5] Jinshan Pan, Zhe Hu, Zhixun Su, and Ming-Hsuan Yang,
% %       Deblurring Text Images via L0-Regularized Intensity and Gradient
% %       Prior, CVPR, 2014. 
% %   [6] Pan, Liyuan,  Yuchao Dai, Miaomiao Liu, and Fatih Porikli. Simultaneous stereo video
% %        deblurring and scene flow estimation. CVPR. 2017.
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Example:
% 1. input: Lenna.  para.needsys = 1; para.motion = 1;
% 2. input: post_blur.  para.needsys = 0; para.motion = 0; para.iter_num = 5;
%%
addpath('./code/');
addpath('./cho_code/');
load('para.mat')
%%
% If input is a clean image which need to be blurred
para.needsys = 1;
% the length of the motion kernel (pixel)
len     = 30;
% the dirction of motion kernel (cloclwise, degree)
theta   = 30;
% For linear kernel,motionb = 1, for non-linear kernel, motionb = 0
para.motion = 1;
% Fast result, may decrease quality
para.fast = 1;
% show figure or not
ifdisply = 0;
% iteration number
para.iter_num = 10;
%% kernel size
para.kernel_size = 35;
% half-kernel size
sn = floor(para.kernel_size/2);
% auto-corralation for show
auto_size = max(30,para.kernel_size);
%% read image and initial
blur_imagec = im2double(imread('./data/Lenna.png'));
% blur - grey im   blurc - color im
[blur, blurc]= data2blurim(blur_imagec, len, theta, para.needsys);
%% Auto-correlation
tic
[p_aut,text_aut,centrh,centrw ]= im2auto_corr(blur,auto_size,ifdisply);
% text_aut is the scaled cross-correlation map
%% find the bright peak point with is direction and length
[blurlen, bluranle] = auto2motion(text_aut);
%% coarse kernel
if para.motion==1
    %%linear motion 
    cPSF = fspecial('motion', blurlen, bluranle);
else
    %%complex motion
    sn = max(sn,floor(blurlen/2));
    cPSF = text_aut(centrh-sn:centrh+sn,centrw-sn:centrw+sn);
end
%% latant image and refined kernel
[PSF, deblurring] = kernelwithLatent(blurc, cPSF, para);
toc
%% result
figure,imshow([blurc,deblurring]);
imwrite(PSF,'./result/ker.png')
imwrite(deblurring,'./result/im.png')


