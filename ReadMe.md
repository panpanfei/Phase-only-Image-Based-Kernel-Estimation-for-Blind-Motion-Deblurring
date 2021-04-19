# For academic use only.
### This package contains an implementation of the image deblurring algorithm described in the paper: 

@inproceedings{pan2019phase,  
  title={Phase-only image based kernel estimation for single image blind deblurring},         
  author={Pan, Liyuan and Hartley, Richard and Liu, Miaomiao and Dai, Yuchao},    
  booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},      
  pages={6034--6043},     
  year={2019}     
} 

### If you use our code, please also cite 
  [1] Li Xu, Cewu Lu, Yi Xu, and Jiaya Jia. Image smoothing via l0 gradient minimization. ACM Trans. Graph., 30(6):174, 2011     
  
  [2] S. Cho, J. Wang, and S. Lee, Handling outliers in non-blind image deconvolution. ICCV 2011.           
  
  [3] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang, Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016.      
  
  [4] Pan, Liyuan,  Yuchao Dai, Miaomiao Liu, and Fatih Porikli. Simultaneous stereo video deblurring and scene flow estimation.          
      CVPR. 2017.     
      
More details can be found in our paper.  

### Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) in an academic publication.

Notes 
----------------
For algorithmic details, please refer to our paper.     

We give the result of our method in file 'result'. Include dataset from 'Levin', 'Kohler', 'Gong' and 'Pan'. For better result, some          
parameters are carefully designed. We also use some strategies like the 'patch-wise process' and 'coarse to fine', etc.         
The current release version is for tackle uniform motion blur. The non-uniform deblurring can be achieved by segment the input blur           
image to overlapped square patches.    


How to use
----------------
The code is tested in MATLAB 2015b(64bit) under the ubuntu 14.04 LTS 64bit version with an Intel Core i7-4790 CPU and 6 GB RAM.

1. unpack the package.      
2. include code/Phase_for_public/ in your Matlab path.      
3. Run "main_uniform.m" to try the examples included in this package.     

User-specified parameter:
----------------
### There are a few parameters that need to be specified by users.

'needsys'    :   1 for synthetic testing data. Blurred the image with a given kernel.     
'motion'     :   1 for linear kernal.       
'fast'       :   1 for fast processing strategy without coarse-to-fine.       
'kernel_size':   The size of blur kernel.       
'auto_size'  :   The scale for autocorrelation.           
'iter_num'   :   Iteration number.          
'lambda_grad' & lambda_l0 & lambda_tv: the weight for the L0/TV regularization.         

IMPORTANT NOTE 
----------------
1. Note that the algorithm sometimes may converge to an incorrect result. When you obtain such an incorrect result, please re-try to  
deblur with slightly changed parameters (e.g., using large blur kernel sizes, fast, or iter_num).  

2. Should you have any questions regarding this code and the corresponding results, please contact Liyuan.Pan@anu.edu.au

3. For the non-uniform blur, please refer to our paper and use a multi-patch strategy.
For example, we segment the image into small patches in size 80*80 with an overlapping of 30 pixels.

### Notes for reproduce. 
1. Following the equitions in our papar
2. Coded some process (e.g., fft or corralation) by yourself instead of using inbuid functions.
