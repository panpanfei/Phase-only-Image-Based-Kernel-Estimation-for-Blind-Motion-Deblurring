# For academic use only
This contains an implementation of the image deblurring algorithm described in: [Phase-only image based kernel estimation for single image blind deblurring, CVPR2019.](https://openaccess.thecvf.com/content_CVPR_2019/papers/Pan_Phase-Only_Image_Based_Kernel_Estimation_for_Single_Image_Blind_Deblurring_CVPR_2019_paper.pdf) 

<pre>
@inproceedings{pan2019phase,  
    title={Phase-only image based kernel estimation for single image blind deblurring},         
    author={Pan, Liyuan and Hartley, Richard and Liu, Miaomiao and Dai, Yuchao},    
    booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},      
    pages={6034--6043},     
    year={2019}     
} 
</pre>

How to use
----------------
The codes are tested in MATLAB 2015b (64bit) under ubuntu 14.04 LTS 64bit version with an Intel Core i7-4790 CPU and 6 GB RAM.

1. Unpack the package.      
2. Run "main_uniform.m".     

### User-specified parameters:

### There are a few parameters that need to be specified by users.

* 'needsys'    :   1 for synthetic testing data. Blurred the image with a given kernel.     
* 'motion'     :   1 for linear kernal.       
* 'fast'       :   1 for fast processing strategy without coarse-to-fine.       
* 'kernel_size':   The size of blur kernel.       
* 'auto_size'  :   The scale for autocorrelation.           
* 'iter_num'   :   Iteration number.          
* 'lambda_grad' & lambda_l0 & lambda_tv: the weight for the L0/TV regularization.         

Notes 
----------------
1. The uploaded version is for tackling uniform motion blur. If you want to handle non-uniform blur, please refer to our paper. You need to segment an input blur image to overlapping square patches. For example, we segment the image into small patches in 80*80 with overlapping at 30 pixels.

2. If you obtain incorrect result (the algorithm converges to a local minimum), please re-try with slightly changed parameters (e.g., using large blur kernel sizes, fast, or iter_num).  

3. The results of our method on datasets ('Levin', 'Kohler', 'Gong' and 'Pan') are saved under the folder 'result'. 

4. Should you have any questions regarding this code and the corresponding results, please contact Liyuan.Pan@anu.edu.au


### Reproducing
1. Exactly following the equations in our paper.
2. Coding (e.g., FFT or correlation) by yourself instead of using in-built functions. 

### Please also cite the following papers if you use the code to generate data (e.g., images, tables of processing times, etc.) in an academic publication. 
  [1] Li Xu, Cewu Lu, Yi Xu, and Jiaya Jia. Image smoothing via l0 gradient minimization. ACM Trans. Graph., 30(6):174, 2011     
  
  [2] S. Cho, J. Wang, and S. Lee, Handling outliers in non-blind image deconvolution. ICCV 2011.           
  
  [3] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang, Blind Image Deblurring Using Dark Channel Prior. CVPR, 2016.      
  
  [4] Liyuan Pan, Yuchao Dai, Miaomiao Liu, and Fatih Porikli. Simultaneous stereo video deblurring and scene flow estimation. CVPR. 2017.     
       
