# LMM_SBD-Hyperspectral-Unmixing
![abstract](https://user-images.githubusercontent.com/56957343/122734524-6915e480-d293-11eb-87e7-026b27f5e388.JPG)
## 
This demo illustrates the unmixing of a 105x128x144 subset of the real Houston dataset with the LMM-SBD method introduced in the following paper:

S. G. Azar, S. Meshgini, S. Beheshti, and T. Y. Rezaii, "Linear Mixing Model with Scaled Bundle Dictionary for Hyperspectral Unmixing with Spectral Variability", Signal Processing, 2021. 
doi: https://doi.org/10.1016/j.sigpro.2021.108214

## Abstract:
Hyperspectral unmixing is the process of separating the signatures of different pure materials in a mixed pixel. For different reasons including the intrinsic variability of materials and variations in data collecting conditions, different forms of spectral signatures can be applied to a particular material. This is referred to as spectral variability. Scaling factors and bundle dictionary are two concepts that address illumination variations and intrinsic variabilities of materials, respectively. In this paper, we propose the linear mixing model with scaled bundle dictionary (LMM-SBD) method which combines both the scaling factors and bundle dictionary to benefit from the advantages of both approaches. Moreover, we use different spatial neighbors to account for the spatial coherence of the neighboring pixels and force their corresponding abundances to have a similar sparsity pattern by adding two mixed norms to the optimization problem. The produced problem is solved using an alternating direction method of multipliers approach. The proposed method is tested on a simulated data and three real hyperspectral data and the results verify that the proposed method is able to overcome the spectral variability of materials, exploit the spatial coherence of the neighboring pixels and yield a successful unmixing procedure compared with several state-of-the-art methods.
 ## 
Please kindly cite our paper if you find the codes useful. If you have any questions, please don't hesitate to contact the author at: sghanbariazar@tabrizu.ac.ir
