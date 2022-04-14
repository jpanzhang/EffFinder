# EffFinder
EffFinder is an R package, for detecting the Normalization efficiencies of single reference gene (RGs) and multi-RG combinations
This package includes two main functions Eff.singleFinder and Eff.indexFinder.

(1) Eff.singleFinder. The program of Eff.singleFinder is very time-consuming, because the large p-value data matrix needs to calculate through a series of computationally intensive procedures. To save time and computational resources, some intermediate variables were kept in global environments for further use.

(2) Eff.indexFinder. The runing time of Eff.singleFinder is very short. You can use Eff.indexFinder to detect the Normalization efficiency of any normalization factor that you want to detect. The serial number of normalization factors can be found from a file named "Index file" generated by Eff.singleFinder function. So the Eff.singleFinder function should be run before Eff.indexFinder.

For more technical details about this algorithm, please see our publication: 
#### Jipan Zhang, Yongju Zhao: Normalization efficiency changes in reference genes (RG): evidence calling for use of multi-RG combination in qPCR experiments, 2022. (Unpublished)

The authors would be glad to hear how EffFinder is employed. You are kindly encouraged to notify Jipan Zhang <jpanzhang@live.com> if you have any trouble with this package.
## ________________________________________________________________________________
## EffFinder version 1.0.0 

