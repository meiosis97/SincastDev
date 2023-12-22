# Knn regression
A log file

## 2023-12-21与2023-12-22
1. 尝试PC了回归。
2. 发现似乎将PC回归的回归系数还原到高维回归系数时，似乎也可得到较好的变量选择。
3. PC会更迅速，因此在每个细胞做random forest也成了可能。
4. 但是PCA作为因子分析的方式，再高维度上的可解释性非常差，所以考虑做sparse PCA，估摸着sPCA能更进一步的提高变量选择的可靠性。
5. 想要可控制每一个loading变量选择的多少。大概会用soft-shresholding的方法，这涉及到了类lasso方法的最优化算法。
6. 学习了Alternating Direction Method of Multipliers (ADMM) 算法，里面有关于L1问题的解法。
7. TODO：细读以下文献：
8. [CeSpGRN: Inferring cell-specific gene regulatory networks from single cell multi-omics and spatial data](https://www.biorxiv.org/content/10.1101/2022.03.03.482887v3.full) (cell-specifc network)。
9. [Cell-specific network constructed by single-cell RNA sequencing data](https://academic.oup.com/nar/article/47/11/e62/5377474?login=false) (cell-specifc network)。
10. [Constructing local cell-specific networks from single-cell data](https://www.pnas.org/doi/10.1073/pnas.2113178118) (cell-specifc network)。
11. [c-CSN: Single-cell RNA Sequencing Data Analysis by Conditional Cell-specific Network](https://www.sciencedirect.com/science/article/pii/S1672022921000589) (cell-specifc network)。
12. [Learning cell-specific networks from dynamics and geometry of single cells](https://www.biorxiv.org/content/10.1101/2023.01.08.523176v3) (cell-specifc network)。
13. [spliceJAC: transition genes and state-specific gene regulation from single-cell transcriptome data](https://pubmed.ncbi.nlm.nih.gov/36321549/) (cell-specifc network)。
14. [Regression Shrinkage and Selection via the Lasso](https://www.jstor.org/stable/2346178) (原始lasso的解法)。
15. [Regularization Paths for Generalized Linear Models via Coordinate Descent](https://pubmed.ncbi.nlm.nih.gov/20808728/) (glmnet里面lasso的解法)。
16. [Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers](https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf) (原版ADMM算法).
17. [With Applications to Machine Learning：chapter 16](https://www.cis.upenn.edu/~cis5150/ws-book-IIb.pdf) (教科书对ADMM的解释，有更多的例子).
18. [Sparse Principal Component Analysis](https://hastie.su.domains/Papers/spc_jcgs.pdf) (原始sparse pca)。
19. [The Bayesian Lasso](http://www.math.chalmers.se/Stat/Grundutb/GU/MSA220/S18/bayeslasso.pdf) (延申阅读)
20. [A Review of the Spike-and-Slab LASSO](https://arxiv.org/abs/2010.06451) (延申阅读)
21. 在convex analysis里面，似乎indicator function都是0或无穷。

## 2023-12-14
1. 尝试了以单细胞的k neighbour为单位进行lasso回归。
2. 尝试了以变量选择是否来对umap上色，尝试了对变量选择是否进行diffusion。
3. 尝试了对变量选择是否矩阵的PCA。细胞可被这个矩阵分开。
4. 和导师们开了会，稍微了解了下大概课题可能存在的问题。
5. 首先导师们对结果的可解释性提出了质疑，即如何将这些结果推销给生物学家，这个方法的生物学应用在哪。我个人对这点上没有什么问题。我解释了可以从tanscriptional factor analysis入手，分析单细胞的转录激活状态，得到了KA的认可。
6. Jarny未提出很多的问题，主要在想大家可能criticize这个课题的理由。其中一点是对数据的Sincast imputation这一步，我也觉得不太妥，因为整个流程用了很多knn这个idea，
7. Jarny提出这个方法的可扩展性，提出了GSEA分析的可能性，我个人认为可能，即在每一个细胞上的变量选择后，对选择的变量进行GSEA分析，在每一个细胞上可计算GSEA score并进行可视化。
8. KA提出了对knn中k的数量，以及变量选择数量的疑惑，这需要后续分析。
9. JD提出了对Sincast imputation的可替代方案。首先，Sincast imputation可以被替换为svd low rank approximation。其次，可以不做任何imputation，直接做PC regression，在回归分析中，PC regression也可以看作是一种对数据的imputation，特别是在线性模型，我们可以通过PC regression的回归系数推导gene的回归系数。
10. 但我还是对PC regression的可解释性存疑，特别要是我们用非线性regression模型的情况下。且无法想象PC regression对如SIENIC之类transcriptional factor analysis的应用。
11. 想要尝试对regression coefficient matrix以细胞为向量空间的PCA，是否得到的PC score为meta regression model? (考虑restricted loading, 如loading为正值且相加的和为1).
12. 想要尝试对lasso regression coefficient matrix的zero-preserving low rank approximation, 是否得到的lra能更好的代表变量选择的结果。

## 2023-12-6
1. 尝试了以单细胞的k neighbour为单位进行多元线性回归。
2. 尝试了对数据整体做scaling，而不是在每个neighbourhood中做独立的scaling。不然不同neighbourhood所得到的结果不可比较。
3. 发现在自变量多的情况下，12月5日提到的unstable fit现象不见了。
4. 我们尝试了对EGR2做多元线性回归，自变量包括CCL3,CCL4,CCL5,EGR3,EGR4,CLEC9A,VCAN等与免疫有关的细胞标志物。
5. 发现当我们仅用CLEC9A (总体极低表达，但在某个单独的细胞簇高表达) 作为自变量 (或多考虑额外的1，2个自变量) 的时候，会出现unstable fit的现象。
6. 但当我们将模型中的自变量数量提升的时候unstable fit的现象消失。
7. 然而，我们发现regression coefficient的大小并无法反映自变量的绝对重要性。
8. 因此，我们尝试了用glmnet进行ridge regression。
9. 在没有tuning的情况下，我们选择了ridge regression path上的第100个regression fit (最小的penalty) 为最终的结果。发现即使在penalty很小的情况下，ridge regression极大的提升了regression coefficient的可解释性。
10. 我们提出了如下的predictivity和stability的测量方案。我们计算了一个基因在所有的细胞（neighbourhood）中regression coefficient拟合的均值为该基因对因变量的predictivity。其coefficient拟合的coefficient of variation (cv) 被用于其stability的测量。cv的计算为标准差处于均值。我们之所以用cv而不是用variance是因为我们发现方差与均值存在强烈的正相关，因此我们想去除均值对dispersion estimation的影响。不然predictivity和stability也会成正相关。
11. 因为在每一个细胞上我们都得到了一个表示regression coefficient的向量，我们尝试了对regression coefficient matrix进行降维，发现细胞簇在这个低维空间上有较好的分离。这代表了cells are clustered by transcriptional regulation on EGR2。
12. 通过将某一个基因的regression coefficient投影在umap上，我们可以判断是哪一个transcriptional regulation的变化 (如启动，或关闭) 导致了细胞的transition。
13. 下一步想要试一下在单个细胞上用binary value来表示regulation的on和off，如使用lasso regression进行变量选择，选中的变量代表on的regulation。尝试用Sincast来impute binary expression。
14. 下一步想要尝试其他非线性模型。特别是genie3 gene regulatory network的应用。

## 2023-12-5
1. 尝试了以单细胞的k neighbour为单位进行单元线性回归。
2. 发现regression coefficient在一些情况下变的非常大 (疑似为unstable fit所造成的outlier，这可能与k太小有关）。
3. 发现在增加k后，还是会出现相同的情况。
4. 考虑是因为基因的scale不同所以导致有些regression coefficient特别大。如scale在0.0001的基因回归分析scale在1的基因。
5. 尝试在每一个neighbourhood中进行回归分析的时候，将自变量和因变量均一化。目的是让regression coefficient的大小反映变量重要性。
6. 问题还是存在。且得考虑上述第5点的合理性：即能否在每一回归分析的时候，对数据做独立的均一化。

## 2023-12-4
1. 尝试了以单细胞的k neighbour为单位计算基因之间的correlation (i.e, EGR2 vs CCL3)。
2. 计算完correlation之后，对每一个单细胞上计算的correlation，使用Sincastd diffusion operator进行降噪处理。
3. 将correlation映射在umap上发现有趣，且好看的pattern，故判断knn regression可行。
4. correlation的问题在于它处于0到1之间，那么降噪之后对数据的scaling是一个需要考虑的问题。
5. 尝试了使用median matching和quantile matching的posting-imputation scaling的方法。
6. 这两种方法存在较大的差异，且都不能保证post-imputation scaling之后的correlation能保持在0到1之间。
7. 当diffusion time (steps of markov random walk) 高的时候，correlation 数据变成了噪音。值得研究。
8. 尝试了将neighbourhood size设为k，但在计算correlation的时候只考虑k/2个细胞。这样做的目的是为了增加计算的稳定性 (减少k对计算的影响)。
