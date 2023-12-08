# Knn regression
A log file

## 2023-12-4
1. 尝试了以单细胞的k neighbour为单位计算基因之间的correlation (i.e, EGR2 vs CCL3)。
2. 计算完correlation之后，对每一个单细胞上计算的correlation，使用Sincastd diffusion operator进行降噪处理。
3. 将correlation映射在umap上发现有趣，且好看的pattern，故判断knn regression可行。
4. correlation的问题在于它处于0到1之间，那么降噪之后对数据的scaling是一个需要考虑的问题。
5. 尝试了使用median matching和quantile matching的posting-imputation scaling的方法。
6. 这两种方法存在较大的差异，且都不能保证post-imputation scaling之后的correlation能保持在0到1之间。
7. 当diffusion time (steps of markov random walk) 高的时候，correlation 数据变成了噪音。值得研究。
8. 尝试了将neighbourhood size设为k，但在计算correlation的时候只考虑k/2个细胞。这样做的目的是为了增加计算的稳定性 (减少k对计算的影响)。

## 2023-12-5
1. 尝试了以单细胞的k neighbour为单位进行单元线性回归。
2. 发现regression coefficient在一些情况下变的非常大 (疑似为unstable fit所造成的outlier，这可能与k太小有关）。
3. 发现在增加k后，还是会出现相同的情况。
4. 考虑是因为基因的scale不同所以导致有些regression coefficient特别大。如scale在0.0001的基因回归分析scale在1的基因。
5. 尝试在每一个neighbourhood中进行回归分析的时候，将自变量和因变量均一化。目的是让regression coefficient的大小反映变量重要性。
6. 问题还是存在。且得考虑上述第5点的合理性：即能否在每一回归分析的时候，对数据做独立的均一化。

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
