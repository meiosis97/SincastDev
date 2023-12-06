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
2. 发现regression coefficient在一些情况下变的非常大。怀疑是k太小造成的。
3. 发现在增加k后，还是会出现相同的情况。
4. 考虑是因为基因的scale不同所以导致有些regression coefficient特别大。如scale在0.0001的基因回归分析scale在1的基因。
5. 尝试在每一个neighbourhood中进行回归分析的时候，将自变量和因变量均一化。
6. 问题还是存在。且得考虑上述第5点的合理性：即能否在每一回归分析的时候，对数据做独立的均一化。

