# SincastDev
A developing version of Sincast

## 2023-12-12

## 2023-12-11
1. We tried to project the Sincast imputed data (`x.imp`) on to the pc space of the original data, and then transform (or project?) it back to the original feature space.
2. Basically we had `x.imp %*% L %*% t(L)`, where `L` is the loading matrix of the PCA performed on the original data. We consider `L %*% t(L) = L %*% inv(t(L)%*%L) %*% t(L)` a projection matrix that projects `x.imp` on to the colunm space of `L`. We hope this step of projection can help adjust for over-smoothing.
4. We tried to calculate `L` on the scaled and unscale original data.
5. We decided to not scale the data before svd since (1) low rank approximation does not have zero mean and standardized variance assumption. (2) After svd, we tried to add back centers and scales of the original data, but the final lra outcome looks weired.
6. In order to benchmark the imputation and the above reconstruction result, we constructed a pseudobulk atlas of the data, and tried to project the imputed data, as well the reconstructed data onto the atlas.
7. We tested three different imputation results, (1) `x.imp` (2) `x.lra`, that is lra on the original data (3) `x.imp.lra`, that is the lra on `x.imp`.
8. The projection of none of the three imputation results onto the pseudobulk atlas looks correct. Especially cell clusters of `x.imp.lra` were projected onto the wrong direction of the atlas.
9. We suspect that it is the quantile or median matching step after either the Sincast or lra imputation that cause the problem. After imputation, we tried to scale (only scale, not mean shift) each imputed feature such that for any cells that has non-zero value on this feature before imputation, the quantile of these cells after imputation should be the same.
10. Without doing this step of quantile matching, all the projections looks mcuh better: cell clusters can be projected onto the correct direction of the atlas, but do not exactly match to the cell clusters of their atlas conterpart (under-projection).
11. Why the quantile matching causes the problem? Using the first feature of the original data as an example, it is almost not expressed by all of the cells, expect a fews that are randomly distributed on the umap. After both the Sincast and lra imputation (without quantile matching), we see that the imputed expression value are very small, close to zero. If we do our proposed quantile matching, we will up-scale the imputed expression to those noise in the original data, making every cell highly express this noisy gene. 
12. We also question the necessaity of quantile matching after lra imputation. Since high rank lra will eventually reconstruct the signal of the original data. If we want to match the imputed data with the original data, then what is the point of denoising? We consider the lra imputation is already in the correct scale.
13. From now on, we skip the quantile matching step.
14. We tried the zero-preserving lra method proposed by [Linderman et al. (2022)](https://www.nature.com/articles/s41467-021-27729-z). The method is developed to preserve biological zeros in the data after imputation.
15. Suprisingly, the zero-preserving lra imputed data makes perfect projection onto the pseudobulk atlas, and it's sparisty matches the sparisty of the pseudobulk data. This could be an evidence that zero-preserving lra correctly recovered biological zeros.
16. Next, we tried zero-preserving lra on the Sincast imputed data. The resulting reconstruction also makes perfrect, but less noiser projection compared to which made by the zero-preserving lra on the original data. This suggests that Sincast imputation also makes imputation on biological zeros.
17. Given what we observed, we propose the following pipeline, (1) Sincast imputation, (2) Zero-preserving lra on the Sincast imputed data to correct for biological zeros (this is similar to the stacking idea in machine learning). (3) Zero-preserving lra on the original data. We assume that the reconstruction correctly preserve the scale of the original data. (4) Scale the Sincast-lra imputed data such that its feature quantiles matche with the lra imputed data. Here, we not only consider quantiles of expressed values, but also zeros since they represent true biology.

## 2023-12-10
1. We tried to perform low rank approximation (lra) on the original data, and treat lra as another imputation result, our goal is to combine two imputation results.
2. The idea is the same as that is in the machine learning meta-anlaysis, where we can combine the result of the two predictions made by two different methods. Here, we can try to combine two different imputation results.
3. The combination method we tried was, aggregation (e.g, take the average of the lra and Sincast imputation result), stacking (e.g, perform Sincast imputation first, and then perform lra on the Sincast imputed data), boosting (it is relatively uneasy to interpret boosting method in the case of imputation, basically, in machine learning, we use the model to predict the outcome, and according to the prediction result, assign weights to samples and perform another iteration of prediction. At the end, the prediction results get from this itretaion are aggregated. We thought about at each diffusion step, perform a lra, and aggreagte all lra).
4. Boosting would take a lot of time to do experiments. Therefore we tested aggreagtion and stacking first.
5. It is hard to image that stacking can in any way correct the scale of Sincast imputed data.
6. Therefore, we tested aggregation.
7. The post scaling pipline strat by performing two seperate imputation by lra and Sincast. After each imputation, we tried to scale each imputed feature such that for any cells that has non-zero value on this feature before imputation, the quantile of these cells after imputation should be the same.
8. Compared to what we did on 2023-12-7, we made the following bayes like assumption: for each observed value $x_{ij}$, we have two imputed values $\theta_{ij}^{lra}$ and $\theta_{ij}^{Sincast}$. We assume that these two imputed values are generated from a normal prior distribution with mean $E[\theta_{ij}] = (\theta_{ij}^{lra} + \theta_{ij}^{Sincast})/2$. And the variance is emperically estimated as by the variance of $\theta_{ij}^{lra}$ and $\theta_{ij}^{Sincast}$. The likelihood variance is given by $(x_{ij} - E[\theta_{ij}])^2$.
9. Though these variance estimation are herustic and **do not follow the rule of probability**, we still can make intutive aggregation of the Sincast and lra imputed data based on these illy defined variances.
10. If there is great discrepency between  $\theta_{ij}^{lra}$ and $\theta_{ij}^{Sincast}$, then $Var[\theta_{ij}]$ should be large, and any of the imputation result should be untrustfull.
11. How to determine whether the discrepency is large or not? We benchmark $(\theta_{ij}^{lra} - \theta_{ij}^{Sincast})^2$ by $(x_{ij} - E[\theta_{ij}])^2$.
12. If $(\theta_{ij}^{lra} - \theta_{ij}^{Sincast})^2$ is much larger than $(x_{ij} - \hat{\theta_{ij}})^2$, we let the final imputed value $\hat{x}\_{ij}$ be closer to $\hat{\theta}\_{ij}$, otherwise it should be closer to the original data, $x_{ij}$.
13. We found that shrink back the imputed data back to the original data may not be a good practice, since the original data can be very noisy and is not the ground truth representation of the biological state of a cell.
14. An interesting feature on which we found both the imputation method give weird result was the first gene (index 1, see details in the R script). This gene lowly, sparsely and irregularly expressed on the umap of the data, and from an objective intepretation, this gene could be completely a noisy gene. However, both lra and Sincast imputation reconstruct a specific expression pattern on this gene. Moerver, after imputation by both the methods, this gene is highly expression in may cells. This makes us wondering whether our imputation method, as well as lra are introducing too much of signal to the data, which should be noise. 

## 2023-12-8
1. We tried to do low rank approximation by PLS, where Y is set to the original data and X is set to the imputed data.
2. It is almost imposible to perform PLS in this case due to computation reason because PLS require svd on `svd(t(Y)%*%X)`, which is quite heavy.
3. Therefore, we first try to perform pca seperately on each of the original and the imputed data, with the same number of pcs.
4. Then we perform PLS on the pc space of the two data matrices.
5. We first tried to reconstruct the pc scores of the imputed matrix using its PLS loadings, and then reconstruct the full imputed matrix by its PC scores. The code was ` pcs.imp %*% pls.imp.L %*% t(pls.imp.L) %*% t(pcs.imp.L)`.
6. The intepertation of this reconstruction is mnot intuitive.
7. Seems like the imputed matrix after reconstruction still is not in the correct scale (does not match with the scale of the original data). Especially, signal of lowly expressed gene is completely buried.
8. As Jiadong suggested, tried to center and scale matrices before low rank approximation. We tried to center and scale a matrix, perform low rank approximation on the centered and scaled matrix, and then add centers and scales back to the low rank reconstruction.
9. We worried that two approximation steps would make the data too noisy.
10. We did the above procedure on the original data, and found that we the reconstructed matrix is not in the correct scale, and got a lots of negative values.
11. We are also not sure about in what way should we reconstruct the matrix, that is, should we reconstruct X, or Y, or $\hat{Y}$ that is predicted by PLS?

## 2023-12-7
1. Try to update Sincast post imputation scaling: PostScale
2. Previously, we calculated differences between actual values and imputed values, based on which we computed squared errors. We also calculated on each imputed value its post-imputation variance. We than perform post-imputation scaling by weighted averaging the actual value and imputed value, where weights are based on squared errors and post-imputation variance.
3. This is not intuitive, instead, we should weight by squared errors and prior-imputation variance. Squared errors should give more weights to the original data, and prior-imputation variance should favour the imputed data.
4. Tried a normal-normal bayesian framework for weights calculation.
5. Let imputed data be the prior normal mean, let original data be the observed value, we want to estimate posterior normal mean as a way to shrink the imputed data back to the observed data.
6. The problem is we don't know what is the prior normal variance, as well the likelihood variance.
7. Propose to emperically estimate these variance.
8. Propose to first emperically estimate the mean of each observation by low rank approximation. Prior normal variance is given by the square error between the low rank approximation and the imputed value, while the likehood variance is given by the square error between the low rank approximation and the original value.
9. Calculate posterior mean on each observation
10. Problem: defining a good rank used for low rank approximation is hard. If the rank is too high, the low rank approximation will be close to the original data and hence the imputed data could be negnect.
11. Jiadong suggested to center and scale the original data before applying low-rank approximation on it.
12. Todo: test low rank approximation by CCA or PLS
   
## 2023-12-6
1. Updated .gitignore.
2. git consists of a local repository and a local folder on which we can directly makes changes and modification.
3. When run `git add`, we can add a modification to a **staging** block, in which modification are waiting to be commited (add to the repository).
5. When run `git commit`, we commit the modification at the staging block to update the repository.
6. local -> stage -> local repository -> remote repository
7.  .gitingore cannot ignoure files that are already committed in the repository, to ignore these files, they need to be deleted from the repository by `git git rm --cached filename`
8. `git push` will return an error if the remote repository contains a file that are not found in the local repository, and the file is not logged to be deleted. 
9. `git push` will return an error if the remote repository is not the same as the local repository before commitment (i.e, untracked changes)???

## 2023-12-5
1. Initialize coding for SincastAtlas R object.
2. Added function PostScale for post imputation scaling.
3. PostScale did not function properly (Something looks wrong with gene-wise mean vs variance estimation)
