# SincastDev
A developing version of Sincast

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

## 2023-12-5
1. Initialize coding for SincastAtlas R object.
2. Added function PostScale for post imputation scaling.
3. PostScale did not function properly (Something looks wrong with gene-wise mean vs variance estimation)
