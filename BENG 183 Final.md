*Final Paper Grading Rubric*

<https://zhonglab.gitbook.io/3dgenome> 3D genome book chapters

# Principal Component Analysis for Single-Cell RNA-seq

**Shuo Qin, Minzhou Ou, Margaret Zhang**

_BENG 183 – Final Project Handout_

---

## Abstract

Modern genomics generates extremely high-dimensional data: thousands of genes with thousands of cells, along with a lot of noise. Principal Component Analysis (PCA) is one of the most widely used dimensionality reduction methods to make such data analyzable. PCA transforms correlated gene expression measurements into a small number of linearly uncorrelated variables, which are now called principal components (PCs), that capture the dominant patterns of variation. In this chapter, we will introduce the intuition and mathematical foundations of PCA, explain why PCA is a natural first step for single-cell RNA-seq (scRNA-seq) analysis, and summarize a benchmarking study that compares PCA against a few other dimensionality reduction methods on 32 scRNA-seq datasets [1]. We then discuss key limitations of PCA and connect it to several bioinformatics applications that build on concepts from this course (normalization, clustering). Finally, we briefly review current software tools that implement PCA in standard scRNA-seq workflows.

---

## 1. Background

### 1.1 Why dimensionality reduction?

In bulk RNA-seq like scRNA-seq, each sample is represented by expression levels for tens of thousands of genes. In scRNA-seq, this becomes even more extreme, where each cell is a point in an approximately 20000-dimensional space with genes, with many zeros due to dropouts and low expression. Similar high-dimensional structures appear in other omics, which is why dimensionality reduction is a recurring theme in genomics.

Working directly with such high dimensions is problematic:

- Distances become hard to interpret (“curse of dimensionality”).
- Noise dominates: technical variation and sequencing depth differences can overshadow biological signals.
- Computation becomes expensive: clustering and graph-based methods scale poorly with the number of features.

Dimensionality reduction addresses these issues by compressing the data into a lower-dimensional representation that preserves the most relevant variation.

### 1.2 What does PCA do, conceptually?

Principal Component Analysis (PCA) finds new axes, the PCs, that are:

- Linear combinations of the original variables (genes).
- Orthogonal to each other (uncorrelated).
- Ordered by variance, so PC1 explains the maximum possible variance, PC2 the next most, and so on.

By projecting each sample or cell onto the first few PCs (for example, 10–50 instead of 20000 genes), we obtain:

- A compressed, denoised representation of the data.
- Low-dimensional coordinates for visualization (PC1 vs. PC2) or for downstream methods like k-means clustering or trajectory inference.

In many scRNA-seq pipelines, PCA sits between the basic data processing (quality control, normalization, highly variable gene selection) and the downstream analyses (clustering, UMAP or t-SNE, lineage reconstruction).

---

## 2. Mathematical Part of PCA

### 2.1 What does PCA do, mathematically?

Here, we briefly review the standard PCA procedure, covariance matrices, eigenvalues, and eigenvectors.

#### Standardization and Centering Data

We first standardize each variable:

$$
z = \frac{\text{value} - \text{mean}}{\text{standard deviation}}
$$

The equation transforms each variable to the same scale to avoid biased results.

#### Covariance Matrix Computation

For three variables \(x, y, z\), the covariance matrix has the form:

$$
\begin{bmatrix}
\mathrm{Cov}(x,x) & \mathrm{Cov}(x,y) & \mathrm{Cov}(x,z) \\
\mathrm{Cov}(y,x) & \mathrm{Cov}(y,y) & \mathrm{Cov}(y,z) \\
\mathrm{Cov}(z,x) & \mathrm{Cov}(z,y) & \mathrm{Cov}(z,z)
\end{bmatrix}
$$

The matrix equation tests the correlation between two variables:

- If positive: the two variables are correlated.
- If negative: the two variables are inversely correlated.

#### Identify Principal Components

If we have two eigenvectors with their corresponding eigenvalues:

$$
v_1 =
\begin{bmatrix}
0.6778736 \\
0.7351785
\end{bmatrix},
\quad \lambda_1 = 1.284028
$$

$$
v_2 =
\begin{bmatrix}
-0.7351785 \\
0.6778736
\end{bmatrix},
\quad \lambda_2 = 0.04908323
$$

If we divide the eigenvalue of each component by the sum

$$
\lambda_1 + \lambda_2
$$

we get the percentage of variance of the data. Therefore, we get each principal component's significance.

#### Create Feature Vector

Since vector 1 is more significant with \(1.284028 > 0.04908323\), the feature vector is formed with vector 1. This result is achieved by reducing dimensionality by 1:

$$
\text{FeatureVector} =
\begin{bmatrix}
0.6778736 \\
0.7351785
\end{bmatrix}
$$

This will be the vector as the principal component.

#### Final Data Projection

Eventually, the feature vector will be used to reorient data from original axes to be represented by principal components:

$$
\text{FinalDataSet} =
\text{FeatureVector}^\top
\cdot
\text{StandardizedOriginalDataSet}^\top
$$

---

## 3. Case Study: Benchmarking PCA in scRNA-seq

To understand how PCA compares to other dimension reduction (DR) methods in realistic single-cell settings, we focus on the benchmarking study by Sun et al. [1].

### 3.1 Evaluation workflow

Sun et al. designed a comprehensive test bench. Below is their project overview.

**Datasets:**

- 30 publicly available scRNA-seq datasets from GEO and the 10X Genomics website, plus 2 additional simulated datasets, for a total of 32 datasets.

**Methods:**

- 18 DR methods, including:
  - Generic linear methods: PCA, factor analysis (FA), independent component analysis (ICA), classical multidimensional scaling (MDS), non-negative matrix factorization (NMF).
  - Non-linear methods: t-SNE, UMAP, Isomap, locally linear embedding (LLE), diffusion maps.
  - And others.

**Downstream tasks evaluated on the low-dimensional space:**

- Neighborhood preservation: Jaccard index to quantify how local neighborhoods are preserved.
- Cell clustering: run clustering algorithms and compare cluster labels with known cell identities using normalized mutual information (NMI) and adjusted Rand index (ARI).
- Lineage or trajectory inference: input the DR results into trajectory methods and measure agreement with known trajectories using the Kendall correlation coefficient.

For each dataset, they also varied the number of components (for example: 2, 6, 14, 20) when applicable, to see how performance depends on dimensionality.

### 3.2 PCA’s role in the workflow

In this benchmark, PCA is treated as a serious baseline rather than a toy.

- PCA is run on normalized, log-transformed counts, often after selecting highly variable genes.
- The resulting PCs are used as input for clustering algorithms (k-means, hierarchical clustering, graph-based methods) and for trajectory inference tools.

This is similar to the pipelines implemented in tools like Seurat v3, where in the paper the author established a pipeline that covers the standard QC, normalization, HVG, PCA, neighbors, clustering, and UMAP in order [5].

### 3.3 Key findings about PCA

Sun et al. reach several important conclusions about PCA in the scRNA-seq context:

- **Competitive clustering performance**  
  When a reasonable number of components is retained, generic linear methods like PCA, FA, ICA, and MDS can achieve a very nice clustering performance compared to specialized DR methods for scRNA-seq, as measured by NMI and ARI.

- **Useful for trajectory inference**  
  For lineage reconstruction, PCA belongs to a group of top performers with high Kendall correlations between inferred and true trajectories.

- **Stability**  
  PCA and other linear methods are highly stable across data splits and much faster than many non-linear or model-based methods, making them practical for large datasets.

Overall, the authors recommend PCA as a robust default choice, especially when analysts want to retain a moderate to large number of components. Non-linear methods like UMAP may be preferable when only a very small number of dimensions is needed or when one cares primarily about visualization. This case study supports the idea that PCA, although very simple conceptually, is still a strong algorithm for scRNA-seq dimensionality reduction.

---

## 4. Assumptions & Limitations

### 4.1 Linearity

Since PCA produces principal components that are linear combinations of the original variables, it is most suitable for modeling linear relationships as it assumes linearity. For example, if we were to model the position of a person on a ferris wheel, the data would be circular, but PCA would still produce two principal components that are linear and orthogonal to each other. As shown in Fig 1, this does not capture the nature of the data well.

<img width="569" height="300" alt="image" src="https://github.com/user-attachments/assets/2b23f519-745f-4b50-90ad-348b063177d2" />
> Fig 1. (A) Modeling a person on a ferris wheel. Red arrows represent PCs from PCA. Due to linearity and orthogonality, the PCs do not model the data well. (B) Data is linear but not orthogonal. Red arrows represent PCs, but do not fit the axes of the data well.

### 4.2 Orthogonality

As mentioned, PCA produces principal components that are orthogonal or uncorrelated. The reason for this is because PCA aims at identifying maximum variance. By assuming orthogonality, PCA identifies principal components uncorrelated with each other, thus modeling the most variance. While this can be beneficial by design, the assumption of orthogonality can also limit the datasets we can model with a strong fit. As shown in Fig 1, a dataset that inherently has axes that are not orthogonal may not be suitably modeled via PCA.

### 4.3 Outliers

PCA is sensitive to outliers, which is a limitation for processing any dataset including extreme outliers. This is due to the fact that PCA looks for maximum variance – since outliers would have variance from other data points, having several extreme outliers would automatically skew the calculation of the principal component.

### 4.4 Sensitivity to Scale

For the same reason of finding variance, PCA is highly sensitive to scale as well. If a dataset contained data using different scales or units, results would be skewed as well due to the bias caused by scaling. In other words, due to the use of different scales, there would be a biased pattern of larger variances that would dominate the results of PCA.

For this reason, typically data should be standardized before PCA. This could be a limitation depending on the nature of the dataset.

### 4.5 Interpretation

Finally, a limitation that PCA poses is the difficulty of interpretation. Since PCA reduces dimensionality, principal components may lose real-world meaning in the process. One example can be using PCA to evaluate data from cancer samples. Suppose our aim is to identify cancer subtypes – what we can see from PCA is groups of genes that have similar patterns across their change in expression levels. However, even if genes show similar patterns in expression levels, they are not necessarily involved in the same biological processes, which leaves us with the question of what the principal components actually stand for in the biological sense.

## 5. Further Applications of PCA

### 5.1 RNA-seq and Microarray Analysis

PCA is often done on RNA-seq and microarray samples. The first reason for using PCA is quality control – we expect that replicate samples should cluster together since they should not have a lot of variance. This is a good way to confirm replicability. Another use of PCA is experimentally identifying patterns across samples from different conditions. PCA gives us a general image of how different data points differ from each other, so it is a good point to start.

### 5.2 Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA is a method that builds a network based on correlation between the expression level of various genes and assigns highly correlated genes to modules. In this process, module eigengenes (MEs) are calculated – this is the PC1 of a module of genes. MEs are used to represent the overall profile of a module of genes. Such modules allow us to further explore biological processes correlated with certain conditions.

### 5.3 Genome Wide Association Study (GWAS)

Ancestry is an important part of GWAS, as samples from the same ancestry tend to exhibit similar patterns. In certain studies, this similarity can cause false or skewed results. Hence, we can separate individuals by ancestries through PCA clusters before continuing further analysis.

## References

Sun, Shiquan, et al. “Accuracy, Robustness and Scalability of Dimensionality Reduction Methods for Single-Cell RNA-Seq Analysis.” *Genome Biology*, vol. 20, Dec. 2019, p. 269. PubMed Central, https://doi.org/10.1186/s13059-019-1898-6.

“Principal Component Analysis (PCA): Explained Step-by-Step.” *Built In*, builtin.com/data-science/step-step-explanation-principal-component-analysis. Accessed 5 Dec. 2025.

“What Is Principal Component Analysis (PCA)?” *IBM*, 17 Nov. 2025, www.ibm.com/think/topics/principal-component-analysis. Accessed 05 Dec. 2025.

Kamperis, Stathis. *Principal Component Analysis Limitations and How to Overcome Them*, 23 Feb. 2021, ekamperi.github.io/mathematics/2021/02/23/pca-limitations.html. Accessed 05 Dec. 2025.

Stuart, Tim, et al. “Comprehensive Integration of Single-Cell Data.” *Cell*, vol. 177, no. 7, June 2019, pp. 1888–1902.e21. PubMed Central, https://doi.org/10.1016/j.cell.2019.05.031.

