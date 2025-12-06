# Principal Component Analysis for Single-Cell RNA-seq

**Shuo Qin, Minzhou Ou, Margaret Zhang**  
BENG 183 – Final Project Handout  

---

## Abstract

Modern genomics generates extremely high-dimensional data: thousands of genes with thousands of cells, along with a lot of noise. Principal Component Analysis (PCA) is one of the most widely used dimensionality reduction methods to make such data analyzable. PCA transforms correlated gene expression measurements into a small number of linearly uncorrelated variables, called principal components (PCs), that capture the dominant patterns of variation.

In this chapter, we introduce the intuition and mathematical foundations of PCA, explain why PCA is a natural first step for single-cell RNA-seq (scRNA-seq) analysis, and summarize a benchmarking study that compares PCA against other dimensionality reduction methods on 32 scRNA-seq datasets [1]. We then discuss key limitations of PCA and connect it to several bioinformatics applications that build on concepts from this course (normalization, clustering). Finally, we briefly review current software tools that implement PCA in standard scRNA-seq workflows.

---

## 1. Background

### 1.1 Why dimensionality reduction?

In bulk RNA-seq and scRNA-seq, each sample is represented by expression levels for tens of thousands of genes. In scRNA-seq, this becomes even more extreme, where each cell is a point in an approximately 20,000-dimensional space, with many zeros due to dropouts and low expression. Similar high-dimensional structures appear in other omics, which is why dimensionality reduction is a recurring theme in genomics.

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

By projecting each sample or cell onto the first few PCs (for example, 10–50 instead of ~20,000 genes), we obtain:

- A compressed, denoised representation of the data.  
- Low-dimensional coordinates for visualization (PC1 vs. PC2) or for downstream methods like k-means clustering or trajectory inference.

In many scRNA-seq pipelines, PCA sits between the basic data processing (quality control, normalization, highly variable gene selection) and the downstream analyses (clustering, UMAP or t-SNE, lineage reconstruction).

---

## 2. Mathematical Part of PCA

Here we briefly review the standard PCA procedure: standardization, covariance matrices, eigenvalues, eigenvectors, and projection.

### 2.1 Standardization and centering data

To avoid biased results due to different scales, we standardize each variable using the z-score:

$$
Z = \frac{\text{value} - \text{mean}}{\text{standard deviation}}
$$

This transforms each variable to the same scale.

### 2.2 Covariance matrix computation

After standardization and centering, we compute the covariance matrix. For three variables \(x, y, z\), the covariance matrix has the form:

$$
\begin{bmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\
\sigma_{yx} & \sigma_{yy} & \sigma_{yz} \\
\sigma_{zx} & \sigma_{zy} & \sigma_{zz}
\end{bmatrix}
$$

The covariance measures how two variables change together:

- If covariance is positive: the two variables are positively correlated.  
- If covariance is negative: the two variables are inversely correlated.

### 2.3 Identify principal components (eigenvalues and eigenvectors)

PCA finds the eigenvectors and eigenvalues of the covariance matrix. For a simple 2D example, suppose we obtain two eigenvectors with their eigenvalues:

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

If we divide each eigenvalue by the sum $(\lambda_1 + \lambda_2)$, we get the percentage of variance explained by each component. This tells us how important each principal component is.

### 2.4 Create feature vector

Since component 1 is more significant $\(\lambda_1 = 1.284028 > 0.04908323\)$, we can reduce dimensionality by keeping only the first eigenvector:

$$
\text{FeatureVector} =
\begin{bmatrix}
0.6778736 \\
0.7351785
\end{bmatrix}
$$

This is our principal component direction.

### 2.5 Final data projection

Finally, we use the feature vector to reorient the data from the original axes into the principal component space:

$$
\text{FinalDataSet} =
\text{FeatureVector}^\top
\times
\text{StandardizedOriginalDataSet}^\top
$$

The projected data in this new coordinate system is lower-dimensional and captures most of the variance of the original dataset.

---

## 3. Case Study: Benchmarking PCA in scRNA-seq

To understand how PCA compares to other dimensionality-reduction (DR) methods in realistic single-cell settings, we focus on the benchmarking study by Sun et al. [1].

### 3.1 Evaluation workflow

Sun et al. designed a comprehensive test bench. Their project overview:

- **Datasets**  
  - 30 publicly available scRNA-seq datasets from GEO and the 10x Genomics website  
  - 2 additional simulated datasets  
  → 32 datasets in total

- **Methods**  
  18 DR methods, including:
  - Generic linear methods: PCA, factor analysis (FA), independent component analysis (ICA), classical multidimensional scaling (MDS), non-negative matrix factorization (NMF)  
  - Non-linear methods: t-SNE, UMAP, Isomap, locally linear embedding (LLE), diffusion maps  
  - And other specialized methods for scRNA-seq

- **Downstream tasks evaluated on the low-dimensional space**  
  - **Neighborhood preservation:** Jaccard index to quantify how local neighborhoods are preserved  
  - **Cell clustering:** clustering algorithms followed by comparison to known cell identities using normalized mutual information (NMI) and adjusted Rand index (ARI)  
  - **Lineage / trajectory inference:** DR results used as input to trajectory methods; agreement with known trajectories measured using the Kendall correlation coefficient  

For each dataset, they also varied the number of components (for example, 2, 6, 14, 20) when applicable, to see how performance depends on dimensionality.

### 3.2 PCA’s role in the workflow

In this benchmark, PCA is treated as a serious baseline rather than a toy method:

- PCA is run on normalized, log-transformed counts, often after selecting highly variable genes.  
- The resulting PCs are used as input for clustering algorithms (k-means, hierarchical clustering, graph-based methods) and for trajectory inference tools.

This is similar to pipelines implemented in tools like Seurat v3, where the authors established a pipeline that covers standard QC, normalization, HVG selection, PCA, neighbor graph construction, clustering, and UMAP in order [5].

### 3.3 Key findings about PCA

Sun et al. reach several important conclusions about PCA in the scRNA-seq context:

- **Competitive clustering performance**  
  When a reasonable number of components is retained, generic linear methods like PCA, FA, ICA, and MDS can achieve very good clustering performance compared to specialized DR methods for scRNA-seq, as measured by NMI and ARI.

- **Useful for trajectory inference**  
  For lineage reconstruction, PCA belongs to a group of top performers with high Kendall correlations between inferred and true trajectories.

- **Stability and efficiency**  
  PCA and other linear methods are highly stable across data splits and much faster than many non-linear or model-based methods, making them practical for large datasets.

Overall, the authors recommend PCA as a robust default choice, especially when analysts want to retain a moderate to large number of components. Non-linear methods like UMAP may be preferable when only a very small number of dimensions is needed or when one cares primarily about visualization. This case study supports the idea that PCA, although simple conceptually, is still a strong algorithm for scRNA-seq dimensionality reduction.

---

## 4. Assumptions & Limitations

### 4.1 Linearity

Since PCA produces principal components that are linear combinations of the original variables, it is most suitable for modeling linear relationships and assumes linearity. For example, if we model the position of a person on a ferris wheel, the data would be circular, but PCA would still produce two principal components that are linear and orthogonal to each other. This does not capture the circular nature of the data well.

### 4.2 Orthogonality

PCA produces principal components that are orthogonal (uncorrelated). PCA aims to identify directions of maximum variance, and orthogonality ensures that each principal component captures unique variance. While this can be beneficial, the assumption of orthogonality can also limit datasets whose natural axes are not orthogonal; PCA may not align well with the intrinsic structure.

### 4.3 Outliers

PCA is sensitive to outliers. Because PCA looks for directions of maximum variance, extreme outliers can strongly influence the principal components and skew their orientation. Datasets with many extreme outliers may yield misleading PCA results unless outliers are removed or robust PCA methods are used.

### 4.4 Sensitivity to scale

For the same reason of maximizing variance, PCA is highly sensitive to scale. If a dataset contains variables measured on very different scales or units, results can be skewed because variables with larger variance dominate. Therefore, data should typically be standardized before PCA.

### 4.5 Interpretation

A final limitation of PCA is the difficulty of interpretation. Since PCA reduces dimensionality, principal components may lose clear real-world meaning. For example, when applying PCA to cancer gene expression data to identify subtypes, PCA reveals groups of genes that change together, but those genes are not necessarily involved in the same biological processes. This leaves us with the question of what each principal component represents biologically.

---

## 5. Further Applications of PCA

### 5.1 RNA-seq and microarray analysis

PCA is often applied to RNA-seq and microarray samples. One important use is **quality control**: replicate samples should cluster together in PCA space if the experiment is reproducible. PCA can also help identify patterns across samples from different conditions, providing a general overview of how samples differ and suggesting hypotheses for further analysis.

### 5.2 Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA builds a network based on correlations between the expression levels of various genes and assigns highly correlated genes to modules. In this process, module eigengenes (MEs) are calculated—this is the first principal component (PC1) of a module of genes. MEs represent the overall expression profile of a module. Such modules allow us to explore biological processes correlated with certain conditions.

### 5.3 Genome-Wide Association Study (GWAS)

Ancestry is an important factor in GWAS, as samples from the same ancestry tend to exhibit similar patterns. This similarity can cause false or skewed results if not accounted for. PCA on genotype dat
