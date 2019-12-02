## quickpath

R package "quickpath" for quick pathway analysis

Author: Jiangyuan Li

Here is a [tutorial](https://rawcdn.githack.com/jiangyuan2li/quickpath/b5188c5ea4a08cf5ef04764cfd8b84036c70f8ae/vignettes/tutorial.html), which demonstrates the usage and convenience of this package.

### Installation
`devtools::install_github("jiangyuan2li/quickpath")`

### Description
Pathway analysis is a common task in bioinformatics research. The data for pathway analysis come from high throughput biology, including high throughput sequencing data and microarray data. A comprehensive pathway analysis consists of extracting differential expression information, doing statistical inference, and making summary figures. By storing frequently-used database locally and parallel computing, this package can speed up pathway analysis a lot. This package contains three common datasets locally (mouse, chicken, human pathway references). This can save querying time and provide stability when doing analysis.

There are two major methods, over-representation analysis and enrichment analysis. Here, we are focusing on over-representation analysis, specifically, hypergeometric distribution test (equivalently one-tailed version of Fisher's exact test) for RNA sequencing data. Over-representation pathway analysis works on a list of "significant" genes and "background" genes. The "significant" genes is usually given from "Cuffdiff" results. This package can extract information easily from "Cuffdiff" output.

After fast pathway analysis, it is easy to extract useful information from results and make summary plots. All the resulting tables and figures can be easily written and saved.

**Dependencies:** biomaRt, parallel, ggplot2, xlsx.
