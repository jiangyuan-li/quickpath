R package "quickpath" for quick pathway analysis

Author: Jiangyuan Li

Pathway analysis is a common task in bioinformatics research. The data for pathway analysis come from high throughput biology, including high throughput sequencing data and microarray data. A comprehensive pathway analysis consists of extracting differential expression information, doing statistical inference, and making summary figures with publication quality. Such process need to be repeated many times, if you are working for a lab that produces data constantly. This package "quickpath" is targeting to automatic the process above and get the results, including figures quickly.

There are two major methods, over-representation analysis and enrichment analysis. Here, we are focusing on over-representation analysis, specifically, hypergeometric distribution test (equivalently one-tailed version of Fisher's exact test) for RNA sequencing data.

Over-representation pathway analysis works on a list of "significant" genes and "background" genes. The "significant" genes is usually given from "Cuffdiff" results. The first function would be extracting the "significant" genes.

The pathway information are usually stored in a data collection such as KEGG, which takes time to query these databases. The most used pathways would be downloaded and distributed with the package to save querying time. The statistical test would be further performed in parallel.

Based on the results, the summary plots with publication quality are usually needed. We would rely on "ggplot2" to make figures about p values and percentages of differential expressed genes in some specific pathways.

**Dependencies:** biomaRt, KEGGREST, parallel, ggplot2.
