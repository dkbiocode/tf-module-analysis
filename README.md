# Network analysis of gene expression data

Weighted Gene Co-expression Network Analysis (WGCNA) is a principled approach for gene cluster identification and inferrence of gene expression modules.

## Background

The embryonic tissue of origin for the *C. elegans* intestine is, like all animals, the endoderm. Unique to *C. elegans*, however, is that the single intestine performs the functions normally taken on by the entire set of organs derived from endoderm in other animals. 

The evolutionarily conserved *Erythroid-Like Transcription* factors (ELTs) propogate specification and differentiation through embryonic development, culminating in ELT-2 and ELT-7 being expressed post-embryonically and, in the case of ELT-2, driving intestinal gene expression. This parallels the vertebrate orthologous gene family, GATAs, which play similar roles in development and interact in pairwise fashion in mammals. 

The nature and role of these two paralogous transcription factors is asymmetric, however, with the experimental perturbation of one (ELT-2, lethal) more devastating to intestinal development than the other (ELT-7, no drastic phenotype). However, the double deletion is far more deleterious than ELT-2d, suggesting a buffering, rather than redundant, role of ELT-7. 

Assymetric paralogues are described elsewhere in literature not only in development, but in cancer progression. The nature of the interaction between transcription factor pairs is critical to understanding the dynamics of tissue-specific traits as they emerge from development into maturity. 

These dynamics drive the gene regulatory network underlying tissue, organ, and tumor realization, and must be understood as a fundamental aspect of developmental processes. 


## Dineen et al.

The study Dineen et al. 2018, made perturbations of elt-2 and elt-7 using a rescue cassette to maintain a normal embryonic development. Following deletion of the cassette, the intestines were dissected and subjected to RNA-seq, providing 3-4 replicates across each of single and double mutants and a wild type control. 

The focus of this study was to compare gene clusters to ELT-2 binding data available at the time. The authors proposed various gene regulation models as a result. 

## Reanalysis with WGCNA 

The following project uses a more [recent method]() to define gene expression clusters, applying scale-free network assumptions to infer modular organization and find "hub genes" underlying these modules. 