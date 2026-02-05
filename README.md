## Problem
Metastasis - primary cause of cancer-related mortality in patients with solid tumors, yet the timing and risk of metastatic 
progression vary widely among individuals (Filipp et al., 2017). Advances in cancer genomics - enabled the collection of high-dimensional molecular data (somatic mutations and copy number 
alterations), offering new opportunities to better understand tumor heterogeneity and disease progression (Schramm et al., 
2025; Zheng & Frost, 2020). However, translating this complex genomic information into accurate predictions of time to metastasis remains a significant 
methodological challenge (Mobadersany et al., 2018).
## Dataset:  MSK-CHORD dataset (Link: https://www.cbioportal.org/study/summary?id=msk_chord_2024)
## Outcome Variable
Time to metastasis -  interval from study entry (baseline) to first documented metastatic even
## Proposed Analysis Design
### Objective 
To identify genomic phenotypes associated with metastatic progression and to evaluate their prognostic relevance.
### Summary 
A two-stage survival-informed modeling framework where Random Survival Forests are used to derive patient similarity and discover 
genomic phenotypes, followed by penalized Cox regression to quantify phenotype-specific metastasis risk.
## Expected Impact
● Identify high-risk patient subtypes 
● Quantify genomic heterogeneity effects on metastasis 
● Integrate ML-based clustering with interpretable survival modeling for precision oncology
