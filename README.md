---
editor_options: 
  markdown: 
    wrap: 72
---

# MSK-CHORD Dashboard: Interactive Exploration of Metastasis and Genomics

## Abstract

Problem : Metastasis is the leading cause of cancer mortality in solid
tumors, with widely variable timing and risk among patients (Filipp et
al., 2017). High-dimensional genomic data, including somatic mutations
and copy number alterations, offer insights into tumor heterogeneity and
progression (Schramm et al., 2025; Zheng & Frost, 2020). Yet, accurately
predicting time to metastasis from such complex data remains challenging
(Mobadersany et al., 2018). \## Methodlogy This project investigates
genomic determinants of metastatic progression using the MSK-CHORD
(Nature 2024) cohort of 25,040 tumors profiled with MSK-IMPACT. The
primary objective is to model **time to metastasis** by integrating
high-dimensional somatic mutation and copy number alteration data with
survival analysis.

A two-stage framework combines Random Survival Forest–based phenotype
discovery with penalized Cox regression to quantify metastasis risk. In
parallel, an interactive Shiny dashboard enables dynamic exploration of
genomic alterations, copy number frequencies, metastatic site
distributions, fusion networks, tumor mutation burden, and survival
outcomes across cancer types.

Together, this work bridges machine learning, survival modeling, and
interactive visualization to support precision oncology and
translational research.

## Dataset

**MSK-CHORD**\
Targeted sequencing of 25,040 tumors from 24,950 patients and their
matched normals via MSK-IMPACT, along with clinical annotations, some
derived from natural language processing (NLP).\
[Dataset link on
cBioPortal](https://www.cbioportal.org/study/summary?id=msk_chord_2024)

**Outcome Variable:**\
Time to metastasis — interval from study entry (baseline) to first
documented metastatic event.

### Expected Impact

-   Identify high-risk patient subtypes.
-   Quantify genomic heterogeneity effects on metastasis.
-   Integrate ML-based clustering with interpretable survival modeling
    for precision oncology.

## MSK-CHORD Dashboard

While the MSK-CHORD cohort is accessible through cBioPortal, this Shiny
dashboard provides enhanced exploratory analysis:

-   Dynamic cancer-type filtering with real-time updates across all
    plots.
-   Interactive CNA frequency visualization, showing precise hoverable
    amplification/deletion percentages.
-   Metastatic site summaries from clinical annotations.
-   Fusion chord diagram visualization, enabling intuitive exploration
    of recurrent gene fusion partners and network relationships (not
    available in the standard portal view).
-   Unified interface combining genomics, clinical variables, and
    survival analysis in one workflow.

Where cBioPortal provides high-level cohort summaries, this dashboard
emphasizes interactive sub-cohort analysis, network-level fusion
insight, scalability, and improved clinical interpretability, making it
well-suited for exploratory and translational research use cases.

#### Dashboard Access

The interactive Shiny dashboard can be deployed locally: \`\`\`R \# Run
locally shiny::runApp("MSK_Chord_dataset\MSK\_Dashboard")

or on a cloud service such as
[**https://indhirav.shinyapps.io/msk_dashboard/**](https://indhirav.shinyapps.io/msk_dashboard/){.uri}.\
