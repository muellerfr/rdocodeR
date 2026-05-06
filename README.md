# rdocodeR

<p align="center">
  <img src="man/figures/logo.png" height="260" alt="rdocodeR hex logo" />
</p>

`rdocodeR` is an R package for decoding cortical overlays into the NIMH RDoC framework and visualizing term-level and domain-level results.

## 🧠 What Is RDoC?

![RDoC Framework](man/figures/rdoc-framework-explainer.png)

Image source: [Boehringer Ingelheim - RDoC Framework Explainer](https://pro.boehringer-ingelheim.com/connecting-psychiatry/news-perspectives/rdoc-framework-explainer)

**Research Domain Criteria (RDoC)** is a research framework developed by the **U.S. National Institute of Mental Health (NIMH)** to study psychopathology in terms of **basic functional dimensions** (for example threat processing, reward learning, cognitive control) that vary along **normal-to-abnormal continua**, rather than starting from symptom-defined categories (DSM/ICD).

RDoC is operationalized through the **RDoC Matrix**, which organizes constructs into major domains and encourages integration across multiple **units of analysis** (genes, molecules, cells, circuits, physiology, behavior, self-report, and task performance) to test mechanistic hypotheses.

NIMH defines six major domains in the matrix. Domain definitions source: [NIMH Definitions of the RDoC domains and constructs](https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/definitions-of-the-rdoc-domains-and-constructs).

- ⚠️ **Negative Valence Systems**: systems supporting responses to aversive situations (for example acute threat/fear, potential threat/anxiety, sustained threat, loss, frustrative nonreward).
- 🌟 **Positive Valence Systems**: systems supporting responses to positive motivational situations (for example reward responsiveness, reward learning, habit).
- 🧠 **Cognitive Systems**: systems for core cognitive operations (for example attention, perception, declarative memory, language, cognitive control, and working memory).
- 👥 **Systems for Social Processes**: systems supporting social functioning (for example affiliation/attachment, social communication, and perception/understanding of self and others).
- 🌙 **Arousal and Regulatory Systems**: systems involved in global regulation (for example arousal, sleep-wake functioning, circadian rhythms, homeostatic regulation).
- 🦾 **Sensorimotor Systems**: systems supporting motor and action-related functions (for example action planning/selection, motor execution, agency/ownership).

Key sources: [Insel et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20595427/) · [NIMH: About RDoC](https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/about-rdoc) · [NIMH: RDoC Matrix](https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/constructs/rdoc-matrix) · [Boehringer Ingelheim: RDoC Framework Explainer](https://pro.boehringer-ingelheim.com/connecting-psychiatry/news-perspectives/rdoc-framework-explainer)

## 🚀 Installation

```r
install.packages("devtools")
devtools::install_github("alegiac95/rdocodeR")
```

> [!NOTE]
> `rdocodeR` does **not** ship precomputed term-specific null maps. Instead, it generates them locally the first time you set up or run the decoder. If no suitable Python environment is provided, `rdocodeR` will create a private virtual environment in the user cache, install the required Python packages there, and reuse that environment on later runs. You can also point the package to an existing Python interpreter with `python = "/path/to/python"` in `rdoc_setup()` or `rdoc_decode()`. The first setup run may take time and disk space because the default workflow generates and caches `1000` null maps per term.

### After Installation: Generate the Null Maps

Recommended explicit setup:

```r
library(rdocodeR)

rdoc_setup()  # generates and caches 1000 null maps per term
```

This will:

1. create or reuse a Python environment if needed
2. install the required Python packages for eigenstrapping
3. generate the cached null maps locally
4. reuse the same cached nulls on later runs

If you already have a Python interpreter you want to use:

```r
rdoc_setup(python = "/path/to/python")
```

If you prefer to skip explicit setup, `rdoc_decode()` can do it automatically on first use:

```r
res <- rdoc_decode(fs_overlay = your_overlay)
```

## ⚡ Quick Workflow

1. Decode your overlay with `rdoc_decode()`.
2. Plot one result with `rdoc_circleplot()`.
3. Compare multiple results with `rdoc_compare_heatplot()` or `rdoc_compare_fanplot()`.

### Minimal Example

```r
library(rdocodeR)

terms <- readRDS(rdoc_terms_file())
overlay <- terms[[1]]

res <- rdoc_decode(fs_overlay = overlay)
head(res)
```

### 1) Decode an Overlay

`rdoc_decode()` uses cached term-specific eigenstrap null maps.
These null maps are **not** shipped precomputed with the package. Instead, on first use,
`rdocodeR` can generate and cache them locally through `rdoc_setup()`.
The default setup generates `1000` nulls per term, saves them in the user cache directory,
and then reuses them automatically on later runs.

```r
library(rdocodeR)

rdoc_setup()   # optional pre-warm; defaults to 1000 nulls

res <- rdoc_decode(
  fs_overlay = your_overlay,
  cor_method = "pearson"   # default
)

# Optional TSV export
rdoc_decode(
  fs_overlay = your_overlay,
  save_results = TRUE,
  results_file = "~/Desktop/rdoc_decode_results.tsv"
)
```

### 2) Plot a Single RDoC Decoding Table

```r
df <- rdoc_example_data()

p_circle <- rdoc_circleplot(
  corr_df = df,
  domain_palette = "Accent",
  show_term_labels = TRUE,
  highlight_significant_terms = TRUE,
  correlation_label = "pearson"
)
p_circle
```

![rdoc circleplot example](man/figures/readme-circleplot.png)

### 3) Compare Multiple Decoding Tables (n >= 2)

```r
df1 <- rdoc_example_data()
df2 <- df1
df3 <- df1
set.seed(1)
df2$r <- pmax(pmin(df2$r + rnorm(nrow(df2), sd = 0.10), 1), -1)
df3$r <- pmax(pmin(df3$r + rnorm(nrow(df3), sd = 0.15), 1), -1)

p_heat <- rdoc_compare_heatplot(
  corr_list = list(Sample_A = df1, Sample_B = df2, Sample_C = df3),
  domain_palette = "Accent",
  show_significance_stars = TRUE,
  show_significant_term_labels = TRUE,
  correlation_label = "pearson"
)
p_heat

p_fan <- rdoc_compare_fanplot(
  corr_list = list(Sample_A = df1, Sample_B = df2, Sample_C = df3),
  domain_palette = "Accent",
  show_significance_stars = TRUE,
  correlation_label = "pearson"
)
p_fan
```

![rdoc heatplot example](man/figures/readme-heatplot.png)

![rdoc fanplot example](man/figures/readme-fanplot.png)

### Standalone Helpers

- `rdoc_available_palettes()`
- `plot_rdoc_legend()`
- `plot_rdoc_heatmap_legend()`
- `rdoc_terms_reference()`, `rdoc_terms_file()`, `rdoc_term_nulls_dir()`, and `rdoc_setup()`

## 🔗 References

### Core framework and NIMH documentation

- Insel, T. R., Cuthbert, B. N., Garvey, M. A., Heinssen, R. K., Pine, D. S., Quinn, K. J., Sanislow, C. A., & Wang, P. S. (2010). *Research domain criteria (RDoC): Toward a new classification framework for research on mental disorders.* *American Journal of Psychiatry, 167*(7), 748-751. https://doi.org/10.1176/appi.ajp.2010.09091379 (PubMed: https://pubmed.ncbi.nlm.nih.gov/20595427/)
- National Institute of Mental Health. (n.d.). About RDoC. https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/about-rdoc
- National Institute of Mental Health. (n.d.). RDoC Matrix. https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/constructs/rdoc-matrix
- National Institute of Mental Health. (n.d.). Definitions of the RDoC domains and constructs. https://www.nimh.nih.gov/research/research-funded-by-nimh/rdoc/definitions-of-the-rdoc-domains-and-constructs
- Boehringer Ingelheim. (n.d.). RDoC framework explainer. https://pro.boehringer-ingelheim.com/connecting-psychiatry/news-perspectives/rdoc-framework-explainer

### Further reading

- Cuthbert, B. N. (2014). *The RDoC framework: Facilitating transition from ICD/DSM to dimensional approaches that integrate neuroscience and psychopathology.* *World Psychiatry, 13*(1), 28-35. https://doi.org/10.1002/wps.20087 (PubMed: https://pubmed.ncbi.nlm.nih.gov/24497240/)
- Kozak, M. J., & Cuthbert, B. N. (2016). *The NIMH Research Domain Criteria Initiative: Background, Issues, and Pragmatics.* *Psychophysiology, 53*(3), 286-297. https://doi.org/10.1111/psyp.12518 (PubMed: https://pubmed.ncbi.nlm.nih.gov/26877115/)
- Casey, B. J., Oliveri, M. E., & Insel, T. (2014). *A neurodevelopmental perspective on the research domain criteria (RDoC) framework.* *Biological Psychiatry, 76*(5), 350-353. https://doi.org/10.1016/j.biopsych.2014.01.006 (PubMed: https://pubmed.ncbi.nlm.nih.gov/25103538/)
- Hyman, S. E. (2010). *The diagnosis of mental disorders: The problem of reification.* *Annual Review of Clinical Psychology, 6*, 155-179. https://doi.org/10.1146/annurev.clinpsy.3.022806.091532 (PubMed: https://pubmed.ncbi.nlm.nih.gov/17716032/)

### 📦 BibTeX

BibTeX entries are provided in: [`inst/references/rdoc_references.bib`](inst/references/rdoc_references.bib)
