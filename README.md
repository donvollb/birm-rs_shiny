# birm-rs_shiny

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![R](https://img.shields.io/badge/R-Shiny-276DC3)

A Shiny app that illustrates the **Beta Item Response Model with Response Styles (BIRM-RS)**. The BIRM-RS extends the Beta Item Response Model (BIRM; Noel & Dauvier, 2007) for responses on visual analogue scales with parameters for extreme and acquiescent response styles (Vollbracht et al., 2026). The app shows how the person, item, and response style parameters affect the distribution of responses.

**App:** https://donvollb.github.io/birm-rs_shiny/

The app was developed as a supplement to the following paper:

> Vollbracht, D., Lischetzke, T., & Henninger, M. (2026). *Detecting response styles in visual analogue scales* [Preprint]. PsyArXiv. https://doi.org/10.31234/osf.io/xteku_v2

The R package [birmrssim](https://github.com/donvollb/birmrssim) provides functions to simulate data from the BIRM-RS and fit the model via Stan.

---

## Tabs

| Tab | Description |
|---|---|
| Standard Beta | Beta distribution with shape parameters $m_{ij}$ and $n_{ij}$ |
| BIRM | Original model with trait, item difficulty, and item dispersion |
| ERS only | Adds the extreme response style parameter |
| ARS only | Adds the acquiescence response style parameter and the item direction |
| BIRM-RS | Combines both response styles |
| BIRM-RS Simulation | Simulated responses of a sample of persons to one item |

In the tabs with response styles, a grey dashed curve shows the same person (or persons) without response styles.

---

## Running the App Locally

The app runs in the browser via [shinylive](https://posit-dev.github.io/r-shinylive/), so no installation is needed. To run it locally in R:

```r
install.packages(c("shiny", "bslib", "shinythemes"))
shiny::runApp("app.R")
```

---

## Repository Structure

| File | Description |
|---|---|
| `app.R` | The Shiny app (single source of the code) |
| `index.qmd` | Web page with explanation; includes `app.R` in a shinylive chunk |
| `docs/` | Rendered web page for GitHub Pages |
| `_extensions/` | Quarto shinylive extension |

If you clone this repository and want to rebuild the web page after changing `app.R` or `index.qmd`, run `quarto render` in the repository folder. This requires [Quarto](https://quarto.org/) and the R package `shinylive`.

---

## References

Noel, Y., & Dauvier, B. (2007). A beta item response model for continuous bounded responses. *Applied Psychological Measurement*, *31*(1), 47–73. https://doi.org/10.1177/0146621605287691

Vollbracht, D., Lischetzke, T., & Henninger, M. (2026). *Detecting response styles in visual analogue scales* [Preprint]. PsyArXiv. https://doi.org/10.31234/osf.io/xteku_v2

---

## License

The code and text in this repository are licensed under the MIT License, see [LICENSE](LICENSE). The folder `docs/index_files/` also contains third-party software bundled by shinylive (e.g., webR and R packages), which remains under its own licenses.
