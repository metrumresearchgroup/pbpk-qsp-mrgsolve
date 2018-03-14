A multiscale systems model of bone health and mineral homeostasis
================
Metrum Research Group, LLC

Background and Motivation
=========================

The model was originally developed to describe the bone marker changes associated with denosumab administration from a then ongoing clinical trial. Associated changes in serum calcium and PTH were also considered of interest at the time and so justified the development of a 'systems' model that included bone remodeling and bone mineral (calcium and phosphate) homeostatic mechanisms. Other therapeutics (e.g., teriparatide) and disease states (kidney failure, parathyroid-related abnormalities) were also considered at the time to further inform the parameterization and estimation of the model (Peterson and Riggs, Bone 2010).

Install from GitHub
===================

``` r
remotes::install_github("metrumresearchgroup/OpenBoneMin")
```

``` r
library(OpenBoneMin)
library(dplyr)
```

Load the Bone / Mineral model
=============================

``` r
mod <- BoneMin()
```

An example simulation of teriparatide administration
====================================================

-   We'll give either 20 or 40 micrograms SQ daily for

``` r
out <- sim_teri(dose = c(20,40), dur = 9)

out
```

    . Model:  OpenBoneMin 
    . Dim:    4804 x 4 
    . Time:   0 to 240 
    . ID:     2 
    .      ID time PTHpm   CaC
    . [1,]  1  0.0  3.85 2.350
    . [2,]  1  0.0  3.85 2.350
    . [3,]  1  0.1 19.59 2.351
    . [4,]  1  0.2 26.29 2.352
    . [5,]  1  0.3 28.62 2.353
    . [6,]  1  0.4 28.87 2.354
    . [7,]  1  0.5 28.17 2.355
    . [8,]  1  0.6 27.04 2.356

PTH profiles for the 20 and 40 microgram doses

``` r
plot(out, PTHpm ~ time|ID, scales = "same")
```

![](img/OpenBoneMin-unnamed-chunk-6-1.png)

Calcium profiles for the 20 and 40 microgram doses

``` r
plot(out, CaC~time)
```

![](img/OpenBoneMin-unnamed-chunk-7-1.png)

Denosumab administration
========================

``` r
out <- sim_denos(dose = c(30,60,210))

plot(out, DENCP ~.,  scales = list(y=list(log = TRUE)), 
     ylim = c(1E-4, 10E5))
```

![](img/OpenBoneMin-unnamed-chunk-8-1.png)

``` r
plot(out, BMDlsDENchange ~ .)
```

![](img/OpenBoneMin-unnamed-chunk-9-1.png)
