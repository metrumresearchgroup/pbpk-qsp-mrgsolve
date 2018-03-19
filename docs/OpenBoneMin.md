A multiscale systems model of bone health and mineral homeostasis
================
Metrum Research Group, LLC

-   [Background and Motivation](#background-and-motivation)
-   [Install from GitHub](#install-from-github)
-   [Load the Bone / Mineral model](#load-the-bone-mineral-model)
-   [Teriparatide example](#teriparatide-example)
-   [Denosumab example](#denosumab-example)
-   [Peterson and Riggs 2012](#peterson-and-riggs-2012)
    -   [Generate the dosing regimens](#generate-the-dosing-regimens)
    -   [Publication figure 3: CTx vs time](#publication-figure-3-ctx-vs-time)
    -   [Publication figure 4: BSAP vs time](#publication-figure-4-bsap-vs-time)
    -   [Publication figure 5: LS BMD vs time](#publication-figure-5-ls-bmd-vs-time)
    -   [Publication figure 6: TGF-*β* vs time](#publication-figure-6-tgf-beta-vs-time)

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
library(tidyverse)
```

Load the Bone / Mineral model
=============================

``` r
mod <- BoneMin()
```

Teriparatide example
====================

-   Teriparatide is a recombinant parathyroid (PTH) hormone, containing only the first 84 amino acid residues
-   Teriparatide doses are administered subcutaneously and are absorbed into the PTH compartment; so teriparatide is considered equivalent to PTH
-   We'll give either 20 or 40 micrograms SQ daily for 10 doses

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

Denosumab example
=================

-   Denosumab is a human monoclonal antibody that binds RANK ligand (RANKL), preventing the formation of RANK/RANKL complex and subsequent simulation of osteoclast maturation
-   The labeled dose for treatment of osteoporosis is 60 mg SQ every 6 months

``` r
out <- sim_denos(dose = c(14,30,60,210))
```

Plot the denosumab concentration versus time

``` r
plot(out, DENCP ~.,  scales = list(y=list(log = TRUE)), ylim = c(1E-4, 10E5))
```

![](img/OpenBoneMin-unnamed-chunk-9-1.png)

Plot the percent change in lumbar spine BMD versus time

``` r
plot(out, BMDlsDENchange ~ .)
```

![](img/OpenBoneMin-unnamed-chunk-10-1.png)

Peterson and Riggs 2012
=======================

**Predicting Nonlinear Changes in Bone Mineral Density Over Time Using a Multiscale Systems Pharmacology Model**

Peterson MC, Riggs MM. CPT Pharmacometrics Syst Pharmacol. 2012 Nov 14;1:e14. doi: 10.1038/psp.2012.15. PubMed PMID: 23835796; PubMed Central PMCID: PMC3600731 [PubMed Citation](https://www.ncbi.nlm.nih.gov/pubmed/23835796) / [Free Full Text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3600731)

-   Bone is constantly being remodeled, constantly turning over. This remodeling process consists of breaking down old bone and building up new bone
    -   *Osteoclasts* are cells that function to disassemble bone (resorption) and *osteoblasts* are cells that function to synthesize or build new bone
    -   One therapeutic strategy for osteoporosis (porous bone) is to inhibit osteoclast activity and thus inhibit the bone resorption activity
    -   **Denosumab** decreases osteoclast activity by preventing maturation of pre-osteoclasts into fully-functional osteoclast cells
-   The goal of the modeling was to predict nonlinear changes in lumbar spine (LS) bone mineral density (BMD) over 48 months of denosumab treatment under different dosing regimens.

-   **LSBMD** was modeled as a function of bone turnover markers
    -   CTx: C-telopeptide, a marker for osteoclast function
    -   BSAP: bone-specific alkaline phosphatase, a marker for osteoblast function
-   Our **simulation objective** is to recreate figures 3, 4, 5, and 6 from the publication. This will show the changes in CTx and BSAP and the subsequent changes in BMD under the model.

Generate the dosing regimens
----------------------------

**60 mg every 6 months x 8**

``` r
month <- 24*28
regi <- function(x) factor(x,labels = c("60 mg Q6M", "14 mg Q6M", 
                                        "30 mg Q3M", "210 mg Q6M"))
xscale <- scale_x_continuous(breaks = c(0,3,6,12,18,24,36,48), limits = c(0,48))
```

``` r
e1 <- ev(amt = 60, ii = 6*month, addl = 7)
e1
```

    . Events:
    .   time cmt amt   ii addl evid
    . 1    0   1  60 4032    7    1

**14 mg every 6 months x4, then 60 mg every 6 months x4**

``` r
e2 <- 
  mutate(e1, amt = 14, addl = 3) %then% 
  mutate(e1, addl = 3)
e2
```

    . Events:
    .    time cmt amt   ii addl evid
    . 1     0   1  14 4032    3    1
    . 2 16128   1  60 4032    3    1

**30 mg every 3 months for 8 doses and changed to 60 mg Q6M starting on month 36**

``` r
e3a <- ev(amt = 30, ii = 3*month, addl = 7)
e3b <- mutate(e1, addl = 1)
e3 <- seq(e3a, wait = 12*month, e3b)
e3
```

    . Events:
    .    time cmt amt   ii addl evid
    . 1     0   1  30 2016    7    1
    . 2 24192   1  60 4032    1    1

**210 mg every 6 months x4 then DC**

``` r
e4 <- mutate(e1, amt = 210, addl = 3)
e4
```

    . Events:
    .   time cmt amt   ii addl evid
    . 1    0   1 210 4032    3    1

**Create one single data set from which to simulate**

``` r
data <- as_data_set(e1,e2,e3,e4)

data
```

    .   ID  time cmt evid amt   ii addl
    . 1  1     0   1    1  60 4032    7
    . 2  2     0   1    1  14 4032    3
    . 3  2 16128   1    1  60 4032    3
    . 4  3     0   1    1  30 2016    7
    . 5  3 24192   1    1  60 4032    1
    . 6  4     0   1    1 210 4032    3

``` r
out <- 
  mod %>%
  mrgsim(data = data, end = 48*month, delta = 24) %>% 
  mutate(ID = regi(ID), Month = time/month)
```

Publication figure 3: CTx vs time
---------------------------------

Shown here is the percent change from baseline value for osteoclasts (via CTx)

``` r
ggplot(out) + 
  geom_line(aes(x = Month, y = OCchange, col = factor(ID)), lwd = 1) + 
  facet_grid(~ID) + geom_hline(yintercept = 100, lty = 2) +
  scale_y_continuous(trans = "log10", breaks = c(10,30,100,300), limits = c(5,300)) + 
  theme_bw() + theme(legend.position = "top") +
   geom_vline(xintercept = c(24,36), lty = 3) + xscale
```

![](img/OpenBoneMin-unnamed-chunk-18-1.png)

Publication figure 4: BSAP vs time
----------------------------------

Shown here is the percent change from baseline value for osteoblasts (via BSAP)

``` r
ggplot(out) + 
  geom_line(aes(x = Month, y = OBchange, col = factor(ID)), lwd = 1) + 
  facet_grid(~ID) + geom_hline(yintercept = 100, lty = 2) +
  scale_y_continuous(trans = "log10", breaks = c(10,30,100,300), limits = c(5,300)) +
  theme_bw() + theme(legend.position = "top") + 
  geom_vline(xintercept = c(24,36), lty = 3) + xscale
```

![](img/OpenBoneMin-unnamed-chunk-19-1.png)

Publication figure 5: LS BMD vs time
------------------------------------

Shown here is the percent change from baseline value for lumbar spine BMD

``` r
ggplot(out) + 
  geom_line(aes(x = Month, y = BMDlsDENchange, col = factor(ID)), lwd = 1) + 
  facet_grid(~ID) + 
  theme_bw() + theme(legend.position = "top") + 
  geom_vline(xintercept = c(24,36), lty = 3) + xscale
```

![](img/OpenBoneMin-unnamed-chunk-20-1.png)

Publication figure 6: TGF-*β* vs time
-------------------------------------

-   Solid lines are **active** TGF-*β*
-   Dashed lines are **latent** TGF-*β*

``` r
out <- 
  out %>%
  group_by(ID) %>%
  mutate(ATGF = 100*TGFBact/first(TGFBact), LTGF = 100*TGFB/first(TGFB)) %>%
  ungroup

ggplot(out) + 
  geom_line(aes(x = Month, y = ATGF, col = factor(ID)), lwd = 1) + 
  geom_line(aes(x = Month, y = LTGF), col = "black", lty = 2, lwd = 0.6) + 
  facet_grid(~ID) + 
  theme_bw() + theme(legend.position = "top") + 
  geom_vline(xintercept = c(24,36), lty = 3) + xscale +
  scale_y_continuous(breaks = c(0,25,50,75,100,150,225), limits = c(0,225)) 
```

![](img/OpenBoneMin-unnamed-chunk-21-1.png)
