Clinical responses to ERK inhibition in BRAF{V600E}-mutant colorectal
cancer
================
Metrum Research Group

- [Reference](#reference)
- [Introduction](#introduction)
  - [Cast of characters](#cast-of-characters)
  - [Translation](#translation)
- [Set up](#set-up)
- [Explore](#explore)
  - [Simulate with ERK inhibitor
    GDC-0944](#simulate-with-erk-inhibitor-gdc-0944)
  - [Sensitivity analysis](#sensitivity-analysis)
- [Predicting clinical outcomes for combination
  therapies](#predicting-clinical-outcomes-for-combination-therapies)
  - [Generate dosing regimens](#generate-dosing-regimens)
  - [Simulate all combination
    therapies](#simulate-all-combination-therapies)
  - [Summarize and plot](#summarize-and-plot)
- [Target populations more likely to
  respond](#target-populations-more-likely-to-respond)
  - [ORR in full population: GDC +/-
    COBI](#orr-in-full-population-gdc---cobi)
  - [ORR in select patients: GDC +/-
    COBI](#orr-in-select-patients-gdc---cobi)

# Reference

**Clinical responses to ERK inhibition in BRAF{V600E}-mutant colorectal
cancer predicted using a computational model**

- Daniel C. Kirouac, Gabriele Schaefer, Jocelyn Chan, Mark Merchant,
  Christine Orr, Shih-Min A. Huang, John Moffat, Lichuan Liu, Kapil
  Gadkar and Saroja Ramanujan

- npj Systems Biology and Applications (2017) 3:14 ; <doi:10.1038> /
  s41540-017-0016-1

# Introduction

(Summarized from Introduction in the reference)

- The V600E/K mutation results in constitutively active BRAF, with
  subsequent signalling through MEK and ERK

- BRAF and MEK inhibitors were found to be effective in V600E mutant
  melanoma, but not so much in colorectal cancer

- Could resistance to BRAF inhibitors be mediated through EGFR
  signalling through RAS and CRAF?

- What about inhibition at ERK?

- Could the effectiveness of different combination therapies be
  predicted with a model characterizing this biology?

## Cast of characters

- **vemurafenib**: BRAF inhibitor (selective for V600E mutant)
- **cobimetinib**: MEK inhibitor
- **cetuximab**: EGFR antibody
- **GDC-0994**: ERK inhibitor (the star)

## Translation

- Model published as SBML
- Translator from previous project work using R bindings to libSBML
- Minor modifications to the translator code to accommodate the MAPK
  model as published

# Set up

``` r
library(mrgsolve)
library(tidyverse)
library(parallel)
library(mrgsim.sa)
library(future.apply)
library(knitr)
source(here("docs/functions.R"))
```

Read in the virtual population

``` r
vp <- readRDS(here("docs/data/s10vpop_pk.RDS")) %>% mutate(VPOP2 = seq(n()))

dim(vp)
```

    . [1] 250 147

Load the model and pick one parameter set from vpop

``` r
mod <- mread("mapk", here("docs/models"), soloc = here("docs/build"), end = 56) 

mod <- param(mod, filter(vp, VPOP2==41))
```

# Explore

## Simulate with ERK inhibitor GDC-0944

``` r
e <- expand.ev(amt = seq(100,600,100), cmt = 12, ii = 1, addl = 20) %>% as.ev()

e <- ev_seq(e, wait = 7, e) %>% as_tibble() %>% arrange(ID)

mod %>% 
  data_set(e) %>%
  mrgsim(delta = 0.25) %>% 
  plot(ERKi + TUMOR ~ time)
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/mapk-unnamed-chunk-5-1.png)<!-- -->

## Sensitivity analysis

The authors note two parameters that are particularly influential with
respect to response rates:

- `wOR`: MAPK pathway dependence parameter
- $\delta_{max}$: the maximum cell death rate

``` r
vp %>% 
  select(wOR,dmax) %>%
  map(quantile, probs = seq(0,1,0.1)) %>% 
  bind_cols() %>% 
  mutate(pctile = seq(0,1,0.1))
```

    . # A tibble: 11 × 3
    .      wOR   dmax pctile
    .    <dbl>  <dbl>  <dbl>
    .  1 0.754 0.0331    0  
    .  2 0.828 0.0369    0.1
    .  3 0.860 0.0380    0.2
    .  4 0.873 0.0400    0.3
    .  5 0.890 0.0416    0.4
    .  6 0.906 0.0426    0.5
    .  7 0.918 0.0440    0.6
    .  8 0.951 0.0449    0.7
    .  9 0.977 0.0464    0.8
    . 10 1     0.0485    0.9
    . 11 1     0.0522    1

**Sensitivity analysis on MAPK pathway dependence**

- Given the other parameters for this virtual patient, whenever `wOR`
  falls below 90% or so, you get no change in tumor size

- The authors discuss the possibility of looking for biomarkers that
  would indicate a patient is showing high dependence on MAPK signalling
  to improve response rates

``` r
ev400 <- filter(e, amt==400) %>% select(-ID) %>% as.ev()

mod %>% 
  ev(ev400) %>% 
  Req(TUMOR) %>%
  parseq_range(wOR = c(0.9,1.05), .n = 8) %>%
  sens_each() %>% 
  sens_plot("TUMOR")
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/mapk-unnamed-chunk-7-1.png)<!-- -->

**Sensitivity analysis on $\delta_{max}$**

- In this analysis, we look between 0.2 and 2 times the nominal
  parameter value for this virtual patient

- The authors do a simulation with a hypothetical agent that increases
  the cell death rate by 10%

``` r
mod %>% 
  ev(ev400) %>% 
  Req(TUMOR) %>%
  parseq_fct(dmax, .n = 8) %>%
  sens_each() %>% 
  sens_plot("TUMOR")
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/mapk-unnamed-chunk-8-1.png)<!-- -->

# Predicting clinical outcomes for combination therapies

- Re-create figure 6B in the publication

## Generate dosing regimens

- **No treatment**

``` r
data0 <- ev(amt = 0, cmt = 8)
```

- **BRAF inhibitor** - vemurafanib (VEMU)
- Compartment 8

``` r
dataV <- ev(amt = 960, cmt = 8, ii = 0.5, addl = 120)
```

- **ERK inhibitor** - GCD-994 (GDC)
- Compartment 12

``` r
dataG <- ev(amt = 400, cmt = 12, ii = 1, addl = 20)

dataG <- seq(dataG, wait = 7, dataG) 

out <- mrgsim(mod, ev = dataG, end = 56, delta = 0.1)

plot(out, ERKi_C~time)
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/mapk-unnamed-chunk-11-1.png)<!-- -->

- **MEK inhibitor** - cobimetinib (COBI)
- Compartment 10

``` r
dataCO <- mutate(dataG, amt = 60, cmt = 10)
```

- **EGFR inihbitor** - cetuximab (CETUX)
- Compartment 7

``` r
dataCE <- ev(cmt = 7, ii = 7, addl = 7, amt = 450)
```

We create two functions: one to combine dosing regimens and the other to
simulate from a dosing regimen

``` r
comb <- function(...) {
  x <- lapply(list(...), as.data.frame)
  bind_rows(x) %>% arrange(time)
}

sim <- function(d, m, v) {
  m %>%
    ev(as.ev(d)) %>%
    mrgsim(idata = v, end = -1, add = 56) %>%
    filter(time==56) 
}
```

For example, to dose cetuximab with vemurafanib, we’d do

``` r
comb(dataCE, dataV)
```

    .   time amt  ii addl cmt evid
    . 1    0 450 7.0    7   7    1
    . 2    0 960 0.5  120   8    1

``` r
sim(d = comb(dataCE,dataV), m = mod, v = slice(vp, seq(10)))
```

    . # A tibble: 10 × 25
    .       ID  time   TD1 CELLS      FB1      FB2      FB3   FB4 RTK1i_blood RAFi_gut
    .    <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>    <dbl> <dbl>       <dbl>    <dbl>
    .  1     1    56 0.996 0.871 0.695    0.695    0.695    0.437        246.    1023.
    .  2     2    56 0.401 0.963 0.000144 0.000144 0.000144 0.547        246.    1316.
    .  3     3    56 0.363 1.18  0.000115 0.000115 0.000115 0.703        246.    1540.
    .  4     4    56 0.186 0.726 0.000152 0.000152 0.000151 0.649        246.    1926.
    .  5     5    56 0.348 0.903 0.0221   0.0220   0.0220   0.662        246.    1027.
    .  6     6    56 0.999 1.06  0.939    0.939    0.939    0.731        246.     960.
    .  7     7    56 0.999 1.42  0.902    0.902    0.902    0.731        246.     960.
    .  8     8    56 0.636 1.02  0.000239 0.000230 0.000227 0.630        246.    1590.
    .  9     9    56 0.952 1.16  0.195    0.195    0.195    0.688        246.     994.
    . 10    10    56 0.964 1.02  0.252    0.252    0.252    0.649        246.    5437.
    . # ℹ 15 more variables: RAFi_blood <dbl>, MEKi_gut <dbl>, MEKi_blood <dbl>,
    . #   ERKi_gut <dbl>, ERKi_blood <dbl>, AKTi_gut <dbl>, AKTi_blood <dbl>,
    . #   MEKi_V3 <dbl>, RTK1i_gut <dbl>, ERKi <dbl>, ERKi_C <dbl>, RAFi <dbl>,
    . #   MEKi <dbl>, TUMOR <dbl>, GDC <dbl>

## Simulate all combination therapies

Generate a data frame of runs to do

``` r
sims <- 
  tribble(
    ~label, ~object, 
    "No Treatment",        data0,
    "CETUX",               dataCE, 
    "VEMU",                dataV,
    "COBI",                dataCO, 
    "GDC",                 dataG,
    "CETUX+VEMU",          comb(dataCE, dataV), 
    "CETUX+COBI",          comb(dataCE, dataCO), 
    "CETUX+GDC",           comb(dataCE, dataG),
    "VEMU+COBI",           comb(dataV, dataG), 
    "VEMU+GDC",            comb(dataV, dataG),
    "COBI+GDC",            comb(dataCO, dataG),
    "CETUX+VEMU+COBI",     comb(dataCE, dataV,dataCO), 
    "CETUX+VEMU+GDC",      comb(dataCE, dataV,dataG), 
    "CETUX+COBI+GDC",      comb(dataCE, dataCO,dataG), 
    "VEMU+COBI+GDC",       comb(dataV, dataCO,dataG),
    "CETUX+VEMU+COBI+GDC", comb(dataCE, dataV, dataCO, dataG)
  ) %>% mutate(object = map(object, as.data.frame))
```

Run the simulation

``` r
plan(multisession, workers = 5L)

sims <- mutate(
  sims, 
  out = future_lapply(X = object, FUN = sim, v = vp, m = mod, future.seed = TRUE),
)
```

## Summarize and plot

Get ready to plot

``` r
sms <- select(sims, label, out) %>% unnest(cols = c(out))

sms <- mutate(
  sms, 
  labelf = fct_inorder(label), 
  gdc = factor(grepl("GDC", label))
)

head(sms)
```

    . # A tibble: 6 × 28
    .   label        ID  time   TD1 CELLS   FB1   FB2   FB3   FB4 RTK1i_blood RAFi_gut
    .   <chr>     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>       <dbl>    <dbl>
    . 1 No Treat…     1    56 0.998 0.872 0.822 0.822 0.822 0.714           0        0
    . 2 No Treat…     2    56 0.998 1.42  0.800 0.800 0.800 0.478           0        0
    . 3 No Treat…     3    56 0.999 1.78  0.916 0.916 0.916 0.708           0        0
    . 4 No Treat…     4    56 0.995 1.51  0.583 0.583 0.583 0.644           0        0
    . 5 No Treat…     5    56 0.999 1.43  0.951 0.951 0.951 0.664           0        0
    . 6 No Treat…     6    56 0.999 1.06  0.940 0.940 0.940 0.731           0        0
    . # ℹ 17 more variables: RAFi_blood <dbl>, MEKi_gut <dbl>, MEKi_blood <dbl>,
    . #   ERKi_gut <dbl>, ERKi_blood <dbl>, AKTi_gut <dbl>, AKTi_blood <dbl>,
    . #   MEKi_V3 <dbl>, RTK1i_gut <dbl>, ERKi <dbl>, ERKi_C <dbl>, RAFi <dbl>,
    . #   MEKi <dbl>, TUMOR <dbl>, GDC <dbl>, labelf <fct>, gdc <fct>

``` r
p1 <- 
  ggplot(data=sms) + 
  geom_point(aes(x=labelf, y=TUMOR),position=position_jitter(width=0.15),col="grey") +
  geom_hline(yintercept=0.7,col="black", lty=1,lwd=0.7)  +
  scale_y_continuous(limits=c(0,2.5),name="Tumor size",breaks=c(0,0.5,1,1.5,2,2.5,3)) +
  scale_x_discrete(name="") + 
  geom_boxplot(aes(x=labelf,y=TUMOR,col=gdc),fill="darkslateblue",alpha=0.2) +
  scale_color_manual(values = c("darkslateblue", "firebrick"), guide = "none") + 
  theme_plain() + rotx(30) + 
  ggtitle("Note: GDC-0944 +/- cobimetinib")

p1
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/mapk-unnamed-chunk-18-1.png)<!-- -->

# Target populations more likely to respond

## ORR in full population: GDC +/- COBI

``` r
sms %>%
  filter(label %in% c("GDC", "COBI+GDC")) %>%
  group_by(label) %>%
  summarise(orr = mean(TUMOR < 0.7)) %>% 
  kable(digits = 3)
```

| label    |   orr |
|:---------|------:|
| COBI+GDC | 0.352 |
| GDC      | 0.140 |

## ORR in select patients: GDC +/- COBI

- We filter down to those with `wOR` greater than the median

``` r
vp_select <- filter(vp, wOR > median(wOR))
```

And resimulate

``` r
plan(sequential)
plan(multisession, workers = 5L)

re_run <- 
  sims %>%
  select(label, object) %>%
  filter(label %in% c("GDC", "COBI+GDC")) %>% 
  mutate(
    out = future_lapply(object, sim, v = vp_select, m = mod, future.seed = TRUE)
  ) %>%
  select(label, out) %>% 
  unnest(cols = c(out))
```

Now, response rates are up

``` r
re_run %>%
  group_by(label) %>%
  summarise(orr = mean(TUMOR < 0.7)) %>% 
  kable(digits = 3)
```

| label    |   orr |
|:---------|------:|
| COBI+GDC | 0.705 |
| GDC      | 0.287 |
