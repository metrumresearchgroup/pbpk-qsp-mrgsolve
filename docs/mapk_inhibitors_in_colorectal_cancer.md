Clinical responses to ERK inhibition in BRAF{V600E}-mutant colorectal cancer
================
Metrum Research Group, LLC

-   [About](#about)
    -   [Cast of characters](#cast-of-characters)
-   [Translation](#translation)
-   [Explore](#explore)
-   [Combination therapies](#combination-therapies)
    -   [Generate dosing regimens](#generate-dosing-regimens)

``` r
library(mrgsolve)
library(tidyverse)
source("functions.R")
```

**Clinical responses to ERK inhibition in BRAF{V600E}-mutant colorectal cancer predicted using a computational model**

-   Daniel C. Kirouac, Gabriele Schaefer, Jocelyn Chan, Mark Merchant, Christine Orr, Shih-Min A. Huang, John Moffat, Lichuan Liu, Kapil Gadkar and Saroja Ramanujan

-   npj Systems Biology and Applications (2017) 3:14 ; <doi:10.1038> / s41540-017-0016-1

About
=====

(Summarized from Introduction in the reference)

-   The V600E/K mutation results in constitutively active BRAF, with subsequent signalling through MEK and ERK

-   BRAF and MEK inhibitors were found to be effective in V600E mutant melanoma, but not so much in colorectal cancer
    -   Could resistance to BRAF inhibitors be mediated through EGFR signalling through RAS and CRAF?
    -   What about inhibition at ERK?
-   Could the effectiveness of different combination therapies be predicted with a model characterizing this biology?

Cast of characters
------------------

-   **vemurafenib**: BFAF inhibitor (selective for V600E mutant)
-   **cobimetinib**: MEK inhibitor
-   **cetuximab**: EGFR antibody
-   **GDC-0994**: ERK inhibitor (the star)

Translation
===========

-   Model published as SBML
-   Translator from previous project work using R bindings to libSBML
-   Minor modifications to the translator code to accommodate the MAPK model as published

Explore
=======

Read in the virtual population

``` r
vp <- readRDS("data/s10vpop_pk.RDS") %>% mutate(VPOP2 = seq(n()))

dim(vp)
```

    . [1] 250 147

Load the model and pick one parameter set from vpop

``` r
mod <- mread_cache("mapk", "models") %>% update(end  = 56)

mod <- param(mod, filter(vp,VPOP2==24))
```

``` r
e <- expand.ev(amt = c(1,3,10,30,60,100), cmt = 10, ii = 1, addl = 20)

e <- ev_seq(e, wait = 7, e) %>% as_data_frame %>%  arrange(ID)

mod %>% 
  data_set(e) %>%
  Req(TUMOR,MEKi) %>%
  mrgsim(delta = 0.25) %>% 
  plot(MEKi+TUMOR ~ time)
```

![](img/mapk-unnamed-chunk-5-1.png)

-   Look at MAPK dependence

``` r
library(mrgsolvetk)

dat <- filter(e, amt==60) %>% select(-ID) %>% as.ev

mod %>% 
  ev(as.ev(dat)) %>% Req(TUMOR) %>%
  sens_spaced(wOR = c(0.9,0.99), .n = 10) %>%
  ggplot(aes(time, TUMOR, col = factor(value), group = value)) + 
  geom_line(lwd = 1) + geom_hline(yintercept = 1, lty = 2) + theme_bw()
```

![](img/mapk-unnamed-chunk-6-1.png)

-   What about *δ*<sub>*m**a**x*</sub> dependence?

``` r
mod %>% 
  ev(as.ev(dat)) %>% Req(TUMOR) %>%
  select(dmax) %>% 
  sens_spaced_factor(.factor = 2, .n = 10) %>%
  mutate(value = signif(value,3)) %>%
  ggplot(aes(time, TUMOR, col = factor(value), group = value)) + 
  geom_line(lwd = 1) + geom_hline(yintercept = 1, lty = 2) + theme_bw()
```

![](img/mapk-unnamed-chunk-7-1.png)

Combination therapies
=====================

-   Re-create figure 6B in the publication

Update the model object with one set of parameters

``` r
mod <- param(mod, slice(vp, 1))
```

``` r
summary(vp$wOR)
```

    .    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    .  0.7540  0.8678  0.9060  0.9108  0.9694  1.0000

``` r
summary(vp$dmax)
```

    .    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    . 0.03306 0.03939 0.04257 0.04264 0.04629 0.05218

Generate dosing regimens
------------------------

-   **No treatment**

``` r
data0 <- ev(amt=0, cmt=8)
```

-   **BRAF inhibitor** - vemurafanib (VEMU)
-   Compartment 8

``` r
dataV <- ev(amt=960,  cmt=8, ii=0.5, addl=120)
```

-   **ERK inhibitor** - GCD-994 (GDC)
-   Compartment 12

``` r
dataG <- ev(amt = 400, cmt = 12, ii = 1, addl = 20)

dataG <- seq(dataG, wait = 7, dataG) 

out <- mrgsim(mod, ev=dataG, end=56)

plot(out, ERKi_C~time)
```

![](img/mapk-unnamed-chunk-12-1.png)

-   **MEK inhibitor** - cobimetinib (COBI)
-   Compartment 10

``` r
dataCO <- mutate(dataG,amt=60,cmt=10)
```

-   **EGFR inihbitor** - cetuximab (CETUX)
-   Compartment 7

``` r
dataCE <- ev(cmt=7,ii=7,addl=7,amt=450)
```

We create two functions: one to combine dosing regimens and the other to simulate from a dosing regimen

``` r
comb <- function(...) {
  x <- lapply(list(...), as.data.frame)
  bind_rows(x) %>% arrange(time)
}

sim <- function(Data,Vp,Mod) {
  Mod %>%
    ev(as.ev(Data)) %>%
    mrgsim(idata=Vp,end=-1, add = 56) %>%
    filter(time==56) 
}
```

For example

``` r
comb(dataCE,dataV)
```

    .   time cmt  ii addl amt evid
    . 1    0   7 7.0    7 450    1
    . 2    0   8 0.5  120 960    1

``` r
sim(comb(dataCE,dataV), Vp = slice(vp,seq(10)), Mod = mod)
```

    . # A tibble: 10 x 25
    .       ID  time   TD1 CELLS      FB1      FB2      FB3   FB4 RTK1i_blood
    .    <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>    <dbl> <dbl>       <dbl>
    .  1    1.   56. 0.996 0.871 0.695    0.695    0.695    0.437        246.
    .  2    2.   56. 0.401 0.963 0.000144 0.000144 0.000144 0.547        246.
    .  3    3.   56. 0.363 1.18  0.000115 0.000115 0.000115 0.703        246.
    .  4    4.   56. 0.186 0.726 0.000152 0.000152 0.000151 0.649        246.
    .  5    5.   56. 0.348 0.903 0.0221   0.0220   0.0220   0.662        246.
    .  6    6.   56. 0.999 1.06  0.939    0.939    0.939    0.731        246.
    .  7    7.   56. 0.999 1.42  0.902    0.902    0.902    0.731        246.
    .  8    8.   56. 0.636 1.02  0.000239 0.000230 0.000227 0.630        246.
    .  9    9.   56. 0.952 1.16  0.195    0.195    0.195    0.688        246.
    . 10   10.   56. 0.964 1.02  0.252    0.252    0.252    0.649        246.
    . # ... with 16 more variables: RAFi_gut <dbl>, RAFi_blood <dbl>,
    . #   MEKi_gut <dbl>, MEKi_blood <dbl>, ERKi_gut <dbl>, ERKi_blood <dbl>,
    . #   AKTi_gut <dbl>, AKTi_blood <dbl>, MEKi_V3 <dbl>, RTK1i_gut <dbl>,
    . #   ERKi <dbl>, ERKi_C <dbl>, RAFi <dbl>, MEKi <dbl>, TUMOR <dbl>,
    . #   GDC <dbl>

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
  ) %>% mutate(object = map(object,as.data.frame))
```

Run the simulation

``` r
sims <- mutate(sims, out = parallel::mclapply(object,sim, Vp = vp, Mod = mod))
```

Get ready to plot

``` r
sms <- select(sims, label, out) %>% unnest()

sms <- mutate(sms, labelf = fct_inorder(label))

head(sms)
```

    . # A tibble: 6 x 27
    .   label           ID  time   TD1 CELLS   FB1   FB2   FB3   FB4 RTK1i_blood
    .   <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>       <dbl>
    . 1 No Treatment    1.   56. 0.998 0.872 0.822 0.822 0.822 0.714          0.
    . 2 No Treatment    2.   56. 0.998 1.42  0.800 0.800 0.800 0.478          0.
    . 3 No Treatment    3.   56. 0.999 1.78  0.916 0.916 0.916 0.708          0.
    . 4 No Treatment    4.   56. 0.995 1.51  0.583 0.583 0.583 0.644          0.
    . 5 No Treatment    5.   56. 0.999 1.43  0.951 0.951 0.951 0.664          0.
    . 6 No Treatment    6.   56. 0.999 1.06  0.940 0.940 0.940 0.731          0.
    . # ... with 17 more variables: RAFi_gut <dbl>, RAFi_blood <dbl>,
    . #   MEKi_gut <dbl>, MEKi_blood <dbl>, ERKi_gut <dbl>, ERKi_blood <dbl>,
    . #   AKTi_gut <dbl>, AKTi_blood <dbl>, MEKi_V3 <dbl>, RTK1i_gut <dbl>,
    . #   ERKi <dbl>, ERKi_C <dbl>, RAFi <dbl>, MEKi <dbl>, TUMOR <dbl>,
    . #   GDC <dbl>, labelf <fct>

``` r
p1 <- 
  ggplot(data=sms) + 
  geom_point(aes(x=labelf, y=TUMOR),position=position_jitter(width=0.15),col="grey") +
  scale_y_continuous(limits=c(0,2.5),name="Tumor size",breaks=c(0,0.5,1,1.5,2,2.5,3)) +
  scale_x_discrete(name="") + 
  geom_hline(yintercept=0.7,col="firebrick", lty=1,lwd=1)  +
  geom_boxplot(aes(x=labelf,y=TUMOR),fill="darkslateblue",col="darkslateblue",alpha=0.2) +
  theme_plain() + rotx(30)
p1
```

![](img/mapk-unnamed-chunk-20-1.png)
