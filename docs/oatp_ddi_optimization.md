Fit parameters in PBPK model
================
Metrum Research Group

- [Packages and setup](#packages-and-setup)
- [Reference](#reference)
- [Data](#data)
- [PBPK model: pitavastatin / CsA
  DDI](#pbpk-model-pitavastatin--csa-ddi)
- [Objective function](#objective-function)
- [Optimize with different methods](#optimize-with-different-methods)
  - [`minqa::newuoa`: minimization without
    derivatives](#minqanewuoa-minimization-without-derivatives)
  - [`optim`: Nelder-Mead](#optim-nelder-mead)
  - [`DEoptim`: differential evolution
    algorithm](#deoptim-differential-evolution-algorithm)
  - [`GenSA`: simulated annealing](#gensa-simulated-annealing)
  - [`hydroPSO`: particle swarm
    optimization](#hydropso-particle-swarm-optimization)
- [Compare optimization methods](#compare-optimization-methods)

# Packages and setup

``` r
library(tidyverse)
library(mrgsolve)
library(pmxTools)
library(minqa)
library(RcppDE)
library(GenSA)
library(hydroPSO)
library(here)
library(knitr)

source(here("docs/functions.R"))

set.seed(10101)

theme_set(theme_bw() + theme(legend.position="top"))

scale_colour_discrete <- function(...) scale_color_brewer(palette="Dark2")
```

# Reference

**Quantitative Analyses of Hepatic OATP-Mediated Interactions Between
Statins and Inhibitors Using PBPK Modeling With a Parameter Optimization
Method**

- T Yoshikado, K Yoshida, N Kotani, T Nakada, R Asaumi, K Toshimoto, K
  Maeda, H Kusuhara and Y Sugiyama

- CLINICAL PHARMACOLOGY & THERAPEUTICS \| VOLUME 100 NUMBER 5 \|
  NOVEMBER 2016

- <https://www.ncbi.nlm.nih.gov/pubmed/27170342>

# Data

- Example taken from figure 4a from the publication
- Using this as example data to fit

``` r
data <- read_csv(here("docs/data/fig4a-fit.csv")) 

data <- 
  mutate(
    data, 
    type = ID, 
    typef = factor(ID, labels = c("Statin", "Statin+CsA"))
  )

head(data)
```

    . # A tibble: 6 × 8
    .      ID  TIME    DV  EVID   AMT   CMT  type typef     
    .   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     
    . 1     2  0     NA       1  2000     2     2 Statin+CsA
    . 2     2  1     NA       1    30     1     2 Statin+CsA
    . 3     2  1.49  73.7     0     0     0     2 Statin+CsA
    . 4     2  1.99 102.      0     0     0     2 Statin+CsA
    . 5     2  2.49  59.9     0     0     0     2 Statin+CsA
    . 6     2  3.00  37.6     0     0     0     2 Statin+CsA

**Important**: the data has `DV` set to `NA` (missing value); this is
important to for the following approach to work.

- The goal is to fit the pitavastatin data either alone (left) or in
  combination with cyclosporin administered 1 hour before the
  pitavastatin

``` r
ggplot(data = data, aes(TIME, DV)) + 
  geom_point(aes(col = typef), size = 3) + 
  geom_line(col = "darkgrey", aes(group = typef)) + 
  scale_y_continuous(trans = "log", limits = c(0.1,300), breaks = logbr()) 
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/oatp_ddi_optimization-unnamed-chunk-4-1.png)<!-- -->

# PBPK model: pitavastatin / CsA DDI

- Check out the model / data with a quick simulation
- We’re only interested in `CP`, the pitavastatin concentration

``` r
mod <- mread_cache("yoshikado", here("docs/models")) 

mod <- update(mod, end = 14, delta = 0.1, outvars = "CP")
```

This model has 43 parameters and 31 compartments

``` r
param(mod)
```

    . 
    .  Model parameters (N=43):
    .  name       value . name     value   . name  value 
    .  beta       0.8   | ifhCLint 0.00839 | PSmus 3.5   
    .  CLr        0     | ika      0.999   | PSski 0.534 
    .  exFadi     0.145 | ikiu     0.0118  | Qadi  0.223 
    .  exFliv     0.278 | iKp_adi  17.3    | Qh    1.2   
    .  exFmus     0.146 | iKp_liv  16.7    | Qmus  0.642 
    .  exFski     0.321 | iKp_mus  2.98    | Qski  0.257 
    .  fafg       1     | iKp_ski  13.6    | Rdiff 0.0345
    .  fb         0.008 | imw      1203    | tlag  1     
    .  fbCLintall 0.737 | itlag    0.254   | Vadi  0.143 
    .  fbile      0.33  | ka       1.06    | Vcent 0.075 
    .  fh         0.035 | Kp_adi   0.086   | Vliv  0.0241
    .  gamma      0.244 | Kp_mus   0.113   | Vmus  0.429 
    .  iClr       0     | Kp_ski   0.481   | Vski  0.111 
    .  ifafg      0.572 | ktr      0.679   | .     .     
    .  ifb        0.06  | PSadi    0.146   | .     .

``` r
init(mod)
```

    . 
    .  Model initial conditions (N=31):
    .  name       value . name         value . name         value
    .  ac (26)    0     | hc4 (18)     0     | iliv3 (29)   0    
    .  adi (5)    0     | hc5 (19)     0     | iliv4 (30)   0    
    .  ae (23)    0     | he1 (10)     0     | iliv5 (31)   0    
    .  cent (3)   0     | he2 (11)     0     | mc (24)      0    
    .  ehc1 (7)   0     | he3 (12)     0     | me (21)      0    
    .  ehc2 (8)   0     | he4 (13)     0     | mus (4)      0    
    .  ehc3 (9)   0     | he5 (14)     0     | sc (25)      0    
    .  gut (1)    0     | icent (20)   0     | se (22)      0    
    .  hc1 (15)   0     | igut (2)     0     | ski (6)      0    
    .  hc2 (16)   0     | iliv1 (27)   0     | . ...        .    
    .  hc3 (17)   0     | iliv2 (28)   0     | . ...        .

Simulate out a quick example … we filter out just the dosing records and
fill in with a smooth set of observations

``` r
doses <- filter(data, EVID==1)

sims <- 
  mod %>% 
  mrgsim(data = doses, obsaug = TRUE) %>% 
  mutate(type = typef(ID))

ggplot(sims, aes(TIME, CP, col = type)) + 
  geom_line(lwd = 1) + 
  scale_x_continuous(breaks = seq(0, 12, 2)) + 
  scale_y_continuous(trans = "log10") + 
  ylab("Pitavastatin concentration") + xlab("Time (hour)") 
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/oatp_ddi_optimization-unnamed-chunk-7-1.png)<!-- -->

``` r
sims %>% 
  group_by(type) %>% 
  pmxTools::get_auc(dv = "CP") %>% 
  mutate(fold_increase = AUC/first(AUC))
```

    .   ID       AUC fold_increase
    . 1  1  44.07871      1.000000
    . 2  2 161.07202      3.654191

# Objective function

- Least squares objective function
- Weighted by the observations

``` r
wss <- function(dv, pred, par = NULL) {
  sum(((dv-pred)/dv)^2, na.rm = TRUE)
}
```

### Prediction function

- Let’s go through step by step what each line is doing for us

``` r
pred <- function(p, m, d, pred = FALSE) {
  
  p <- lapply(p, exp)
  
  names(p) <- names(theta)
  
  m <- param(m, p)
  
  if(pred) {
    out <- mrgsim_df(m, data = d, recover = "type")
    return(out)
  }
  
  out <- mrgsim_q(m, data = d, output = "df")
  
  return(wss(d$DV, out$CP))
  
  #-1*sum(dnorm(log(yobs),log(out$CP),.par$sigma,log=TRUE), na.rm=TRUE)
}
```

# Optimize with different methods

## `minqa::newuoa`: minimization without derivatives

- These appear to be the parameters that the authors are fitting

``` r
theta <- log(c(fbCLintall = 1.2, ikiu = 1.2, fbile = 0.9, ka = 0.1, ktr = 0.1))

control <- list(iprint = 25)

fit1 <- newuoa(theta, pred, m = mod, d = data, control = control)
```

    . npt = 7 , n =  5 
    . rhobeg =  0.460517 , rhoend =  4.60517e-07 
    . start par. =  0.1823216 0.1823216 -0.1053605 -2.302585 -2.302585 fn =  11.02654 
    . rho:    0.046 eval:   8 fn:      7.75648 par:-0.278195 0.182322 -0.105361 -2.30259 -2.30259 
    .  25:     5.9761213: -0.384231 0.0362186 0.00414868 -1.21552 -2.32950
    .  50:     5.3741109: -0.308723 0.648357 0.292933 -0.905966 -1.96759
    . rho:   0.0046 eval:  51 fn:      5.32508 par:-0.353296 0.655989 0.284855 -0.909043 -1.96654 
    .  75:     5.3059744: -0.350731 0.459433 0.217035 -0.850124 -1.92008
    . 100:     4.3785352: -0.370051 -2.44139 -0.294690 -0.577846 -1.60233
    . 125:     1.3293983: -0.247092 -4.51333 -0.870589 -0.331008 -1.15243
    . 150:    0.93251296: -0.202694 -4.66897 -1.12121 -0.210335 -0.727268
    . 175:    0.68900360: -0.201159 -4.52256 -1.07453 -0.0104370 -0.399597
    . rho:  0.00046 eval: 194 fn:     0.686165 par:-0.205710 -4.51109 -1.07050 -0.0167822 -0.374735 
    . 200:    0.68611383: -0.206148 -4.51161 -1.06892 -0.0115912 -0.374251
    . rho:  4.6e-05 eval: 220 fn:     0.686078 par:-0.204591 -4.51402 -1.06778 -0.0113134 -0.371530 
    . 225:    0.68607711: -0.204374 -4.51398 -1.06767 -0.0112089 -0.371226
    . rho:  4.6e-06 eval: 237 fn:     0.686077 par:-0.204382 -4.51408 -1.06765 -0.0111493 -0.371452 
    . rho:  4.6e-07 eval: 246 fn:     0.686077 par:-0.204382 -4.51408 -1.06765 -0.0111493 -0.371452 
    . 250:    0.68607681: -0.204382 -4.51408 -1.06765 -0.0111494 -0.371452
    . At return
    . eval: 255 fn:     0.68607672 par: -0.204382 -4.51408 -1.06765 -0.0111493 -0.371452

``` r
fit1$par <- setNames(fit1$par, names(theta))
```

### Get some predictions to look at how the fit went

- Predictions with the final estimates
- Predictions with the initial estimates
- Observed data to overlay

``` r
df_pred <- pred(fit1$par, mod, doses, pred = TRUE) %>% mutate(type = typef(type))
df_init <- pred(theta, mod, doses, pred = TRUE) %>% mutate(type = typef(type))
df_obs <- mutate(data, type = typef(type))
```

### Plot

``` r
ggplot(data = df_pred, lwd = 0.7) + 
  geom_line(data = df_init,aes(TIME, CP, lty = "A")) +
  geom_line(aes(TIME, CP, lty = "B")) + 
  geom_point(data = df_obs, aes(TIME, DV, col = type), size = 3) + 
  facet_wrap(~type) + 
  scale_y_log10(limits=c(0.1, 100)) +
  ylab("Pitavastatin concentration (ng/mL)") +
  scale_x_continuous(breaks = seq(0, 14, 2), name = "Time (hours)") +
  scale_linetype_manual(
    values = c(2,1), 
    labels = c("Initial estimates", "Final estimates"), 
    guide = "none",
    name = ""
  )
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/oatp_ddi_optimization-unnamed-chunk-12-1.png)<!-- -->

### The final objective function value and estimates

``` r
pred(fit1$par, m = mod, d = data)
```

    . [1] 0.6860767

``` r
exp(fit1$par)
```

    . fbCLintall       ikiu      fbile         ka        ktr 
    . 0.81515103 0.01095373 0.34381394 0.98891267 0.68973229

## `optim`: Nelder-Mead

This optimizer comes from the `stats` package in base R.

``` r
fit1b <- stats::optim(theta, pred, m = mod, d = data)
```

## `DEoptim`: differential evolution algorithm

“Performs evolutionary global optimization via the Differential
Evolution algorithm.”

``` r
lower <- rep(-6, length(theta)) %>% setNames(names(theta))
upper <- rep( 4, length(theta)) %>% setNames(names(theta))

set.seed(330303)

decontrol <- DEoptim.control(
  NP = 10*length(theta), 
  CR = 0.925, 
  F = 0.85,
  itermax = 90, 
  storepopfrom = 0, 
  trace = 10
)

fit2 <- DEoptim(
  fn = pred, 
  lower = lower, 
  upper = upper, 
  control = decontrol,
  m = mod, 
  d = data
)
```

    . Iteration: 10 bestvalit: 3.516849 bestmemit:   -0.632696   -3.722511   -1.360311   -0.036983   -0.612196
    . Iteration: 20 bestvalit: 2.350025 bestmemit:    0.000703   -4.042717   -0.719919    0.230071   -0.511443
    . Iteration: 30 bestvalit: 1.984471 bestmemit:    0.068825   -4.628240   -0.675047   -0.255022   -0.762239
    . Iteration: 40 bestvalit: 0.875897 bestmemit:   -0.125607   -4.415018   -0.963055   -0.076342   -0.074339
    . Iteration: 50 bestvalit: 0.702213 bestmemit:   -0.217592   -4.554969   -1.094951   -0.058136   -0.450317
    . Iteration: 60 bestvalit: 0.702213 bestmemit:   -0.217592   -4.554969   -1.094951   -0.058136   -0.450317
    . Iteration: 70 bestvalit: 0.687438 bestmemit:   -0.196411   -4.523280   -1.054508   -0.014670   -0.361633
    . Iteration: 80 bestvalit: 0.686394 bestmemit:   -0.207140   -4.517563   -1.073372   -0.011188   -0.382893
    . Iteration: 90 bestvalit: 0.686113 bestmemit:   -0.206033   -4.512065   -1.069557   -0.010903   -0.373121

### Plot the history

We can plot each population from the DE optimization over time as well
as the best member for each iteration.

``` r
hx <- lapply(fit2$member$storepop, as.data.frame) %>% bind_rows()
hx <- mutate(hx, iteration = rep(1:decontrol$itermax, each = decontrol$NP))
hx <- mutate(hx, pop = rep(1:decontrol$NP, time = decontrol$itermax))
hxm <- pivot_longer(hx, fbCLintall:ktr) %>% mutate(value = exp(value))
best <- 
  fit2$member$bestmemit %>% 
  as_tibble() %>% 
  mutate(iteration = seq(decontrol$itermax))
bestm <- pivot_longer(best, fbCLintall:ktr) %>% mutate(value = exp(value))
```

``` r
ggplot(data = hxm) + 
  geom_line(aes(iteration, value, group = pop), col = "darkslateblue") + 
  geom_line(data = bestm, aes(iteration ,value), col = "orange", lwd = 1) + 
  scale_y_log10(name = "Parameter value") + 
  facet_wrap(~name, ncol = 2, scales = "free_y")
```

![](/Users/kyleb/git/metrumresearchgroup/pbpk-qsp-mrgsolve/docs/img/oatp_ddi_optimization-unnamed-chunk-15-1.png)<!-- -->

## `GenSA`: simulated annealing

``` r
set.seed(11001)

sacontrol <- list(maxit = 100, nb.stop.improvement = 20, verbose = TRUE)

fit3 <- GenSA(NULL, pred, lower, upper, m = mod, d = data, control = sacontrol)
```

    . Initializing par with random data inside bounds
    . It: 1, obj value: 6.303744151
    . It: 24, obj value: 0.6860791496

## `hydroPSO`: particle swarm optimization

``` r
set.seed(2202201)

fit4 <- hydroPSO(
  theta, 
  fn = pred,
  lower = lower, 
  upper = upper, 
  control = list(maxit = 100, REPORT = 5),
  m = mod, d = data
)
```

# Compare optimization methods

``` r
results <- list(theta, fit1$par, fit1b$par, fit2$optim$bestmem, fit3$par, fit4$par)

results <- map(results, exp)

tibble(
  method = c("initial", "newuoa", "nelder", "RcppDE", "SA", "PSO"),
  fbCLintall = map_dbl(results, "fbCLintall"), 
  ikiu = map_dbl(results, "ikiu"), 
  fbile = map_dbl(results, "fbile"), 
  ka = map_dbl(results, "ka"), 
  ktr = map_dbl(results, "ktr")
) %>% kable(digits = 4)
```

| method  | fbCLintall |   ikiu |  fbile |     ka |    ktr |
|:--------|-----------:|-------:|-------:|-------:|-------:|
| initial |     1.2000 | 1.2000 | 0.9000 | 0.1000 | 0.1000 |
| newuoa  |     0.8152 | 0.0110 | 0.3438 | 0.9889 | 0.6897 |
| nelder  |     0.8387 | 0.0106 | 0.3572 | 0.9904 | 0.7114 |
| RcppDE  |     0.8138 | 0.0110 | 0.3432 | 0.9892 | 0.6886 |
| SA      |     0.8155 | 0.0110 | 0.3439 | 0.9885 | 0.6904 |
| PSO     |     0.8150 | 0.0110 | 0.3436 | 0.9883 | 0.6897 |
