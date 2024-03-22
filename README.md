# Calculation of the ROC-GLM and AUC confidence interval for distributed non-disclosive analysis in DataSHIELD

<!-- README.md is generated from README.Rmd. Please edit that file -->

- [About the repository](#about-the-repository)
- [Reproduce the results from the paper](#reproduce-the-results-from-the-paper)
- [Inspect results using the Docker](#inspect-results-using-the-docker)
- [Analyse results](#analyse-results) 
  - [Setup](#setup)
  - [Simulation Results](#simulation-results)
  - [Figures](#figures)
  - [Table of AUC values](#table-of-auc-values)
  - [Visualization of Gaussian mechanism](#visualization-of-gaussian-mechanism)
- [Session Info](#session-info)

## About the repository

This repository contains all code that is needed to run/reproduce the
simulation study for the respective paper. It also contains the code
used to create the figures and tables.

### Structure

- File to run the simulations: `simulation.R`
- `R` code: 
  - Helper function for calculating the Probit regression, AUC values, and confidence intervals: `R/helper.R`
  - Code to define and add the experiments using `batchtools` (called by `simulation.R`): `R/add-experiments.R`
  - Installation script to get all required package with the version used for the benchmark: `R/install-pkgs-versions.R`
  - Setup variables like repetitions or the grid of values for the simulation for the benchmark: `R/setup.R`
- Batchtools registry and results: `batchtools` (see [below](#simulation-results) how to load the results)
- Folder for the paper figures `*.pdf` and README figures `*.png`: `figures`
- Folder containing `*.tex` files of tables: `tables`

## Reproduce the results from the paper

To fully reproduce the results run the `simulation.R` script. To
reproduce the figures of the paper, render the README with
`rmarkdown::render("README.Rmd")`. When rendering the README, all
figures are created and stored in `figures` while the table is stored in
`tables`.

## Inspect results using the Docker

Running the
[docker](https://hub.docker.com/r/schalkdaniel/simulations-distr-auc)
provides an RStudio API in your browser with all packages pre-installed
and data to inspect the results. Therefore do:

1. Get the docker:

- **Build the docker manually:** Run `sudo docker build -t schalkdaniel/simulations-distr-auc .` (You can use whatever tag you like, but for consistency we use `schalkdaniel/simulations-distr-auc`)
- **Pull the docker:** Run `sudo docker pull schalkdaniel/simulations-distr-auc`

1. Run the docker: `sudo docker run -d -p 8787:8787 -e PASSWORD=test schalkdaniel/simulations-distr-auc`
2. Open your browser and visit `localhost:8787`
3. Login with `rstudio` as user and `test` as password

## Analyse results

The following code was used to load and analyse the results. It can be
used to fully reproduce the figures and tables of the paper.

### Setup

```r
# source(here::here("R/install-pkgs-versions.R"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(ggridges)
library(knitr)
library(batchtools)


### Theme setup:
## ==================================================================

# Set font (if available):
font = "Tinos"
sysfonts::font_add_google(font, font)
extrafont::loadfonts()
ft = extrafont::fonttable()

#if (all(! grepl(font, x = ft$FamilyName))) {
if (TRUE) {
  my_theme = theme_minimal()
} else {
  my_theme = theme_minimal(base_family = font)
}

theme_set(
  my_theme +
  theme(
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 9),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", face = "bold", size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),

    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
)

# to determine width, use the latex package `layouts` and put
# `\printinunitsof{mm}\prntlen{\textwidth}` into the tex file.
textwidth = 148.92

## Helper functions:
## =====================================================================

source(here::here("R/helper.R"))

#' Open pdf files using evince on linux machines
#' @param file to the pdf
evince = function(file) system(paste0("evince ", file, " &"))
```

### Simulation Results

```r
## Load data
loadRegistry(here::here("batchtools/"))
#> Experiment Registry
#>   Backend   : Interactive
#>   File dir  : /nfsmb/koll/rrehms/auc_paper/paper_commits/simulations-distr-auc/batchtools
#>   Work dir  : /nfsmb/koll/rrehms/auc_paper/paper_commits/simulations-distr-auc
#>   Jobs      : 126
#>   Problems  : 1
#>   Algorithms: 1
#>   Seed      : 31415
#>   Writeable : FALSE

### Load ROC-GLM approximation::
results = reduceResultsList()

# Split results:
# - 1. element of list contains the ordinary ROC-GLM.
# - All other elements contain the distributed ROC-GLM with differential privacy.
aucs_emp = do.call(rbind, lapply(results[[1]], function(ll) ll$aucs))
aucs_dp0 = do.call(rbind, lapply(results[-1], function(ll) {
  do.call(rbind, lapply(ll, function(l) l$aucs))
}))

## Manual inspection:
#aucs_dp0 %>%
#  filter(l2sens == 0.07, epsilon == 0.5, delta == 0.5)


### Distributed ROC-GLM approximations:
#lbreaks = c((0.8 - 0.6) * 0.01, 1 * 0.01)
lbreaks = c(0.01, 0.05)
app_acc = rbind(
  data.frame(lower = -Inf,       upper = lbreaks[1], qual = paste0("$\\Delta AUC < ", lbreaks[1], "$")),
  #data.frame(lower = lbreaks[1], upper = lbreaks[2], qual = paste0("$\\Delta AUC < ", lbreaks[2], "$")),
  data.frame(lower = lbreaks[2], upper = Inf,        qual = "unacceptable")
)

aucs_dp = aucs_dp0 %>%
  mutate(
    auc_diff = abs(auc_emp - auc_roc),
    st_auc_diff = abs(auc_emp - auc_roc) / auc_emp,
    qual = cut(auc_diff, breaks = c(app_acc$lower, Inf), labels = app_acc$qual))
```

### Figures

#### AUC densities, Figure 3, Section 5.3.1

```r
## Comparison for different n intervals:
#nb      = c(100, 200, 400, 800, 1600, 2500)
nb      = seq(100, 2500, length.out = 7)
clabels = paste0("n in (", nb[-length(nb)], ", ", nb[-1], "]")
ncut    = cut(aucs_emp$n, breaks = nb, labels = clabels)
clabels = paste0(clabels, " (Count: ", table(ncut), ")")
ncut    = cut(aucs_emp$n, breaks = nb, labels = clabels)

aucs_emp_plt = aucs_emp %>%
  mutate(n_cat = ncut) %>%
  select(auc_emp, auc_roc, n_cat) %>%
  na.omit()

gg_den_both_facet = ggplot(data = aucs_emp_plt) +
  geom_histogram(aes(x = auc_emp, color = "Empirical AUC", fill = "Empirical AUC", y = ..density..),
     color = "white", size = 0.2, alpha = 0.4) +
  geom_histogram(aes(x = auc_roc, color = "AUC (ROC-GLM)", fill = "AUC (ROC-GLM)", y = ..density..),
     color = "white", size = 0.2, alpha = 0.2) +
  geom_histogram(aes(x = auc_emp, color = "Empirical AUC", fill = "Empirical AUC", y = ..density..),
     color = "white", size = 0.2, alpha = 0) +
  geom_density(aes(x = auc_emp, color = "Empirical AUC", fill = "Empirical AUC"),
     fill = "transparent", size = 0.6, alpha = 0.4) +
  geom_density(aes(x = auc_roc, color = "AUC (ROC-GLM)", fill = "AUC (ROC-GLM)"),
     fill = "transparent", size = 0.6, alpha = 0.2) +
  theme(legend.position = "bottom", legend.key.size = unit(0.2, "cm")) +
  xlab("AUC") +
  ylab("Density") +
  labs(fill = "", color = "") +
  scale_color_uchicago() +
  scale_fill_uchicago() +
  facet_wrap(~ n_cat)

gg_den_both_facet
```

![](figures/unnamed-chunk-4-1.png)<!-- -->

```r
ggsave(plot = gg_den_both_facet,
  filename = here::here("figures/auc-emp-density-facets.pdf"),
  width = textwidth,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/auc-emp-density-facets.pdf"))
```

#### Discrepance as Area between empirical AUC and ROC GLM

```r
gg_l1_dist = ggplot(data = aucs_emp,aes(x = as.numeric(discr), y = ..density..))+
  geom_histogram(color = "white", size = 0.2, alpha = 0.9, fill = "grey60")+
    geom_density(color = "grey40")
gg_l1_dist
```

![](figures/unnamed-chunk-5-1.png)<!-- -->

```r
ggsave(plot = gg_l1_dist,
  filename = here::here("figures/l1_disc.pdf"),
  width = textwidth * 1,
  height = textwidth * 0.5,
  units = "mm")
```

### Approximation errors

#### Figure 4, Section 5.3.1

```r
## AUC VALUES
## =======================================================

l2senss  = sort(unique(aucs_dp$l2sens))
auc_cuts = seq(0.5, 1, length.out = 21L)
ggs_auc  = list()
for (i in seq_along(l2senss)) {
  df_dp = aucs_dp %>% filter(l2sens == l2senss[i])

  tab = df_dp %>%
    mutate(
      auc_emp_cut = cut(x = auc_emp, breaks = auc_cuts)
    ) %>%
    group_by(auc_emp_cut, epsilon, delta) %>%
    summarize(
      Min.      = min(auc_diff, na.rm = TRUE),
      "1st Qu." = quantile(auc_diff, 0.25, na.rm = TRUE),
      Median    = median(auc_diff, na.rm = TRUE),
      Mean      = mean(auc_diff, na.rm = TRUE),
      "3rd Qu." = quantile(auc_diff, 0.75, na.rm = TRUE),
      Max.      = max(auc_diff, na.rm = TRUE),
      Sd.       = sd(auc_diff, na.rm = TRUE),
      Count     = n(),
      MAE       = mean(auc_diff, na.rm = TRUE),
      RMSE      = sqrt(mean(auc_diff^2)),
      stMAE     = mean(st_auc_diff, na.rm = TRUE)) %>%
    mutate(qual = cut(abs(MAE), breaks = c(app_acc$lower, Inf), labels = app_acc$qual)) %>%
    select(qual, auc_emp_cut, delta, epsilon) %>%
    na.omit()

  tab$qual = factor(tab$qual, levels = rev(app_acc$qual))

  ggs_auc[[i]] = ggplot(tab, aes(x = "1", y = auc_emp_cut, color = qual, fill = qual)) +
    geom_tile() +
    scale_y_discrete(limits = rev, breaks = levels(tab$auc_emp_cut)[c(5, 11, 17)],
      labels = auc_cuts[c(5, 11, 17)]) +
    xlab("") +
    theme(axis.text.x = element_blank()) +
    scale_color_npg(labels = unname(latex2exp::TeX(rev(app_acc$qual)))) +
    scale_fill_npg(labels = unname(latex2exp::TeX(rev(app_acc$qual)))) +
    labs(color = "Accuracy", fill = "Accuracy") +
    facet_grid(delta ~ epsilon, labeller = label_bquote(delta == .(delta), epsilon == .(epsilon))) +
    ggtitle(latex2exp::TeX(paste0("$\\Delta_2(f) = ", l2senss[i], "$")))

  if (i == 1) {
    ggs_auc[[i]] = ggs_auc[[i]] +
      ylab("AUC bin") +
      theme(
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        strip.text = element_text(color = "black", face = "bold", size = 5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 270)
      )
  } else {
    ggs_auc[[i]] = ggs_auc[[i]] +
      ylab("") +
      theme(
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        strip.text = element_text(color = "black", face = "bold", size = 5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"),
        axis.text.y = element_blank(),
        strip.text.x = element_text(angle = 270),
        axis.title.y = element_blank()
      )
  }
}

gg_auc_dp = do.call(ggpubr::ggarrange, c(ggs_auc, list(nrow = 1,
  common.legend = TRUE, legend = "bottom",
  widths = c(1, rep(0.75, length(ggs_auc) - 1)))))

gg_auc_dp
```

![](figures/unnamed-chunk-6-1.png)<!-- -->

```r
ggsave(plot = gg_auc_dp,
  filename = here::here("figures/auc-diff-priv.pdf"),
  width = textwidth * 1,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/auc-diff-priv.pdf"))
```

#### Figure 5, Section 5.3.2

```r
## CI BOUNDARIES
## =======================================================

ggs_ci = list()
for (i in seq_along(l2senss)) {

  df_dp = aucs_dp %>%
    filter(l2sens == l2senss[i])

  tab = df_dp %>%
    mutate(
      auc_app_cut = cut(x = auc_roc, breaks = auc_cuts),
      auc_emp_cut = cut(x = auc_emp, breaks = auc_cuts)
    ) %>%
    group_by(auc_emp_cut, epsilon, delta) %>%
    summarize(
      Min.      = min(delta_ci, na.rm = TRUE),
      "1st Qu." = quantile(delta_ci, 0.25, na.rm = TRUE),
      Median    = median(delta_ci, na.rm = TRUE),
      Mean      = mean(delta_ci, na.rm = TRUE),
      "3rd Qu." = quantile(delta_ci, 0.75, na.rm = TRUE),
      Max.      = max(delta_ci, na.rm = TRUE),
      Sd.       = sd(delta_ci, na.rm = TRUE),
      Count     = n()) %>%
    mutate(qual_ci = cut(Mean, breaks = c(-Inf, 0.01, Inf),
      labels = c("$\\Delta ci_{\\alpha} < 0.01", "unacceptable")))

  tab$qual_ci = factor(tab$qual_ci, levels = rev(levels(tab$qual_ci)))

  ggs_ci[[i]] = ggplot(na.omit(tab), aes(x = "1", y = auc_emp_cut, color = qual_ci, fill = qual_ci)) +
    geom_tile() +
    scale_y_discrete(limits = rev, breaks = levels(tab$auc_emp_cut)[c(5, 11, 17)],
      labels = auc_cuts[c(5, 11, 17)]) +
    xlab("") +
    theme(axis.text.x = element_blank()) +
    scale_color_npg(labels = unname(latex2exp::TeX(levels(tab$qual_ci)))) +
    scale_fill_npg(labels = unname(latex2exp::TeX(levels(tab$qual_ci)))) +
    labs(fill = "Accuracy") +
    guides(color = "none") +
    facet_grid(delta ~ epsilon, labeller = label_bquote(delta == .(delta), epsilon == .(epsilon))) +
    ggtitle(latex2exp::TeX(paste0("$\\Delta_2(f) = ", l2senss[i], "$")))

  if (i == 1) {
    ggs_ci[[i]] = ggs_ci[[i]] +
      ylab("AUC bin") +
      theme(
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        strip.text = element_text(color = "black", face = "bold", size = 5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 270)
      )
  } else {
    ggs_ci[[i]] = ggs_ci[[i]] +
      ylab("") +
      theme(
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        strip.text = element_text(color = "black", face = "bold", size = 5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"),
        axis.text.y = element_blank(),
        strip.text.x = element_text(angle = 270),
        axis.title.y = element_blank()
      )
  }
}

gg_ci_dp = do.call(ggpubr::ggarrange, c(ggs_ci, list(nrow = 1,
  common.legend = TRUE, legend = "bottom",
  widths = c(1, rep(0.75, length(ggs_ci) - 1)))))

gg_ci_dp
```

![](figures/unnamed-chunk-7-1.png)<!-- -->

```r
ggsave(plot = gg_ci_dp,
  filename = here::here("figures/cis-diff-priv.pdf"),
  width = textwidth,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/cis-diff-priv.pdf"))
```

### Table of AUC values

#### Table 1, Section 5.3.1

```r
# Discretize empirical AUC values into 20 bins between 0.5 and 1.
tab = aucs_emp %>% mutate(
    auc_diff = auc_emp - auc_roc,
    auc_app_cut = cut(x = auc_roc, breaks = seq(0.5, 1, length.out = 21L)),
    auc_emp_cut = cut(x = auc_emp, breaks = seq(0.5, 1, length.out = 21L))
  ) %>%
  group_by(auc_emp_cut) %>%
  # Calculate summary statistics per bin:
  summarize(
    "Min."    = min(auc_diff, na.rm = TRUE),
    "1st Qu." = quantile(auc_diff, 0.25, na.rm = TRUE),
    "Median"  = median(auc_diff, na.rm = TRUE),
    "Mean"    = mean(auc_diff, na.rm = TRUE),
    "3rd Qu." = quantile(auc_diff, 0.75, na.rm = TRUE),
    "Max."    = max(auc_diff, na.rm = TRUE), Sd. = sd(auc_diff, na.rm = TRUE), Count = n()) %>%
  na.omit()

tab0 = tab

tab$auc_emp_cut = paste0("$", tab$auc_emp_cut, "$")
tab$Count = paste0("$", tab$Count, "$")
tnames = names(tab)
for (tn in tnames[-c(1, length(tnames))]) {
  tab[[tn]] = ifelse((abs(tab[[tn]]) > 0.01) & (tn %in% c("Median", "Mean")),
    paste0("$\\mathbf{", sprintf("%.4f", round(tab[[tn]], 4)), "}$"),
    paste0("$", sprintf("%.4f", round(tab[[tn]], 4)), "$"))
}

tab_latex = tab %>% kable(format = "latex", escape = FALSE)
writeLines(tab_latex, here::here("tables/auc-approximations.tex"))

tab0 %>% kable()
```

| auc_emp_cut  |       Min. |    1st Qu. |     Median |       Mean |    3rd Qu. |      Max. |       Sd. | Count |
|:-------------|-----------:|-----------:|-----------:|-----------:|-----------:|----------:|----------:|------:|
| (0.5,0.525\]  | \-0.0044383 | \-0.0001307 |  0\.0004448 |  0\.0038600 |  0\.0013542 | 0\.0495491 | 0\.0098073 |   431 |
| (0.525,0.55\] | \-0.0052108 |  0\.0000685 |  0\.0006418 |  0\.0027560 |  0\.0011252 | 0\.0985955 | 0\.0122949 |   505 |
| (0.55,0.575\] | \-0.0031457 |  0\.0003188 |  0\.0008907 |  0\.0014210 |  0\.0014852 | 0\.1298428 | 0\.0080411 |   465 |
| (0.575,0.6\]  | \-0.0018138 |  0\.0005808 |  0\.0011924 |  0\.0015039 |  0\.0017489 | 0\.1567206 | 0\.0071699 |   482 |
| (0.6,0.625\]  | \-0.0043680 |  0\.0008882 |  0\.0015116 |  0\.0014311 |  0\.0020180 | 0\.0063648 | 0\.0010403 |   485 |
| (0.625,0.65\] | \-0.0038982 |  0\.0011677 |  0\.0017045 |  0\.0017158 |  0\.0022246 | 0\.0069136 | 0\.0010315 |   501 |
| (0.65,0.675\] | \-0.0031417 |  0\.0012548 |  0\.0017989 |  0\.0017817 |  0\.0023333 | 0\.0068008 | 0\.0010851 |   503 |
| (0.675,0.7\]  | \-0.0022132 |  0\.0012076 |  0\.0017841 |  0\.0017879 |  0\.0023216 | 0\.0063895 | 0\.0010377 |   465 |
| (0.7,0.725\]  | \-0.0081956 |  0\.0010002 |  0\.0016267 |  0\.0016282 |  0\.0022721 | 0\.0070282 | 0\.0011775 |   523 |
| (0.725,0.75\] | \-0.0030765 |  0\.0008437 |  0\.0014636 |  0\.0014205 |  0\.0021148 | 0\.0086789 | 0\.0012128 |   485 |
| (0.75,0.775\] | \-0.0058233 |  0\.0004062 |  0\.0010563 |  0\.0010395 |  0\.0018143 | 0\.0053310 | 0\.0012823 |   501 |
| (0.775,0.8\]  | \-0.0052752 | \-0.0003056 |  0\.0003994 |  0\.0004837 |  0\.0012293 | 0\.0087512 | 0\.0014813 |   523 |
| (0.8,0.825\]  | \-0.0060755 | \-0.0012676 | \-0.0002452 | \-0.0003585 |  0\.0004864 | 0\.0045295 | 0\.0015518 |   476 |
| (0.825,0.85\] | \-0.0124976 | \-0.0022868 | \-0.0013400 | \-0.0014091 | \-0.0002653 | 0\.0059060 | 0\.0019481 |   484 |
| (0.85,0.875\] | \-0.0110969 | \-0.0036791 | \-0.0026149 | \-0.0025210 | \-0.0013951 | 0\.0074169 | 0\.0020318 |   520 |
| (0.875,0.9\]  | \-0.0136367 | \-0.0056322 | \-0.0044412 | \-0.0043459 | \-0.0030480 | 0\.0076091 | 0\.0022857 |   534 |
| (0.9,0.925\]  | \-0.0195223 | \-0.0079785 | \-0.0064726 | \-0.0065343 | \-0.0052044 | 0\.0065504 | 0\.0026299 |   515 |
| (0.925,0.95\] | \-0.0192918 | \-0.0105496 | \-0.0090891 | \-0.0089337 | \-0.0076055 | 0\.0056324 | 0\.0029665 |   481 |
| (0.95,0.975\] | \-0.0226558 | \-0.0137940 | \-0.0112812 | \-0.0113273 | \-0.0092846 | 0\.0066963 | 0\.0037026 |   503 |
| (0.975,1\]    | \-0.0180451 | \-0.0093118 | \-0.0061871 | \-0.0064037 | \-0.0034321 | 0\.0012870 | 0\.0039401 |   529 |

### Visualization of Gaussian mechanism

```r
set.seed(31415)
x = runif(10, 0, 1)

l2s = c(0.01, 0.05, 0.1)
epsilons = seq(0.1, 0.5, length.out = 3L)
deltas = seq(0.1, 0.5, length.out = 3L)

getSD = function(l2s, epsilon, delta) sqrt(2 * log(1.25 / delta)) * l2s / epsilon

pars = expand.grid(x = x, l2s = l2s, epsilon = epsilons, delta = deltas)
pars$sd = getSD(pars$l2s, pars$epsilon, pars$delta)
pars$xnew = rnorm(n = nrow(pars), pars$x, sd = pars$sd)
pars$xnew = ifelse(pars$xnew > 1, 1, pars$xnew)
pars$xnew = ifelse(pars$xnew < 0, 0, pars$xnew)
pars$id = rep(seq_along(x), times = nrow(pars) / length(x))

ll = list()
for (i in seq_len(nrow(pars))) {
  sdd = pars$sd[i]
  xx = seq(pars$x[i] - 4 * sdd, pars$x[i] + 4 * sdd, length.out = 100L)
  yx = dnorm(x = xx, mean = pars$x[i], sd = sdd)
  ll[[i]] = data.frame(x = pars$x[i], xx = xx, yx = yx, l2s = pars$l2s[i],
    epsilon = pars$epsilon[i], delta = pars$delta[i], group = as.character(i),
    xnew = pars$xnew[i], id = pars$id[i], sdd = sdd)
}
df_plt = do.call(rbind, ll)

library(ggplot2)
library(dplyr)
library(ggsci)

df_text = df_plt %>%
  group_by(l2s, delta) %>%
  mutate(ymax = max(yx)) %>%
  group_by(l2s, epsilon) %>%
  mutate(xleft = min(xx)) %>%
  group_by(group, id, l2s) %>%
  filter(row_number() == 1) %>%
  mutate(y = 0.1 * ymax) %>%
  group_by(l2s, epsilon, delta) %>%
  mutate(order_original = order(order(x)), order_new = order(order(xnew)))
  #mutate(y = (ymax - yx) * 0.05)

order(df_text[1:10, ]$x)
#>  [1]  5  3  7  2  6  8  4 10  9  1

df_tau = df_text %>%
  group_by(l2s, epsilon, delta) %>%
  filter(row_number() == 1) %>%
  mutate(sd_label = latex2exp::TeX(paste0("$\\tau = ",  round(sdd, 3), "$"), output = "character"))

ggs_gm = list()
for (l2s0 in l2s) {
  gg_dp = ggplot() +
    geom_line(data = df_plt %>% filter(l2s == l2s0),
      mapping = aes(x = xx, y = yx, group = group, color = as.factor(id)), alpha = 0.6, show.legend = FALSE) +
    geom_segment(data = df_text %>% filter(l2s == l2s0), aes(x = x, y = 0, xend = xnew, yend = -2 * y)) +
    geom_label(data = df_text %>% filter(l2s == l2s0), aes(x = x, y = 0, label = order_original, fill = as.factor(id)),
      size = 1.5, show.legend = FALSE, color = "white", #family = font,
      fontface = "bold") +
    geom_label(data = df_text %>% filter(l2s == l2s0), aes(x = xnew, y = -2 * y, label = order_new, fill = as.factor(id)),
      size = 1.5, show.legend = FALSE, color = "white", #family = font,
      fontface = "bold") +
    geom_label(data = df_tau %>% filter(l2s == l2s0, epsilon == 0.1), aes(x = xleft, y = 0, label = "f(x)"),
      #family = font,
               size = 2, hjust = 0, fill = "white", label.size = 0) +
    geom_label(data = df_tau %>% filter(l2s == l2s0, epsilon == 0.1), aes(x = xleft, y = -2 * y, label = "f(x) + r"),
      #family = font,
               size = 2, hjust = -0, fill = "white", label.size = 0) +
    geom_label(data = df_tau %>% filter(l2s == l2s0), aes(x = xleft, y = y * 9, label = sd_label),
      parse = TRUE, fill = "white", hjust = -0, size = 3, label.size = 0#, family = font
      ) +
    facet_grid(delta ~ epsilon, scales = "free", labeller = label_bquote(delta == .(delta), epsilon == .(epsilon))) +
    scale_color_npg() +
    scale_fill_npg() +
    xlab(latex2exp::TeX("Score values $f(x)$")) +
    ylab(latex2exp::TeX("Density of $f(x) + r$ with $r \\sim N(0, \\tau^2)$ distribution")) +
    ggtitle(latex2exp::TeX(paste0("Gaussian Mechanism for $\\Delta_2(f) = ", l2s0, "$")))

  ggs_gm[[paste0("l2s", l2s0)]] = gg_dp

  ggsave(plot = gg_dp,
    filename = here::here(paste0("figures/gaussian-mechanism", l2s0, ".pdf")),
    width = textwidth,
    height = textwidth,
    units = "mm")
}
```

#### Appendix A.2, Figure 9

```r
ggs_gm[[1]]
```

![](figures/unnamed-chunk-10-1.png)<!-- -->

#### Appendix A.2, Figure 10

```r
ggs_gm[[2]]
```

![](figures/unnamed-chunk-11-1.png)<!-- -->

## Session Info

```r
sessionInfo()
#> R version 4.2.2 Patched (2022-11-10 r83330)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Debian GNU/Linux 12 (bookworm)
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so
#> 
#> locale:
#>  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] pROC_1.18.0       checkmate_2.0.0   batchtools_0.9.15 knitr_1.36       
#>  [5] ggridges_0.5.3    ggsci_2.9         gridExtra_2.3     ggplot2_3.3.5    
#>  [9] tidyr_1.1.4       dplyr_1.0.7      
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_1.7.3    carData_3.0-5     here_0.1          assertthat_0.2.1 
#>  [5] highr_0.10        base64url_1.4     yaml_2.3.7        progress_1.2.2   
#>  [9] Rttf2pt1_1.3.12   pillar_1.9.0      backports_1.4.1   glue_1.6.2       
#> [13] extrafontdb_1.0   digest_0.6.31     ggsignif_0.6.4    colorspace_2.1-0 
#> [17] cowplot_1.1.1     htmltools_0.5.4   plyr_1.8.8        pkgconfig_2.0.3  
#> [21] broom_1.0.3       sysfonts_0.8.1    purrr_1.0.1       scales_1.2.1     
#> [25] brew_1.0-8        tibble_3.2.1      generics_0.1.3    farver_2.1.1     
#> [29] car_3.1-1         ellipsis_0.3.2    ggpubr_0.3.0      withr_2.5.0      
#> [33] cli_3.6.0         magrittr_2.0.3    crayon_1.5.2      evaluate_0.20    
#> [37] fs_1.6.1          fansi_1.0.4       rstatix_0.7.2     textshaping_0.3.6
#> [41] tools_4.2.2       data.table_1.14.8 prettyunits_1.1.1 hms_1.1.2        
#> [45] lifecycle_1.0.3   stringr_1.5.0     munsell_0.5.0     compiler_4.2.2   
#> [49] systemfonts_1.0.4 rlang_1.1.2       grid_4.2.2        rappdirs_0.3.3   
#> [53] labeling_0.4.2    rmarkdown_2.11    gtable_0.3.1      abind_1.4-5      
#> [57] DBI_1.1.3         curl_5.0.0        R6_2.5.1          fastmap_1.1.1    
#> [61] extrafont_0.17    utf8_1.2.3        rprojroot_2.0.3   latex2exp_0.5.0  
#> [65] ragg_1.2.5        stringi_1.7.12    Rcpp_1.0.10       vctrs_0.6.4      
#> [69] tidyselect_1.2.0  xfun_0.41
```