Calculation of the ROC-GLM and AUC confidence interval for distributed
non-disclosive analysis in DataSHIELD
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

  - [About the repository](#about-the-repository)
  - [Reproduce the results from the
    paper](#reproduce-the-results-from-the-paper)
  - [Inspect results using the
    Docker](#inspect-results-using-the-docker)
  - [Analyse results](#analyse-results)
      - [Setup](#setup)
      - [Simulation Results](#simulation-results)
      - [Figures](#figures)
      - [Table of AUC values](#table-of-auc-values)
      - [Visualization of Gaussian
        mechanism](#visualization-of-gaussian-mechanism)
  - [Session Info](#session-info)

## About the repository

This repository contains all code that is needed to run/reproduce the
simulation study for the respective paper. It also contains the code
used to create the figures and tables.

### Structure

  - File to run the simulations: `simulation.R`
  - `R` code:
      - Helper function for calculating the Probit regression, AUC
        values, and confidence intervals: `R/helper.R`
      - Code to define and add the experiments using `batchtools`
        (called by `simulation.R`): `R/add-experiments.R`
      - Installation script to get all required package with the version
        used for the benchmark: `R/install-pkgs-versions.R`
      - Setup variables like repetitions or the grid of values for the
        simulation for the benchmark: `R/setup.R`
  - Batchtools registry and results: `batchtools` (see
    [below](#simulation-results) how to load the results)
  - Folder for the paper figures `*.pdf` and README figures `*.png`:
    `figures`
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

1.  Get the docker:

<!-- end list -->

  - **Build the docker manually:** Run `sudo docker build -t
    schalkdaniel/simulations-distr-auc .` (You can use whatever tag you
    like, but for consistency we use
    `schalkdaniel/simulations-distr-auc`)
  - **Pull the docker:** Run `sudo docker pull
    schalkdaniel/simulations-distr-auc`

<!-- end list -->

2.  Run the docker: `sudo docker run -d -p 8787:8787 -e PASSWORD=test
    schalkdaniel/simulations-distr-auc`
3.  Open your browser and visit <localhost:8787>
4.  Login with `rstudio` as user and `test` as password

## Analyse results

The following code was used to load and analyse the results. It can be
used to fully reproduce the figures and tables of the paper.

### Setup

``` r
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

``` r
## Load data
loadRegistry(here::here("batchtools/"))
#> Experiment Registry
#>   Backend   : Interactive
#>   File dir  : /home/daniel/repos/simulations-distr-auc/batchtools
#>   Work dir  : /home/daniel/repos/simulations-distr-auc
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

``` r
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

``` r

ggsave(plot = gg_den_both_facet,
  filename = here::here("figures/auc-emp-density-facets.pdf"),
  width = textwidth,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/auc-emp-density-facets.pdf"))
```

### Approximation errors

#### Figure 4, Section 5.3.1

``` r
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

![](figures/unnamed-chunk-5-1.png)<!-- -->

``` r

ggsave(plot = gg_auc_dp,
  filename = here::here("figures/auc-diff-priv.pdf"),
  width = textwidth * 1,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/auc-diff-priv.pdf"))
```

#### Figure 5, Section 5.3.2

``` r
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

![](figures/unnamed-chunk-6-1.png)<!-- -->

``` r

ggsave(plot = gg_ci_dp,
  filename = here::here("figures/cis-diff-priv.pdf"),
  width = textwidth,
  height = textwidth * 0.5,
  units = "mm")

#evince(here::here("figures/cis-diff-priv.pdf"))
```

### Table of AUC values

#### Table 1, Section 5.3.1

``` r
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

| auc\_emp\_cut |     Min. |  1st Qu. |   Median |     Mean |  3rd Qu. |   Max. |    Sd. | Count |
| :------------ | -------: | -------: | -------: | -------: | -------: | -----: | -----: | ----: |
| (0.5,0.525\]  | \-0.0044 | \-0.0001 |   0.0005 |   0.0040 |   0.0014 | 0.0506 | 0.0100 |   431 |
| (0.525,0.55\] | \-0.0052 |   0.0001 |   0.0006 |   0.0027 |   0.0011 | 0.0986 | 0.0123 |   505 |
| (0.55,0.575\] | \-0.0031 |   0.0003 |   0.0009 |   0.0014 |   0.0015 | 0.1298 | 0.0080 |   465 |
| (0.575,0.6\]  | \-0.0018 |   0.0006 |   0.0012 |   0.0015 |   0.0017 | 0.1567 | 0.0072 |   482 |
| (0.6,0.625\]  | \-0.0044 |   0.0009 |   0.0015 |   0.0014 |   0.0020 | 0.0064 | 0.0010 |   485 |
| (0.625,0.65\] | \-0.0039 |   0.0012 |   0.0017 |   0.0017 |   0.0022 | 0.0069 | 0.0010 |   501 |
| (0.65,0.675\] | \-0.0031 |   0.0013 |   0.0018 |   0.0018 |   0.0023 | 0.0068 | 0.0011 |   503 |
| (0.675,0.7\]  | \-0.0022 |   0.0012 |   0.0018 |   0.0018 |   0.0023 | 0.0064 | 0.0010 |   465 |
| (0.7,0.725\]  | \-0.0082 |   0.0010 |   0.0016 |   0.0016 |   0.0023 | 0.0070 | 0.0012 |   523 |
| (0.725,0.75\] | \-0.0031 |   0.0008 |   0.0015 |   0.0014 |   0.0021 | 0.0087 | 0.0012 |   485 |
| (0.75,0.775\] | \-0.0058 |   0.0004 |   0.0011 |   0.0010 |   0.0018 | 0.0053 | 0.0013 |   501 |
| (0.775,0.8\]  | \-0.0053 | \-0.0003 |   0.0004 |   0.0005 |   0.0012 | 0.0088 | 0.0015 |   523 |
| (0.8,0.825\]  | \-0.0061 | \-0.0013 | \-0.0002 | \-0.0004 |   0.0005 | 0.0045 | 0.0016 |   476 |
| (0.825,0.85\] | \-0.0125 | \-0.0023 | \-0.0013 | \-0.0014 | \-0.0003 | 0.0059 | 0.0019 |   484 |
| (0.85,0.875\] | \-0.0111 | \-0.0037 | \-0.0026 | \-0.0025 | \-0.0014 | 0.0074 | 0.0020 |   520 |
| (0.875,0.9\]  | \-0.0136 | \-0.0056 | \-0.0044 | \-0.0043 | \-0.0030 | 0.0076 | 0.0023 |   534 |
| (0.9,0.925\]  | \-0.0195 | \-0.0080 | \-0.0065 | \-0.0065 | \-0.0052 | 0.0066 | 0.0026 |   515 |
| (0.925,0.95\] | \-0.0193 | \-0.0105 | \-0.0091 | \-0.0089 | \-0.0076 | 0.0056 | 0.0030 |   481 |
| (0.95,0.975\] | \-0.0227 | \-0.0138 | \-0.0113 | \-0.0113 | \-0.0093 | 0.0067 | 0.0037 |   503 |
| (0.975,1\]    | \-0.0180 | \-0.0093 | \-0.0062 | \-0.0064 | \-0.0034 | 0.0013 | 0.0039 |   529 |

### Visualization of Gaussian mechanism

``` r
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

``` r
ggs_gm[[1]]
```

![](figures/unnamed-chunk-9-1.png)<!-- -->

#### Appendix A.2, Figure 10

``` r
ggs_gm[[2]]
```

![](figures/unnamed-chunk-10-1.png)<!-- -->

## Session Info

``` r
sessionInfo()
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Arch Linux
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/libblas.so.3.10.0
#> LAPACK: /usr/lib/liblapack.so.3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_GB.UTF-8   
#>  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] pROC_1.18.0       checkmate_2.0.0   batchtools_0.9.15 knitr_1.36        ggridges_0.5.3    ggsci_2.9         gridExtra_2.3    
#>  [8] ggplot2_3.3.5     tidyr_1.2.0       dplyr_1.0.8      
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_1.8.0    carData_3.0-4     here_0.1          assertthat_0.2.1  highr_0.8         prettycode_1.1.0  base64url_1.4    
#>  [8] cellranger_1.1.0  yaml_2.2.1        progress_1.2.2    Rttf2pt1_1.3.8    pillar_1.7.0      backports_1.4.1   glue_1.6.2       
#> [15] extrafontdb_1.0   digest_0.6.29     ggsignif_0.6.0    colorspace_2.0-2  cowplot_1.0.0     htmltools_0.4.0   plyr_1.8.6       
#> [22] pkgconfig_2.0.3   broom_0.7.1       haven_2.4.3       sysfonts_0.8.1    purrr_0.3.4       scales_1.1.1      brew_1.0-6       
#> [29] openxlsx_4.2.2    rio_0.5.16        tibble_3.1.6      generics_0.1.2    farver_2.1.0      car_3.0-10        ellipsis_0.3.2   
#> [36] ggpubr_0.3.0      withr_2.4.3       cli_3.2.0         readxl_1.3.1      magrittr_2.0.2    crayon_1.5.0      evaluate_0.14    
#> [43] fs_1.5.0          fansi_1.0.2       forcats_0.5.1     rstatix_0.5.0     foreign_0.8-81    textshaping_0.3.6 tools_4.1.2      
#> [50] data.table_1.14.2 prettyunits_1.1.1 hms_1.1.1         lifecycle_1.0.1   stringr_1.4.0     munsell_0.5.0     zip_2.1.1        
#> [57] compiler_4.1.2    systemfonts_1.0.3 rlang_1.0.1       grid_4.1.2        rappdirs_0.3.3    labeling_0.4.2    rmarkdown_2.11   
#> [64] gtable_0.3.0      abind_1.4-5       DBI_1.1.0         curl_4.3.2        R6_2.5.1          extrafont_0.17    utf8_1.2.2       
#> [71] rprojroot_2.0.2   latex2exp_0.5.0   ragg_1.2.0        stringi_1.7.6     Rcpp_1.0.8        vctrs_0.3.8       tidyselect_1.1.2 
#> [78] xfun_0.27
```
