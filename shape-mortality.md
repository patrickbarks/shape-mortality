The shape of mortality trajectories in plants
================
Patrick Barks
2019-01-30

Whereas most animal species exhibit age-related increases in mortality and declines in reproduction (i.e. senescence), most plant species seem to do the opposite. Previous analyses with COMPADRE suggest that only about 1/5 of plant species exhibit age-related increases in mortality. Below I replicate those analyses.

Preliminaries
-------------

#### Install required packages

``` r
devtools::install_github("jschoeley/pash")
devtools::install_github("jonesor/Rcompadre")
devtools::install_github("patrickbarks/Rage", ref = "devel")
install.packages("tidyverse")
```

#### Load required libraries

``` r
library(Rcompadre)
library(Rage)
library(pash)
library(tidyverse)
```

#### Function to calculate the shape of a mortality trajectory

The life table coefficient of variation is one way to measure the 'shape' of a mortality trajectory. Values greater than zero indicate that mortality tends to increase with age (i.e. senescence), whereas values less than zero indicate that mortality tends to decline with age (i.e. negative senescence).

``` r
# calculate life table coefficient of variation
shape_mortality <- function(matU, startLife = 1, N = 1000) {
  lx <- mpm_to_lx(matU, startLife = startLife, N = N)
  len <- max(which(lx > 1e-06))
  if (len < 5) return(NA_real_)
  x <- seq_along(lx) - 1
  ps <- try(pash::Inputlx(x[1:len], lx[1:len]), silent = TRUE)
  if (class(ps) == 'try-error') return(NA_real_)
  out <- pash::GetShape(ps, type = 'cv', harmonized = TRUE)
  return(out)
}
```

Prepare data
------------

#### Fetch COMPADRE and subset to matrices of interest

``` r
# fetch compadre db from web
comp <- cdb_fetch("compadre")

# subset to matrices of interest, and collapse to mean matrix by spp
comp_sub <- comp %>% 
  cdb_flag(c("check_NA_U", "check_zero_U")) %>% 
  filter(check_NA_U == FALSE,
         check_zero_U == FALSE,
         MatrixCaptivity == 'W',
         MatrixTreatment == "Unmanipulated",
         MatrixComposite != "Seasonal",
         MatrixPopulation != "FRANK",
         AnnualPeriodicity == 1,
         OrganismType %in% c('Herbaceous perennial', 'Tree', 'Shrub',
                             'Succulent', 'Palm'),
         MatrixDimension > 2,
         SurvivalIssue <= 1.05) %>% 
  mutate(id_stage = cdb_id_stages(.))
```

#### Collapse to mean matrix population model by species

``` r
comp_mean <- comp_sub %>% 
  cdb_collapse(c("SpeciesAuthor", "id_stage")) %>% 
  cdb_unnest()
```

Analyze the shape of mortality trajectories
-------------------------------------------

``` r
comp_shape <- comp_mean %>% 
  mutate(start = mpm_first_active(.)) %>% 
  mutate(shape = map2_dbl(matU, start, shape_mortality)) %>%
  mutate(id = as.factor(1:n())) %>% 
  filter(!is.na(shape))
```

#### Example mortality trajectories

Below I show examples of two mortality trajectories: a trajectory of declining mortality from the coniferous tree species *Juniperus procera*, and a trajectory of increasing mortality from the tropical tree species *Oxandra asbeckii*.

``` r
example <- comp_shape %>% 
  filter(SpeciesAuthor %in% c("Juniperus_procera", "Oxandra_asbeckii")) %>% 
  mutate(age = list(0:100)) %>% 
  mutate(hazard = map2(matU, start, mpm_to_hx, N = 100)) %>% 
  CompadreData() %>% 
  select(SpeciesAuthor, age, hazard) %>% 
  unnest()

example_cv <- comp_shape %>% 
  filter(SpeciesAuthor %in% c("Juniperus_procera", "Oxandra_asbeckii")) %>% 
  CompadreData() %>% 
  select(SpeciesAuthor, shape) %>% 
  mutate(label = paste0("shape_mortality==", formatC(shape, format = "f", digits = 2)))

ggplot(example) +
  geom_line(aes(age, hazard)) +
  geom_text(data = example_cv, aes(x = 10, y = Inf, label = label),
            hjust = 0, vjust = 1.8, parse = TRUE) +
  facet_wrap(~ SpeciesAuthor, scales = "free_y") +
  labs(x = "Age (years)", y = "Mortality hazard")
```

![](shape-mortality_files/figure-markdown_github/unnamed-chunk-8-1.png)

#### The shape of mortality trajectories across COMPADRE

Here's a plot showing the distribution of mortality trajectory shapes across all species in COMPADRE.

``` r
ggplot(comp_shape, aes(fct_reorder(id, shape))) +
  geom_point(aes(y = shape), size = 2, alpha = 0.3, pch = 21) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_x_discrete(expand = c(0.015, 0)) +
  scale_y_continuous(breaks = seq(-5, 1, 1)) +
  coord_cartesian(ylim = c(-5, 1)) +
  labs(x = "Species (ranked)", y = "Shape mortality (life table CV)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](shape-mortality_files/figure-markdown_github/unnamed-chunk-9-1.png)
