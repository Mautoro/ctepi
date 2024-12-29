# `ctepi`
> Cumulative Treatment Effect (CTE) with partial identification (PI) on R.


Partial identification bounds for the Average Treatment Effect (ATE) using
the Cumulative Treatment Effect (CTE) estimator. Empirical cumulative distribution
function (ECDF) with partial identification bounds under different assumptions.
Bounds of Williamson and Downs (1990) for the distribution of the difference between
two random variables.

- <a href="#arrow_down-installation"
  id="toc-arrow_down-installation">:arrow_down: Installation</a>
- <a href="#link-quick-links"
  id="toc-link-quick-links">:link: Quick links
  </a>
- <a href="#blue_book-license"
  id="toc-blue_book-license"><strong>:blue_book:</strong> License</a>
- <a href="#rocket-examples" id="toc-rocket-examples">:rocket: Examples</a>


## :arrow_down: Installation

### Install from GitHub

#### Dependencies:

``` r
install.packages("devtools")
install.packages("Rcpp")
```
For Windows users: You must also install `Rtools`. [install Rtools.](https://cran.r-project.org/bin/windows/Rtools/)

#### Installation:

``` r
devtools::install_github("Mautoro/ctepi")
```

## :link: Quick links

### I would like to report a bug

Head to the [ctepi issue tracker](https://github.com/Mautoro/ctepi/issues).


## **:blue_book:** License

GPL 3.0

## :rocket: Examples

For this illustration, we will use the National Health and Nutrition Examination Survey Data | Epidemiologic Follow-up Study (NHEFS). The dataset will be loaded from the `causaldata` package.

``` r
install.packages("causaldata")
nhefs <- causaldata::nhefs
```

The outcome of interest is the patients' weight, measured in the 1st questionnaire (1971) (`wt71`) and in 1982 (`wt82`). The variable `qsmk` is 1 if patients quit smoking between the 1st questionnaire and 1982, and 0 otherwise.

``` r
nhefs$wt71[nhefs$qsmk==0] |> ecdfPI() |> plot(col="red", col.bounds = "red", 
                                              main="Weight (kg), 1971", 
                                              xlim = range(c(nhefs$wt71,nhefs$wt82),na.rm = T) )
nhefs$wt71[nhefs$qsmk==1] |> ecdfPI() |> plot(col="blue",add=T)
legend( "bottom", legend = c("Quit smoking","Keep smoking"), col = c("blue", "red"), 
        lty=1, cex = 0.8 ,inset = c(0,1) , xpd=TRUE, horiz=T, bty="n")
```

``` r
nhefs$wt82[nhefs$qsmk==0] |> ecdfPI() |> plot(col="red", col.bounds = "red", main="Weight (kg), 1982",
                                              xlim = range(c(nhefs$wt71,nhefs$wt82),na.rm = T) )
nhefs$wt82[nhefs$qsmk==1] |> ecdfPI() |> plot(col="blue",add=T)
legend( "bottom", legend = c("Quit smoking","Keep smoking"), col = c("blue", "red"), 
        lty=1, cex = 0.8 ,inset = c(0,1) , xpd=TRUE, horiz=T, bty="n")
```






