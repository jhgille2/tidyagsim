# Extract the parameter set of a spec

The parameter set of the arguments marked
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html), with
default ranges from the spec. Change ranges with
[`update()`](https://rdrr.io/r/stats/update.html), exactly as in
tidymodels, and pass the result to a tuner as `param_info`.

## Usage

``` r
# S3 method for class 'sim_spec'
extract_parameter_set_dials(x, ..., grid = NULL)
```

## Arguments

- x:

  A simulation spec.

- ...:

  Not used.

- grid:

  Optional data frame of designs; tuned arguments without a default
  range take theirs from its values.

## Value

A dials `parameters` object.

## Examples

``` r
spec <- mapping_population(type = "RIL", lines = tune(), generation = tune())
params <- extract_parameter_set_dials(spec)
params
#> Collection of 2 parameters for tuning
#> 
#>  identifier       type    object
#>       lines      lines nparam[+]
#>  generation generation nparam[+]
#> 
update(params, lines = n_lines(c(100, 400)))
#> Collection of 2 parameters for tuning
#> 
#>  identifier       type    object
#>       lines      lines nparam[+]
#>  generation generation nparam[+]
#> 
```
