# A quantitative parameter for any argument

A shortcut for
[`dials::new_quant_param()`](https://dials.tidymodels.org/reference/new-param.html)
used by
[`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md)
and
[`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md),
where arguments have no built-in parameter function.

## Usage

``` r
sim_param(range, type = NULL, label = "value")

sim_param_values(values, label = "value")
```

## Arguments

- range:

  Inclusive lower and upper bounds.

- type:

  `"integer"` or `"double"`.

- label:

  Display label.

- values:

  Allowed values.

## Value

A dials `param` object.

## Examples

``` r
sim_param(c(1, 10))
#> value (quantitative)
#> Range: [1, 10]
```
