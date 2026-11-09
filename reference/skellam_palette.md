# Colours and Theme for skellam Plots

The colour roles and the ggplot2 theme used by the package's plots. The
colours come from the Okabe-Ito palette, which stays distinguishable for
the common forms of colour blindness and in greyscale; colour is never
the only cue (observed and fitted values also differ in shape and line
type).

## Usage

``` r
skellam_palette()

theme_skellam(base_size = 11, base_family = "")

scale_colour_skellam(...)

scale_fill_skellam(...)
```

## Arguments

- base_size:

  Base font size, in points.

- base_family:

  Base font family.

- ...:

  Passed to
  [`scale_colour_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  or
  [`scale_fill_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

`skellam_palette()` returns a named character vector of colours
(`observed`, `fitted`, `reference`, `secondary`, `neutral`, and `band`,
a light tint of the reference colour for reference bands; it is opaque
so figures print and export the same on every graphics device);
`theme_skellam()` a ggplot2 theme; the scale functions ggplot2 scales.

## Examples

``` r
skellam_palette()
#>  observed    fitted reference secondary   neutral      band 
#> "#4D4D4D" "#0072B2" "#D55E00" "#009E73" "#999999" "#F7DCCB" 
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  d <- data.frame(x = -6:6, p = dskellam(-6:6, 3, 2))
  ggplot(d, aes(x, p)) +
    geom_col(fill = skellam_palette()["fitted"]) +
    theme_skellam()
}
```
