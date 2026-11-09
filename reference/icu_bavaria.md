# COVID-19 Intensive Care Occupancy in Bavaria

Daily number of COVID-19 patients in adult intensive care units in
Bavaria, Germany, with the reported first admissions and the infections
reported in the previous seven days, from 20 March 2020 to 31 December
2023.

## Usage

``` r
icu_bavaria
```

## Format

A data frame with 1382 rows, one per day, and 4 columns:

- `date`:

  Date of the report.

- `occupancy`:

  COVID-19 patients in adult intensive care (`faelle_covid_aktuell` in
  the register).

- `admissions`:

  First admissions of COVID-19 patients to intensive care that day
  (`faelle_covid_erstaufnahmen`); `NA` before 29 July 2021, when the
  register started recording them. Patients moved between intensive care
  units are not first admissions.

- `cases7`:

  Laboratory-confirmed SARS-CoV-2 infections in Bavaria, all ages,
  reported in the seven days up to and including `date`
  (`Faelle_7-Tage`).

## Source

Robert Koch-Institut (2026). Intensivkapazitaeten und
COVID-19-Intensivbettenbelegung in Deutschland, version of 8 November
2026. Zenodo,
[doi:10.5281/zenodo.24053879](https://doi.org/10.5281/zenodo.24053879) .
Robert Koch-Institut (2026). COVID-19 7-Tage-Inzidenz in Deutschland,
version of 9 November 2026. Zenodo,
[doi:10.5281/zenodo.24074082](https://doi.org/10.5281/zenodo.24074082) .
Both are licensed under CC BY 4.0
(<https://creativecommons.org/licenses/by/4.0/>). Changes: rows for
Bavaria and all ages only, dates up to 31 December 2023, four columns
kept and renamed. The script `data-raw/icu_bavaria.R` in the package
repository rebuilds the data from the sources.

## Details

The day-to-day change in `occupancy` is the number of patients who
entered intensive care minus the number who left it (discharged, moved
to a normal ward, or died), two counts that were not published for most
of the pandemic. Modelled as independent Poisson counts, their
difference has a Skellam distribution (Rave and Kauermann, 2025). Since
29 July 2021 the register also reports first admissions, which can be
used to check the admissions estimated from the occupancy alone; see the
article [Disease
surveillance](https://monty-se.github.io/skellam/articles/surveillance.html)
on the package website.

## References

Rave, M. and Kauermann, G. (2025). The Skellam distribution revisited:
Estimating the unobserved incoming and outgoing ICU COVID-19 patients on
a regional level in Germany. *Statistical Modelling* **25**(3), 270-288.
[doi:10.1177/1471082X241235024](https://doi.org/10.1177/1471082X241235024)

## Examples

``` r
head(icu_bavaria)
#>         date occupancy admissions cases7
#> 1 2020-03-20        61         NA   2967
#> 2 2020-03-21        67         NA   3580
#> 3 2020-03-22        82         NA   4459
#> 4 2020-03-23        94         NA   4851
#> 5 2020-03-24       123         NA   5305
#> 6 2020-03-25       132         NA   6048

# the daily change in occupancy: admissions minus departures
change <- diff(icu_bavaria$occupancy)
c(mean = mean(change), variance = var(change))
#>         mean     variance 
#>   0.07892831 157.83072233 

# both rates over the first quarter of 2022, from the changes alone
z <- subset(icu_bavaria, date >= as.Date("2022-01-01") & date < as.Date("2022-04-01"))
fit <- skellam.fit(diff(z$occupancy))
coef(fit)
#>  lambda1  lambda2 
#> 102.5985 105.2389 
mean(z$admissions)   # reported first admissions per day
#> [1] 37.36667
# The estimated rates are far higher: without covariates, every source of
# day-to-day variation (transfers between units, reporting delays) counts as
# admissions and departures. The article shows how a model in which
# departures are proportional to occupancy recovers the reported admissions.
```
