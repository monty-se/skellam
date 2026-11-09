## --------------------------------------------------------------------------
## icu_bavaria.R -- builds data/icu_bavaria.rda
##
## Sources (Robert Koch-Institut, licence CC BY 4.0, pinned to the versions
## of 8 and 9 November 2026):
##   - Intensivkapazitaeten und COVID-19-Intensivbettenbelegung in Deutschland
##     (DIVI-Intensivregister), file Intensivregister_Bundeslaender_Kapazitaeten.csv
##     https://github.com/robert-koch-institut/Intensivkapazitaeten_und_COVID-19-Intensivbettenbelegung_in_Deutschland
##   - COVID-19 7-Tage-Inzidenz in Deutschland,
##     file COVID-19-Faelle_7-Tage-Inzidenz_Bundeslaender.csv
##     https://github.com/robert-koch-institut/COVID-19_7-Tage-Inzidenz_in_Deutschland
##
## Changes made: rows for Bavaria (state 09) and all ages only, dates from
## 20 March 2020 to 31 December 2023, four columns kept and renamed.
##
## Run from the root of the clone: Rscript data-raw/icu_bavaria.R
##
## Authors: Montasser Ghachem and Oguz Ersan, 2026
## License: GPL (>= 2)
## --------------------------------------------------------------------------

raw <- "https://raw.githubusercontent.com/robert-koch-institut/"
icu_url <- paste0(raw, "Intensivkapazitaeten_und_COVID-19-Intensivbettenbelegung_in_Deutschland/",
                  "ac3da9bd6b7d42f0f644fe21a96c50f966e5153d/Intensivregister_Bundeslaender_Kapazitaeten.csv")
inz_url <- paste0(raw, "COVID-19_7-Tage-Inzidenz_in_Deutschland/",
                  "8427099a9bb7411918ae9c180a9ecb95de5bf7e4/COVID-19-Faelle_7-Tage-Inzidenz_Bundeslaender.csv")

icu <- utils::read.csv(icu_url, colClasses = c(bundesland_id = "character"))
inz <- utils::read.csv(inz_url, colClasses = c(Bundesland_id = "character"), check.names = FALSE)

icu <- icu[icu$bundesland_id == "09", ]
inz <- inz[inz$Bundesland_id == "09" & inz$Altersgruppe == "00+", ]

icu_bavaria <- data.frame(date = as.Date(icu$datum),
                          occupancy = as.integer(icu$faelle_covid_aktuell),
                          admissions = as.integer(icu$faelle_covid_erstaufnahmen))
icu_bavaria$cases7 <- as.integer(inz[["Faelle_7-Tage"]][match(icu_bavaria$date, as.Date(inz$Meldedatum))])
icu_bavaria <- icu_bavaria[icu_bavaria$date <= as.Date("2023-12-31"), ]
icu_bavaria <- icu_bavaria[order(icu_bavaria$date), ]
rownames(icu_bavaria) <- NULL

stopifnot(!anyNA(icu_bavaria[c("date", "occupancy", "cases7")]),
          all(diff(icu_bavaria$date) == 1))

save(icu_bavaria, file = "data/icu_bavaria.rda", compress = "xz")
