test_that("icu_bavaria has one row per day and the documented columns", {
  expect_named(icu_bavaria, c("date", "occupancy", "admissions", "cases7"))
  expect_s3_class(icu_bavaria$date, "Date")
  expect_true(all(diff(icu_bavaria$date) == 1))
  expect_false(anyNA(icu_bavaria[c("occupancy", "cases7")]))
  expect_true(all(is.na(icu_bavaria$admissions[icu_bavaria$date < as.Date("2021-07-29")])))
  expect_false(anyNA(icu_bavaria$admissions[icu_bavaria$date >= as.Date("2021-07-29")]))
})
