# scripts/R/utils/factor_order.R — the one display order of named factors in the figures
# (9e5dd0: region-plot series and the climate correlogram's axes came text-sorted, bio_1,
# bio_10, ..., bio_2).

test_that("numbered names come first in numeric order, the others after them as given", {
    x <- c("bio_10", "climate_multivariate", "bio_2", "bio_1", "height", "bio_19")
    expect_identical(x[numbered_order(x)],
                     c("bio_1", "bio_2", "bio_10", "bio_19", "climate_multivariate", "height"))
})

test_that("names without a number keep their order", {
    x <- c("height", "flowering_time", "disease_score")
    expect_identical(x[numbered_order(x)], x)
    expect_identical(numbered_order(character(0)), integer(0))
})

test_that("ties keep the input order, and the number must end the name", {
    x <- c("trait_2", "bio_2", "PC1", "bio_2x")
    expect_identical(x[numbered_order(x)], c("trait_2", "bio_2", "PC1", "bio_2x"))
})
