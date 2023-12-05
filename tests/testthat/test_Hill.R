# Prepare object for area-based test
tmat <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
simpsonIndex <- (4/9)^2 + (3/9)^2 + (2/9)^2

test_that("Test equality H'== Hill alpha-->1 sequential", {
	expect_equal(
		as.numeric(Shannon(tmat, na.tolerance=0, np=1)),
		as.numeric(log(Hill(tmat, alpha=1, na.tolerance=0, np=1)$alpha.1))
		)
})

test_that("Test equality H'== Hill alpha-->1 parallel", {
	expect_equal(
		as.numeric(Shannon(tmat, na.tolerance=0, np=2)),
		as.numeric(log(Hill(tmat, alpha=1, na.tolerance=0, np=2)$alpha.1))
		)
})

test_that("Test equality 1/Simpson index == Hill alpha-->2 sequential", {
	expect_equal(
		Hill(tmat, alpha=2, na.tolerance=0, np=1)$alpha.2[2,2],
		1/simpsonIndex
		)
})

test_that("Test equality 1/Simpson index == Hill alpha-->2 parallel", {
	expect_equal(
		Hill(tmat, alpha=2, na.tolerance=0, np=2)$alpha.2[2,2],
		1/simpsonIndex
		)
})