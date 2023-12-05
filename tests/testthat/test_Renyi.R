# Prepare object for area-based test
tmat <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)

test_that("Test equality BP == Renyi alpha-->+Inf sequential", {
	expect_equal(
		as.numeric(BergerParker(tmat, na.tolerance=0, np=1)),
		as.numeric(Renyi(tmat, alpha=Inf, na.tolerance=0, np=1)$alpha.Inf)
		)
})

test_that("Test equality BP == Renyi alpha-->+Inf parallel", {
	expect_equal(
		as.numeric(Shannon(tmat, na.tolerance=0, np=2)),
		as.numeric(Renyi(tmat, alpha=1, na.tolerance=0, np=2)$alpha.1)
		)
})