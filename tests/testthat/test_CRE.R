tmat <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)

test_that("CRE uni and multicore", {
	expect_equal(
		CRE(tmat,window=3,na.tolerance=1,np=1),
		CRE(tmat,window=3,na.tolerance=1,np=2)
		)
})