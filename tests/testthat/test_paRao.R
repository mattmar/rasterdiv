# Prepare object for area-based test
testpol <- as.polygons(terra::rast(terra::rast(matrix(c(2,2,2,4,4,4,2,3,4),ncol=3))))
testras <- terra::rast(matrix(c(2,2,2,4,4,4,2,3,4),ncol=3))
testpol[["id"]] <- 1

paramrao <- function(x,size,alpha) 
{
	n <- size^2
	d <- dist(x)
	prop <- 1/(n^2)
	dalpha <- d^(alpha)
	dprop <- dalpha*prop+(dalpha*prop)
	par <- sum(dprop)
	prao <- par^(1/alpha)
	prao
}

tmat <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
tvec <- c(10,10,10,20,20,20,20,30,30)

test_that("paRao uni and multicore many alphas", {
	expect_equal(
		paRao(tmat,dist_m="euclidean",na.tolerance=1,alpha=c(1:5,Inf),np=1),
		paRao(tmat,dist_m="euclidean",na.tolerance=1,alpha=c(1:5,Inf),np=2)
		)
})

test_that("paRao manual and rasterdiv alpha=1", {
	expect_equal(
		paramrao(tvec,size=3,alpha=1),
		as.numeric(unlist(paRao(tmat,alpha=1))[1])
		)
})

test_that("paRao manual and rasterdiv alpha=5", {
	expect_equal(
		paramrao(tvec,size=3,alpha=5),
		as.numeric(unlist(paRao(tmat,alpha=5))[1])
		)
})

test_that("paRao alpha=0 is 0", {
	expect_equal(
		sum(unlist(paRao(tmat,window=3,alpha=0))),
		0
		)
})

test_that("Negative alpha is an error", {
	expect_error(
		paRao(tmat,window=3,alpha=-1))
})

test_that("Rao alpha==2 > Rao alpha==1", {
	expect_gt(
		sum(unlist(paRao(tmat,window=3,alpha=2))),
		sum(unlist(paRao(tmat,window=3,alpha=1)))
		)
})

test_that("multi paRao uni and multicore many alphas", {
	expect_equal(
		paRao(list(tmat,tmat),method="multidimension",na.tolerance=1,alpha=c(1:5,Inf),np=1),
		paRao(list(tmat,tmat),method="multidimension",na.tolerance=1,alpha=c(1:5,Inf),np=2)
		)
})

test_that("Rao alpha==2^31 is equal to Rao alpha==(2^31)-1", {
	expect_equal(
		sum(unlist(paRao(tmat,window=3,alpha=(.Machine$integer.max)-1,rasterOut=FALSE))),
		sum(unlist(paRao(tmat,window=3,alpha=(.Machine$integer.max),rasterOut=FALSE)))
		)
})

test_that("Area-based Rao against human-deriverd results...", {
	expect_equal(
		as.numeric(paRao(x=testras, area = testpol, window = 3, field = 'id', alpha=c(1:5), debugging=FALSE)$alpha.1[1]),
		c(((4/9*4/9)*2 + (1/9*4/9)*1 + (1/9*4/9)*1)*2)
		)
})