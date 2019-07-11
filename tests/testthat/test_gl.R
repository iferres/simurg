context('Simulating gene gain and loss.')

# library(ape)

ne <- 1e10
ou <- 1e-8
ov <- 1e-10
norg <- 10
theta <- 2000
rho <- 20
mrca_acc_o <- round(theta / rho)
phy <- structure(list(edge = structure(c(11L, 13L, 19L, 19L, 13L, 18L,
                                         18L, 11L, 12L, 15L, 15L, 16L, 16L, 12L, 14L, 14L, 17L, 17L, 13L,
                                         19L, 1L, 2L, 18L, 3L, 4L, 12L, 15L, 5L, 16L, 6L, 7L, 14L, 8L,
                                         17L, 9L, 10L), .Dim = c(18L, 2L)),
                      edge.length = c(0.570743614638502, 0.743533392747197, 0.00315730974266033, 0.00315730974266033,
                                      0.718087960984797, 0.028602741505061, 0.028602741505061, 0.226254712324589,
                                      0.961297501178038, 0.129882103625732, 0.0568438422459096, 0.0730382613798221,
                                      0.0730382613798221, 0.924327430521649, 0.166852174282121, 0.116060928134785,
                                      0.0507912461473357, 0.0507912461473357),
                      tip.label = c("genome10","genome6", "genome5", "genome7", "genome9", "genome4", "genome1",
                                    "genome2", "genome8", "genome3"), Nnode = 9L),
                 class = "phylo", order = "cladewise")

m <- as.data.frame(phy$edge)
m$length <- phy$edge.length

depth <- coalescent.intervals(phy)$total.depth
brti <- c(structure(rep(0, norg), names=1:norg), branching.times(phy))

pm <- simba:::.sim_gl(phy = phy,
                       m = m,
                       brti = brti,
                       ne = ne,
                       norg = norg,
                       ou = ou,
                       ov = ov,
                       mrca_acc_o = mrca_acc_o)

test_that('.sim_gl works',{
  expect_is(pm, 'matrix')
  expect_equivalent(dim(pm)[1], norg)
  expect_identical(dimnames(pm)[[1]], paste0('genome', 1:norg))
  expect_true(all(sapply(pm, is.numeric)))
})
