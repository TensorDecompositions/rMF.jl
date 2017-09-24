info("Test rMF")

info("Test23 with concentrations only")
rMF.loaddata("test23")
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test23 with concentrations and deltas")
rMF.loaddata("test23delta")
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test23 with concentrations and ratios")
rMF.loaddata("test23ratio")
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test27 with ratios only")
rMF.loaddata("test27ratiosonly")
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test with concentrations only")
rMF.loaddata("test", nw=5, nc=3, ns=2)
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test with concentrations and deltas")
rMF.loaddata("test", nw=5, nc=3, ns=2, nd=1)
@time rMF.execute(2:3; maxiter=100, tol=1e-2)

info("Test with concentrations and ratios but ignore ratios")
rMF.loaddata("test", nw=5, nc=3, ns=2, nr=1)
@time rMF.execute(2:3; maxiter=100, tol=1e-2, ignoreratios=true)