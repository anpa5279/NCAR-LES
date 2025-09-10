FUNCTION rlim(d1, d2, d3)
! CEE'S KAPPA = 1/3 SCHEME

    r = (d1 - d2 + 1.e-100) / (d2 - d3 - 1.e-100)
    rlim = (d2 - d3) * AMAX1(0., AMIN1(r, AMIN1(1./6.+1./3.*r, 1.)))

    RETURN
END FUNCTION
