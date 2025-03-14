FUNCTION rlim(d1,d2,d3)
! CEE'S KAPPA = 1/3 SCHEME
! ensures continuity
  r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
  rlim = (d2-d3)*AMAX1(0.,AMIN1(r,AMIN1(1./6.+1./3.*r,1.)))
  !AMAX1/AMIN1 returns the maximum/min value among two real numbers, respectively. operates on single-precision (real) numbers.
  !AMAX1 ensures that the value is bounded between 0 and 1/3.
  RETURN
END FUNCTION
