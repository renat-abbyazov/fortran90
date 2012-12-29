module f90_vs_numpy_arrays_5_qpck

implicit none
private
public qags, qextr, qk21, qsort, timestamp, test04, f03


contains


  subroutine qags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

    !*****************************************************************************80
    !
    !! QAGS estimates the integral of a function.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !
    !    Output, real ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer IER, error flag.
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                         = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more sub-
    !                             divisions by increasing the data value of
    !                             limit in qags (and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement
    !                             it is advised to analyze the integrand
    !                             in order to determine the integration
    !                             difficulties. if the position of a
    !                             local difficulty can be determined (e.g.
    !                             singularity, discontinuity within the
    !                             interval) one will probably gain from
    !                             splitting up the interval at this point
    !                             and calling the integrator on the sub-
    !                             ranges. if possible, an appropriate
    !                             special-purpose integrator should be used,
    !                             which is designed for handling the type
    !                             of difficulty involved.
    !                         = 2 the occurrence of roundoff error is detec-
    !                             ted, which prevents the requested
    !                             tolerance from being achieved.
    !                             the error may be under-estimated.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some  points of the integration
    !                             interval.
    !                         = 4 the algorithm does not converge. roundoff
    !                             error is detected in the extrapolation
    !                             table. it is presumed that the requested
    !                             tolerance cannot be achieved, and that the
    !                             returned result is the best which can be
    !                             obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier.
    !                         = 6 the input is invalid, because
    !                             epsabs < 0 and epsrel < 0,
    !                             result, abserr and neval are set to zero.
    !
    !  Local Parameters:
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least limexp+2 containing
    !                       the part of the epsilon table which is still
    !                       needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !                       (before that subdivision has taken place)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left interval
    !           *****2    - variable for the right interval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements currently in rlist2. if an
    !                       appropriate approximation to the compounded
    !                       integral has been obtained it is put in
    !                       rlist2(numrl2) after numrl2 has been increased
    !                       by one.
    !           small     - length of the smallest interval considered
    !                       up to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine is
    !                       attempting to perform extrapolation i.e. before
    !                       subdividing the smallest interval we try to
    !                       decrease the value of erlarg.
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true value)
    !
    implicit none

    integer, parameter :: limit = 500

    real a
    real abseps
    real abserr
    real alist(limit)
    real area
    real area1
    real area12
    real area2
    real a1
    real a2
    real b
    real blist(limit)
    real b1
    real b2
    real correc
    real defabs
    real defab1
    real defab2
    real dres
    real elist(limit)
    real epsabs
    real epsrel
    real erlarg
    real erlast
    real errbnd
    real errmax
    real error1
    real error2
    real erro12
    real errsum
    real ertest
    logical extrap
    real, external :: f
    integer id
    integer ier
    integer ierro
    integer iord(limit)
    integer iroff1
    integer iroff2
    integer iroff3
    integer jupbnd
    integer k
    integer ksgn
    integer ktmin
    integer last
    logical noext
    integer maxerr
    integer neval
    integer nres
    integer nrmax
    integer numrl2
    real resabs
    real reseps
    real result
    real res3la(3)
    real rlist(limit)
    real rlist2(52)
    real small

!interface
!    real function f(x)
!    implicit none
!    real, intent(in) :: x
!    end function
!end interface

    !
    !  The dimension of rlist2 is determined by the value of
    !  limexp in QEXTR (rlist2 should be of dimension
    !  (limexp+2) at least).
    !
    !  Test on validity of parameters.
    !
    ier = 0
    neval = 0
    last = 0
    result = 0.0e+00
    abserr = 0.0e+00
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0e+00
    elist(1) = 0.0e+00

    if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
       ier = 6
       return
    end if
    !
    !  First approximation to the integral.
    !
    ierro = 0
    call qk21 ( f, a, b, result, abserr, defabs, resabs )
    !
    !  Test on accuracy.
    !
    dres = abs ( result )
    errbnd = max ( epsabs, epsrel * dres )
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1

    if ( abserr <= 1.0e+02 * epsilon ( defabs ) * defabs .and. &
         abserr > errbnd ) then
       ier = 2
    end if

    if ( limit == 1 ) then
       ier = 1
    end if

    if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
         abserr == 0.0e+00 ) go to 140
    !
    !  Initialization.
    !
    rlist2(1) = result
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    abserr = huge ( abserr )
    nrmax = 1
    nres = 0
    numrl2 = 2
    ktmin = 0
    extrap = .false.
    noext = .false.
    iroff1 = 0
    iroff2 = 0
    iroff3 = 0

    if ( dres >= (1.0e+00-5.0e+01* epsilon ( defabs ) ) * defabs ) then
       ksgn = 1
    else
       ksgn = -1
    end if

    do last = 2, limit
       !
       !  Bisect the subinterval with the nrmax-th largest error estimate.
       !
       a1 = alist(maxerr)
       b1 = 5.0e-01 * ( alist(maxerr) + blist(maxerr) )
       a2 = b1
       b2 = blist(maxerr)
       erlast = errmax
       call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
       call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
       !
       !  Improve previous approximations to integral and error
       !  and test for accuracy.
       !
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)

       if ( defab1 == error1 .or. defab2 == error2 ) go to 15

       if ( abs ( rlist(maxerr) - area12) > 1.0e-05 * abs(area12) &
            .or. erro12 < 9.9e-01 * errmax ) go to 10

       if ( extrap ) then
          iroff2 = iroff2+1
       else
          iroff1 = iroff1+1
       end if

10     continue

       if ( last > 10 .and. erro12 > errmax ) then
          iroff3 = iroff3+1
       end if

15     continue

       rlist(maxerr) = area1
       rlist(last) = area2
       errbnd = max ( epsabs, epsrel*abs(area) )
       !
       !  Test for roundoff error and eventually set error flag.
       !
       if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
          ier = 2
       end if

       if ( iroff2 >= 5 ) then
          ierro = 3
       end if
       !
       !  Set error flag in the case that the number of subintervals
       !  equals limit.
       !
       if ( last == limit ) then
          ier = 1
       end if
       !
       !  Set error flag in the case of bad integrand behavior
       !  at a point of the integration range.
       !
       if ( max ( abs(a1),abs(b2)) <= (1.0e+00+1.0e+03* epsilon ( a1 ) )* &
            (abs(a2)+1.0e+03* tiny ( a2 ) ) ) then
          ier = 4
       end if
       !
       !  Append the newly-created intervals to the list.
       !
       if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
       else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
       end if
       !
       !  Call QSORT to maintain the descending ordering
       !  in the list of error estimates and select the subinterval
       !  with nrmax-th largest error estimate (to be bisected next).
       !
       call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

       if ( errsum <= errbnd ) go to 115

       if ( ier /= 0 ) then
          exit
       end if

       if ( last == 2 ) go to 80
       if ( noext ) go to 90

       erlarg = erlarg-erlast

       if ( abs(b1-a1) > small ) then
          erlarg = erlarg+erro12
       end if
       !
       !  Test whether the interval to be bisected next is the
       !  smallest interval.
       !
       if ( .not. extrap ) then
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
          extrap = .true.
          nrmax = 2
       end if

       !40  continue
       !
       !  The smallest interval has the largest error.
       !  Before bisecting decrease the sum of the errors over the
       !  larger intervals (erlarg) and perform extrapolation.
       !
       if ( ierro /= 3 .and. erlarg > ertest ) then

          id = nrmax
          jupbnd = last

          if ( last > (2+limit/2) ) then
             jupbnd = limit+3-last
          end if

          do k = id, jupbnd
             maxerr = iord(nrmax)
             errmax = elist(maxerr)
             if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
                go to 90
             end if
             nrmax = nrmax+1
          end do

       end if
       !
       !  Perform extrapolation.
       !
       !60  continue

       numrl2 = numrl2+1
       rlist2(numrl2) = area
       call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
       ktmin = ktmin+1

       if ( ktmin > 5 .and. abserr < 1.0e-03 * errsum ) then
          ier = 5
       end if

       if ( abseps < abserr ) then

          ktmin = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs,epsrel*abs(reseps))

          if ( abserr <= ertest ) then
             exit
          end if

       end if
       !
       !  Prepare bisection of the smallest interval.
       !
       if ( numrl2 == 1 ) then
          noext = .true.
       end if

       if ( ier == 5 ) then
          exit
       end if

       maxerr = iord(1)
       errmax = elist(maxerr)
       nrmax = 1
       extrap = .false.
       small = small * 5.0e-01
       erlarg = errsum
       go to 90

80     continue

       small = abs ( b - a ) * 3.75e-01
       erlarg = errsum
       ertest = errbnd
       rlist2(2) = area

90     continue

    end do
    !
    !  Set final result and error estimate.
    !
    if ( abserr == huge ( abserr ) ) then
       go to 115
    end if

    if ( ier + ierro == 0 ) then
       go to 110
    end if

    if ( ierro == 3 ) then
       abserr = abserr + correc
    end if

    if ( ier == 0 ) then
       ier = 3
    end if

    if ( result /= 0.0e+00.and.area /= 0.0e+00 ) then
       go to 105
    end if

    if ( abserr > errsum ) go to 115
    if ( area == 0.0e+00 ) go to 130
    go to 110

105 continue

    if ( abserr/abs(result) > errsum/abs(area) ) go to 115
    !
    !  Test on divergence.
    !
110 continue

    if ( ksgn == (-1).and.max ( abs(result),abs(area)) <=  &
         defabs*1.0e-02 ) go to 130

    if ( 1.0e-02 > (result/area) .or. (result/area) > 1.0e+02 &
         .or. errsum > abs(area) ) then
       ier = 6
    end if

    go to 130
    !
    !  Compute global integral sum.
    !
115 continue

    result = sum ( rlist(1:last) )

    abserr = errsum

130 continue

    if ( 2 < ier ) then
       ier = ier - 1
    end if

140 continue

    neval = 42*last-21

    return
  end subroutine qags


  subroutine qk21 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real F, the name of the function routine, of the form
    !      function f ( x )
    !      real f
    !      real x
    !    which evaluates the integrand function.
    !
    !    Input, real A, B, the limits of integration.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 21-point Kronrod rule (resk) 
    !    obtained by optimal addition of abscissae to the 10-point Gauss 
    !    rule (resg).
    !
    !    Output, real ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    implicit none

    real a
    real absc
    real abserr
    real b
    real centr
    real dhlgth
    real, external :: f
    real fc
    real fsum
    real fval1
    real fval2
    real fv1(10)
    real fv2(10)
    real hlgth
    integer j
    integer jtw
    integer jtwm1
    real resabs
    real resasc
    real resg
    real resk
    real reskh
    real result
    real wg(5)
    real wgk(11)
    real xgk(11)
! interface
!    real function f(x)
!    implicit none
!    real, intent(in) :: x
!    end function
!end interface
   !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 21-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 10-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 10-point Gauss rule
    !
    !           wgk    - weights of the 21-point Kronrod rule
    !
    !           wg     - weights of the 10-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11)/ &
         9.956571630258081e-01,     9.739065285171717e-01, &
         9.301574913557082e-01,     8.650633666889845e-01, &
         7.808177265864169e-01,     6.794095682990244e-01, &
         5.627571346686047e-01,     4.333953941292472e-01, &
         2.943928627014602e-01,     1.488743389816312e-01, &
         0.000000000000000e+00/
    !
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11)/ &
         1.169463886737187e-02,     3.255816230796473e-02, &
         5.475589657435200e-02,     7.503967481091995e-02, &
         9.312545458369761e-02,     1.093871588022976e-01, &
         1.234919762620659e-01,     1.347092173114733e-01, &
         1.427759385770601e-01,     1.477391049013385e-01, &
         1.494455540029169e-01/
    !
    data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
         6.667134430868814e-02,     1.494513491505806e-01, &
         2.190863625159820e-01,     2.692667193099964e-01, &
         2.955242247147529e-01/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 10-point Gauss formula
    !           resk   - result of the 21-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 21-point Kronrod approximation to the
    !  integral, and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr)
    resk = wgk(11)*fc
    resabs = abs(resk)

    do j = 1, 5
       jtw = 2*j
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 5
       jtwm1 = 2*j-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(11)*abs(fc-reskh)

    do j = 1, 10
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk21

  subroutine qextr ( n, epstab, result, abserr, res3la, nres )

    !*****************************************************************************80
    !
    !! QEXTR carries out the Epsilon extrapolation algorithm.
    !
    !  Discussion:
    !
    !    The routine determines the limit of a given sequence of approximations, 
    !    by means of the epsilon algorithm of P. Wynn.  An estimate of the 
    !    absolute error is also given.  The condensed epsilon table is computed.
    !    Only those elements needed for the computation of the next diagonal
    !    are preserved.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer N, indicates the entry of EPSTAB which contains
    !    the new element in the first column of the epsilon table.
    !
    !    Input/output, real EPSTAB(52), the two lower diagonals of the triangular
    !    epsilon table.  The elements are numbered starting at the right-hand 
    !    corner of the triangle.
    !
    !    Output, real RESULT, the estimated value of the integral.
    !
    !    Output, real ABSERR, estimate of the absolute error computed from
    !    RESULT and the 3 previous results.
    !
    !    ?, real RES3LA(3), the last 3 results.
    !
    !    Input/output, integer NRES, the number of calls to the routine.  This
    !    should be zero on the first call, and is automatically updated
    !    before return.
    !
    !  Local Parameters:
    !
    !           e0     - the 4 elements on which the
    !           e1       computation of a new element in
    !           e2       the epsilon table is based
    !           e3                 e0
    !                        e3    e1    new
    !                              e2
    !           newelm - number of elements to be computed in the new
    !                    diagonal
    !           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
    !           result - the element in the new diagonal with least value
    !                    of error
    !           limexp is the maximum number of elements the epsilon table
    !           can contain. if this number is reached, the upper diagonal
    !           of the epsilon table is deleted.
    !
    implicit none

    real abserr
    real delta1
    real delta2
    real delta3
    real epsinf
    real epstab(52)
    real error
    real err1
    real err2
    real err3
    real e0
    real e1
    real e1abs
    real e2
    real e3
    integer i
    integer ib
    integer ib2
    integer ie
    integer indx
    integer k1
    integer k2
    integer k3
    integer limexp
    integer n
    integer newelm
    integer nres
    integer num
    real res
    real result
    real res3la(3)
    real ss
    real tol1
    real tol2
    real tol3

    nres = nres+1
    abserr = huge ( abserr )
    result = epstab(n)

    if ( n < 3 ) then
       abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
       return
    end if

    limexp = 50
    epstab(n+2) = epstab(n)
    newelm = (n-1)/2
    epstab(n) = huge ( epstab(n) )
    num = n
    k1 = n

    do i = 1, newelm

       k2 = k1-1
       k3 = k1-2
       res = epstab(k1+2)
       e0 = epstab(k3)
       e1 = epstab(k2)
       e2 = res
       e1abs = abs(e1)
       delta2 = e2-e1
       err2 = abs(delta2)
       tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
       delta3 = e1-e0
       err3 = abs(delta3)
       tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
       !
       !  If e0, e1 and e2 are equal to within machine accuracy, convergence 
       !  is assumed.
       !
       if ( err2 <= tol2 .and. err3 <= tol3 ) then
          result = res
          abserr = err2+err3
          abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
          return
       end if

       e3 = epstab(k1)
       epstab(k1) = e1
       delta1 = e1-e3
       err1 = abs(delta1)
       tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
       !
       !  If two elements are very close to each other, omit a part
       !  of the table by adjusting the value of N.
       !
       if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

       ss = 1.0e+00/delta1+1.0e+00/delta2-1.0e+00/delta3
       epsinf = abs ( ss*e1 )
       !
       !  Test to detect irregular behavior in the table, and
       !  eventually omit a part of the table adjusting the value of N.
       !
       if ( epsinf > 1.0e-04 ) go to 30

20     continue

       n = i+i-1
       exit
       !
       !  Compute a new element and eventually adjust the value of RESULT.
       !
30     continue

       res = e1+1.0e+00/ss
       epstab(k1) = res
       k1 = k1-2
       error = err2+abs(res-e2)+err3

       if ( error <= abserr ) then
          abserr = error
          result = res
       end if

    end do
    !
    !  Shift the table.
    !
    if ( n == limexp ) then
       n = 2*(limexp/2)-1
    end if

    if ( (num/2)*2 == num ) then
       ib = 2
    else
       ib = 1
    end if

    ie = newelm+1

    do i = 1, ie
       ib2 = ib+2
       epstab(ib) = epstab(ib2)
       ib = ib2
    end do

    if ( num /= n ) then

       indx = num-n+1

       do i = 1, n
          epstab(i)= epstab(indx)
          indx = indx+1
       end do

    end if

    if ( nres < 4 ) then
       res3la(nres) = result
       abserr = huge ( abserr )
    else
       abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
            +abs(result-res3la(1))
       res3la(1) = res3la(2)
       res3la(2) = res3la(3)
       res3la(3) = result
    end if

    abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))

    return
  end subroutine qextr

  subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

    !*****************************************************************************80
    !
    !! QSORT maintains the order of a list of local error estimates.
    !
    !  Discussion:
    !
    !    This routine maintains the descending ordering in the list of the 
    !    local error estimates resulting from the interval subdivision process. 
    !    At each call two error estimates are inserted using the sequential 
    !    search top-down for the largest error estimate and bottom-up for the
    !    smallest error estimate.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer LIMIT, the maximum number of error estimates the list can
    !    contain.
    !
    !    Input, integer LAST, the current number of error estimates.
    !
    !    Input/output, integer MAXERR, the index in the list of the NRMAX-th 
    !    largest error.
    !
    !    Output, real ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
    !
    !    Input, real ELIST(LIMIT), contains the error estimates.
    !
    !    Input/output, integer IORD(LAST).  The first K elements contain 
    !    pointers to the error estimates such that ELIST(IORD(1)) through
    !    ELIST(IORD(K)) form a decreasing sequence, with
    !      K = LAST 
    !    if 
    !      LAST <= (LIMIT/2+2), 
    !    and otherwise
    !      K = LIMIT+1-LAST.
    !
    !    Input/output, integer NRMAX.
    !
    implicit none

    integer last

    real elist(last)
    real ermax
    real errmax
    real errmin
    integer i
    integer ibeg
    integer iord(last)
    integer isucc
    integer j
    integer jbnd
    integer jupbn
    integer k
    integer limit
    integer maxerr
    integer nrmax
    !
    !  Check whether the list contains more than two error estimates.
    !
    if ( last <= 2 ) then
       iord(1) = 1
       iord(2) = 2
       go to 90
    end if
    !
    !  This part of the routine is only executed if, due to a
    !  difficult integrand, subdivision increased the error
    !  estimate. in the normal case the insert procedure should
    !  start after the nrmax-th largest error estimate.
    !
    errmax = elist(maxerr)

    do i = 1, nrmax-1

       isucc = iord(nrmax-1)

       if ( errmax <= elist(isucc) ) then
          exit
       end if

       iord(nrmax) = isucc
       nrmax = nrmax-1

    end do
    !
    !  Compute the number of elements in the list to be maintained
    !  in descending order.  This number depends on the number of
    !  subdivisions still allowed.
    !
    jupbn = last

    if ( (limit/2+2) < last ) then
       jupbn = limit+3-last
    end if

    errmin = elist(last)
    !
    !  Insert errmax by traversing the list top-down, starting
    !  comparison from the element elist(iord(nrmax+1)).
    !
    jbnd = jupbn-1
    ibeg = nrmax+1

    do i = ibeg, jbnd
       isucc = iord(i)
       if ( elist(isucc) <= errmax ) then
          go to 60
       end if
       iord(i-1) = isucc
    end do

    iord(jbnd) = maxerr
    iord(jupbn) = last
    go to 90
    !
    !  Insert errmin by traversing the list bottom-up.
    !
60  continue

    iord(i-1) = maxerr
    k = jbnd

    do j = i, jbnd
       isucc = iord(k)
       if ( errmin < elist(isucc) ) then
          go to 80
       end if
       iord(k+1) = isucc
       k = k-1
    end do

    iord(i) = last
    go to 90

80  continue

    iord(k+1) = last
    !
    !  Set maxerr and ermax.
    !
90  continue

    maxerr = iord(nrmax)
    ermax = elist(maxerr)

    return
  end subroutine qsort

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp

subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests QAGS.
!
!  Discussion:
!
!    QAGS is an adaptive integrator for endpoint singularities.
!
!    integrate log(x)/sqrt(x) from 0 to 1.
!
!    The exact answer is -4.
!
  implicit none

  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f03
  integer ier
  integer neval
  real result
  real, parameter :: true = -4.0E+00

  call qags ( f03, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test QAGS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is LOG(X)/SQRT(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end subroutine test04

real function f03 ( x )

!*****************************************************************************80
!
!! F03 is the integrand function LOG(X)/SQRT(X).
!
  implicit none

    real, intent(in) :: x

  if ( x <= 0.0E+00 ) then
    f03 = 0.0E+00
  else
    f03 = log ( x ) / sqrt ( x )
  end if

  return
end function f03

end module f90_vs_numpy_arrays_5_qpck
