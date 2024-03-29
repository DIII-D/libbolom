c     algorithm 587, collected algorithms from acm.
c     algorithm appeared in acm-trans. math. software, vol.8, no. 3,
c     sep., 1982, p. 323.


C 20060607 tbt Changed Dimension of in put arguments from (1) with (*)


      subroutine lsei(w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml,
     *                mode, ws, ip)


c     dimension w(mdw,n+1),prgopt(*),x(n),
c     ws(2*(me+n)+k+(mg+2)*(n+7)),ip(mg+2*n+2)
c     above, k=max(ma+mg,n).
c
c     abstract
c
c     this subprogram solves a linearly constrained least squares
c     problem with both equality and inequality constraints, and, if the
c     user requests, obtains a covariance matrix of the solution
c     parameters.
c
c     suppose there are given matrices e, a and g of respective
c     dimensions me by n, ma by n and mg by n, and vectors f, b and h of
c     respective lengths me, ma and mg.  this subroutine solves the
c     linearly constrained least squares problem
c
c                   ex = f, (e me by n) (equations to be exactly
c                                       satisfied)
c                   ax = b, (a ma by n) (equations to be
c                                       approximately satisfied,
c                                       least squares sense)
c                   gx.ge.h,(g mg by n) (inequality constraints)
c
c     the inequalities gx.ge.h mean that every component of the product
c     gx must be .ge. the corresponding component of h.
c
c     in case the equality constraints cannot be satisfied, a
c     generalized inverse solution residual vector length is obtained
c     for f-ex. this is the minimal length possible for f-ex.
c
c
c     any values me.ge.0, ma.ge.0, or mg.ge.0 are permitted.  the
c     rank of the matrix e is estimated during the computation. we call
c     this value kranke. it is an output parameter in ip(1) defined
c     below. using a generalized inverse solution of ex=f, a reduced
c     least squares problem with inequality constraints is obtained.
c     the tolerances used in these tests for determining the rank
c     of e and the rank of the reduced least squares problem are
c     given in sandia tech. rept. sand 78-1290. they can be
c     modified by the user if new values are provided in
c     the option list of the array prgopt(*).
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start editing at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/sasum/dasum/,/sdot/ddot/,
c     /snrm2/dnrm2/,/ sqrt/ dsqrt/,/ abs/ dabs/,/amax1/dmax1/,
c     /scopy/dcopy/,/sscal/dscal/,/saxpy/daxpy/,/sswap/dswap/,/e0/d0/,
c     /, dummy/,sngl(dummy)/,/srelpr/drelpr/
c
c     written by r. j. hanson and k. h. haskell.  for further math.
c     and algorithmic details see sandia laboratories tech. repts.
c     sand 77-0552, (1978), sand 78-1290, (1979), and
c     math. programming, vol. 21, (1981), p.98-118.
c
c     the user must dimension all arrays appearing in the call list..
c     w(mdw,n+1),prgopt(*),x(n),ws(2*(me+n)+k+(mg+2)*(n+7)),ip(mg+2*n+2)
c     where k=max(ma+mg,n).  this allows for a solution of a range of
c     problems in the given working space.  the dimension of ws(*)
c     given is a necessary overestimate.  once a particular problem
c     has been run, the output parameter ip(3) gives the actual
c     dimension required for that problem.
c
c     the parameters for lsei( ) are
c
c     input..
c
c     w(*,*),mdw,   the array w(*,*) is doubly subscripted with
c     me,ma,mg,n    first dimensioning parameter equal to mdw.
c                   for this discussion let us call m = me+ma+mg.  then
c                   mdw must satisfy mdw.ge.m.  the condition
c                   mdw.lt.m is an error.
c
c                   the array w(*,*) contains the matrices and vectors
c
c                                  (e  f)
c                                  (a  b)
c                                  (g  h)
c
c                   in rows and columns 1,...,m and 1,...,n+1
c                   respectively.
c
c                   the integers me, ma, and mg are the
c                   respective matrix row dimensions
c                   of e, a and g. each matrix has n columns.
c
c     prgopt(*)    this array is the option vector.
c                  if the user is satisfied with the nominal
c                  subprogram features set
c
c                  prgopt(1)=1 (or prgopt(1)=1.0)
c
c                  otherwise prgopt(*) is a linked list consisting of
c                  groups of data of the following form
c
c                  link
c                  key
c                  data set
c
c                  the parameters link and key are each one word.
c                  the data set can be comprised of several words.
c                  the number of items depends on the value of key.
c                  the value of link points to the first
c                  entry of the next group of data within
c                  prgopt(*).  the exception is when there are
c                  no more options to change.  in that
c                  case link=1 and the values key and data set
c                  are not referenced. the general layout of
c                  prgopt(*) is as follows.
c
c               ...prgopt(1)=link1 (link to first entry of next group)
c               .  prgopt(2)=key1 (key to the option change)
c               .  prgopt(3)=data value (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...prgopt(link1)=link2 (link to the first entry of
c               .                       next group)
c               .  prgopt(link1+1)=key2 (key to the option change)
c               .  prgopt(link1+2)=data value
c               ...     .
c               .       .
c               .       .
c               ...prgopt(link)=1 (no more options to change)
c
c                  values of link that are nonpositive are errors.
c                  a value of link.gt.nlink=100000 is also an error.
c                  this helps prevent using invalid but positive
c                  values of link that will probably extend
c                  beyond the program limits of prgopt(*).
c                  unrecognized values of key are ignored.  the
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  to prevent cycling in the
c                  processing of the option array a count of the
c                  number of options changed is maintained.
c                  whenever this count exceeds nopt=1000 an error
c                  message is printed and the subprogram returns.
c
c                  options..
c
c                  key=1
c                         compute in w(*,*) the n by n
c                  covariance matrix of the solution variables
c                  as an output parameter.  nominally the
c                  covariance matrix will not be computed.
c                  (this requires no user input.)
c                  the data set for this option is a single value.
c                  it must be nonzero when the covariance matrix
c                  is desired.  if it is zero, the covariance
c                  matrix is not computed.  when the covariance matrix
c                  is computed, the first dimensioning parameter
c                  of the array w(*,*) must satisfy mdw.ge.max0(m,n).
c
c                  key=2
c                         scale the nonzero columns of the
c                         entire data matrix.
c                  (e)
c                  (a)
c                  (g)
c
c                  to have length one.  the data set for this
c                  option is a single value.  it must be
c                  nonzero if unit length column scaling
c                  is desired.
c
c                  key=3
c                         scale columns of the entire data matrix
c                  (e)
c                  (a)
c                  (g)
c
c                  with a user-provided diagonal matrix.
c                  the data set for this option consists
c                  of the n diagonal scaling factors, one for
c                  each matrix column.
c
c                  key=4
c                         change the rank determination tolerance for
c                  the equality constraint equations from
c                  the nominal value of sqrt(srelpr).  this quantity can
c                  be no smaller than srelpr, the arithmetic-
c                  storage precision.  the quantity srelpr is the
c                  largest positive number such that t=1.+srelpr
c                  satisfies t.eq.1.  the quantity used
c                  here is internally restricted to be at
c                  least srelpr.  the data set for this option
c                  is the new tolerance.
c
c                  key=5
c                         change the rank determination tolerance for
c                  the reduced least squares equations from
c                  the nominal value of sqrt(srelpr).  this quantity can
c                  be no smaller than srelpr, the arithmetic-
c                  storage precision.  the quantity used
c                  here is internally restricted to be at
c                  least srelpr.  the data set for this option
c                  is the new tolerance.
c
c                  for example, suppose we want to change
c                  the tolerance for the reduced least squares
c                  problem, compute the covariance matrix of
c                  the solution parameters, and provide
c                  column scaling for the data matrix.  for
c                  these options the dimension of prgopt(*)
c                  must be at least n+9.  the fortran statements
c                  defining these options would be as follows.
c
c                  prgopt(1)=4 (link to entry 4 in prgopt(*))
c                  prgopt(2)=1 (covariance matrix key)
c                  prgopt(3)=1 (covariance matrix wanted)
c
c                  prgopt(4)=7 (link to entry 7 in prgopt(*))
c                  prgopt(5)=5 (least squares equas. tolerance key)
c                  prgopt(6)=... (new value of the tolerance)
c
c                  prgopt(7)=n+9 (link to entry n+9 in prgopt(*))
c                  prgopt(8)=3 (user-provided column scaling key)
c
c                  call scopy(n,d,1,prgopt(9),1)  (copy the n
c                    scaling factors from the user array d(*)
c                    to prgopt(9)-prgopt(n+8))
c
c                  prgopt(n+9)=1 (no more options to change)
c
c                  the contents of prgopt(*) are not modified
c                  by the subprogram.
c                  the options for wnnls( ) can also be included
c                  in this array.  the values of key recognized
c                  by wnnls( ) are 6, 7 and 8.  their functions
c                  are documented in the usage instructions for
c                  subroutine wnnls( ).  normally these options
c                  do not need to be modified when using lsei( ).
c
c     ip(1),       the amounts of working storage actually
c     ip(2)        allocated for the working arrays ws(*) and
c                  ip(*), respectively.  these quantities are
c                  compared with the actual amounts of storage
c                  needed by lsei( ).  insufficient storage
c                  allocated for either ws(*) or ip(*) is an
c                  error.  this feature was included in lsei( )
c                  because miscalculating the storage formulas
c                  for ws(*) and ip(*) might very well lead to
c                  subtle and hard-to-find execution errors.
c
c                  the length of ws(*) must be at least
c
c                  lw = 2*(me+n)+k+(mg+2)*(n+7)
c
c
c                  where k = max(ma+mg,n)
c                  this test will not be made if ip(1).le.0.
c
c                  the length of ip(*) must be at least
c
c                  lip = mg+2*n+2
c                  this test will not be made if ip(2).le.0.
c
c     output..
c
c     x(*),rnorme,  the array x(*) contains the solution parameters
c     rnorml        if the integer output flag mode = 0 or 1.
c                   the definition of mode is given directly below.
c                   when mode = 0 or 1, rnorme and rnorml
c                   respectively contain the residual vector
c                   euclidean lengths of f - ex and b - ax.  when
c                   mode=1 the equality constraint equations ex=f
c                   are contradictory, so rnorme.ne.0. the residual
c                   vector f-ex has minimal euclidean length. for
c                   mode.ge.2, none of these parameters are
c                   defined.
c
c     mode          integer flag that indicates the subprogram
c                   status after completion.  if mode.ge.2, no
c                   solution has been computed.
c
c                   mode =
c
c                   0  both equality and inequality constraints
c                      are compatible and have been satisfied.
c
c                   1  equality constraints are contradictory.
c                      a generalized inverse solution of ex=f was used
c                      to minimize the residual vector length f-ex.
c                      in this sense, the solution is still meaningful.
c
c                   2  inequality constraints are contradictory.
c
c                   3  both equality and inequality constraints
c                      are contradictory.
c
c                   the following interpretation of
c                   mode=1,2 or 3 must be made.  the
c                   sets consisting of all solutions
c                   of the equality constraints ex=f
c                   and all vectors satisfying gx.ge.h
c                   have no points on common.  (in
c                   particular this does not say that
c                   each individual set has no points
c                   at all, although this could be the
c                   case.)
c
c                   4  usage error occurred.  the value
c                      of mdw is .lt. me+ma+mg, mdw is
c                      .lt. n and a covariance matrix is
c                      requested, the option vector
c                      prgopt(*) is not properly defined,
c                      or the lengths of the working arrays
c                      ws(*) and ip(*), when specified in
c                      ip(1) and ip(2) respectively, are not
c                      long enough.
c
c     w(*,*)        the array w(*,*) contains the n by n symmetric
c                   covariance matrix of the solution parameters,
c                   provided this was requested on input with
c                   the option vector prgopt(*) and the output
c                   flag is returned with mode = 0 or 1.
c
c     ip(*)         the integer working array has three entries
c                   that provide rank and working array length
c                   information after completion.
c
c                      ip(1) = rank of equality constraint
c                              matrix.  define this quantity
c                              as kranke.
c
c                      ip(2) = rank of reduced least squares
c                              problem.
c
c                      ip(3) = the amount of storage in the
c                              working array ws(*) that was
c                              actually used by the subprogram.
c                              the formula given above for the length
c                              of ws(*) is a necessary overestimate.
c     user designated
c     working arrays..
c
c     ws(*),ip(*)              these are resp. type floating point
c                              and type integer working arrays.
c                              their required minimal lengths are
c                              given above.
c
c
c     subroutines called
c
c     lsi           part of this package.  solves a
c                   constrained least squares problem with
c                   inequality constraints.
c
c++
c     sdot,sscal,   subroutines from the blas package.
c     saxpy,sasum,  see trans. math. soft., vol. 5, no. 3, p. 308.
c     scopy,snrm2,
c     sswap
c
c     h12           subroutine to construct and apply a
c                   householder transformation.
c
c     xerror        from slatec error processing package.
c                   this is documented in sandia tech. rept.,
c                   sand78-1189.
c
c     subroutine lsei(w,mdw,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,ws,
c    1 ip)
c
c     revised oct. 1, 1981.
c

C 20060607 tbt CHange (1) with (*)
C     real             w(mdw,1), prgopt(1), x(1), ws(1), rnorme, rnorml

      real             w(mdw,*), prgopt(*), x(*), ws(*), rnorme, rnorml
      integer ip(3)
      real             dummy, enorm, srelpr, fnorm, gam, half, one, rb,
     *rn, rnmax, size, sn, snmax, t, tau, uj, up, vj, xnorm, xnrme, zero
      real             sasum, sdot, snrm2, sqrt, abs, amax1
      logical cov
c     data zero /0.e0/, srelpr /1.0e-12/, half /0.5e0/, one /1.e0/
      data zero /0.e0/, srelpr /0.e0/, half /0.5e0/, one /1.e0/
c
c     check that enough storage was allocated in ws(*) and ip(*).
      if (.not.(ip(1).gt.0)) go to 20
      lchk = 2*(me+n) + max0(ma+mg,n) + (mg+2)*(n+7)
      if (.not.(ip(1).lt.lchk)) go to 10
      mode = 4
      nerr = 2
      iopt = 1
      call xerrwv(67hlsei( ), insufficient storage allocated for ws(*),
     *need lw=i1 below, 67, nerr, iopt, 1, lchk, 0,
     * 0, dummy, dummy)
      return
   10 continue
   20 if (.not.(ip(2).gt.0)) go to 40
      lchk = mg + 2*n + 2
      if (.not.(ip(2).lt.lchk)) go to 30
      mode = 4
      nerr = 2
      iopt = 1
      call xerrwv(68hlsei( ), insufficient storage allocated for ip(*),
     *need lip=i1 below, 68, nerr, iopt, 1, lchk, 0,
     * 0, dummy, dummy)
      return
   30 continue
c
c     compute machine precision=srelpr only when necessary.
   40 if (.not.(srelpr.eq.zero)) go to 70
      srelpr = one
   50 if (one+srelpr.eq.one) go to 60
      srelpr = srelpr*half
      go to 50
   60 srelpr = srelpr + srelpr
c
c     compute number of possible right multiplying householder
c     transformations.
   70 m = me + ma + mg
      mode = 0
      if (n.le.0 .or. m.le.0) return
      if (.not.(mdw.lt.m)) go to 80
      nerr = 1
      iopt = 1
      call xerror(36hlsei( ), mdw.lt.me+ma+mg is an error, 36, nerr,
     * iopt)
      mode = 4
      return
   80 np1 = n + 1
      kranke = min0(me,n)
      n1 = 2*kranke + 1
      n2 = n1 + n
c
c     process-option-vector
      assign 90 to igo990
      go to 480
   90 if (.not.(cov .and. mdw.lt.n)) go to 100
      nerr = 2
      iopt = 1
      call xerror(
     * 54hlsei( ), mdw.lt.n, when cov matrix needed, is an error, 54,
     * nerr, iopt)
      mode = 4
      return
  100 l = kranke
c
c     compute norm of equality constraint matrix and rt side.
      enorm = zero
      do 110 j=1,n
        enorm = amax1(enorm,sasum(me,w(1,j),1))
  110 continue
      fnorm = sasum(me,w(1,np1),1)
      if (.not.(l.gt.0)) go to 190
      snmax = zero
      rnmax = zero
      do 180 i=1,l
c
c     compute maximum ratio of vector lengths. partition
c     is at col. i.
        do 150 k=i,me
          sn = sdot(n-i+1,w(k,i),mdw,w(k,i),mdw)
          rn = sdot(i-1,w(k,1),mdw,w(k,1),mdw)
          if (.not.(rn.eq.zero .and. sn.gt.snmax)) go to 120
          snmax = sn
          imax = k
          go to 140
  120     if (.not.(k.eq.i .or. (sn*rnmax.gt.rn*snmax))) go to 130
          snmax = sn
          rnmax = rn
          imax = k
  130     continue
  140     continue
  150   continue
c
c     interchange rows if necessary.
        if (i.ne.imax) call sswap(np1, w(i,1), mdw, w(imax,1), mdw)
        if (.not.(snmax.gt.tau**2*rnmax)) go to 160
c
c     eliminate elems i+1,...,n in row i.
        call h12(1, i, i+1, n, w(i,1), mdw, ws(i), w(i+1,1), mdw, 1,
     *   m-i)
        go to 170
  160   kranke = i - 1
        go to 200
  170   continue
  180 continue
  190 continue
  200 continue
c
c     save diag. terms of lower trap. matrix.
      call scopy(kranke, w, mdw+1, ws(kranke+1), 1)
c
c     use householder trans from left to achieve kranke by kranke upper
c     triangular form.
      if (.not.(kranke.gt.0 .and. kranke.lt.me)) go to 220
      do 210 kk=1,kranke
        k = kranke + 1 - kk
c
c     apply tranformation to matrix cols. 1,...,k-1.
        call h12(1, k, kranke+1, me, w(1,k), 1, up, w, 1, mdw, k-1)
c
c     apply to rt side vector.
        call h12(2, k, kranke+1, me, w(1,k), 1, up, w(1,np1), 1, 1, 1)
  210 continue
  220 if (.not.(kranke.gt.0)) go to 240
c
c     solve for variables 1,...,kranke in new coordinates.
      call scopy(kranke, w(1,np1), 1, x, 1)
      do 230 i=1,kranke
        x(i) = (x(i)-sdot(i-1,w(i,1),mdw,x,1))/w(i,i)
  230 continue
c
c     compute residuals for reduced problem.
  240 mep1 = me + 1
      rnorml = zero
      if (.not.(me.lt.m)) go to 270
      do 260 i=mep1,m
        w(i,np1) = w(i,np1) - sdot(kranke,w(i,1),mdw,x,1)
        sn = sdot(kranke,w(i,1),mdw,w(i,1),mdw)
        rn = sdot(n-kranke,w(i,kranke+1),mdw,w(i,kranke+1),mdw)
        if (.not.(rn.le.tau**2*sn .and. kranke.lt.n)) go to 250
        w(i,kranke+1) = zero
        call scopy(n-kranke, w(i,kranke+1), 0, w(i,kranke+1), mdw)
  250   continue
  260 continue
c
c     compute equal. constraint equas. residual length.
  270 rnorme = snrm2(me-kranke,w(kranke+1,np1),1)
c
c     move reduced problem data upward if kranke.lt.me.
      if (.not.(kranke.lt.me)) go to 290
      do 280 j=1,np1
        call scopy(m-me, w(me+1,j), 1, w(kranke+1,j), 1)
  280 continue
c
c     compute soln of reduced problem.
  290 call lsi(w(kranke+1,kranke+1), mdw, ma, mg, n-kranke, prgopt,
     * x(kranke+1), rnorml, mode, ws(n2), ip(2))
      if (.not.(me.gt.0)) go to 330
c
c     test for consistency of equality constraints.
      mdeqc = 0
      xnrme = sasum(kranke,w(1,np1),1)
      if (rnorme.gt.tau*(enorm*xnrme+fnorm)) mdeqc = 1
      mode = mode + mdeqc
c
c     check if soln to equal. constraints satisfies inequal.
c     constraints when there are no degrees of freedom left.
      if (.not.(kranke.eq.n .and. mg.gt.0)) go to 320
      xnorm = sasum(n,x,1)
      mapke1 = ma + kranke + 1
      mend = ma + kranke + mg
      do 310 i=mapke1,mend
        size = sasum(n,w(i,1),mdw)*xnorm + abs(w(i,np1))
        if (.not.(w(i,np1).gt.tau*size)) go to 300
        mode = mode + 2
        go to 450
  300   continue
  310 continue
  320 continue
  330 if (.not.(kranke.gt.0)) go to 420
c
c     replace diag. terms of lower trap. matrix.
      call scopy(kranke, ws(kranke+1), 1, w, mdw+1)
c
c     reapply trans to put soln in original coordinates.
      do 340 ii=1,kranke
        i = kranke + 1 - ii
        call h12(2, i, i+1, n, w(i,1), mdw, ws(i), x, 1, 1, 1)
  340 continue
c
c     compute cov matrix of equal. constrained problem.
      if (.not.(cov)) go to 410
      do 400 jj=1,kranke
        j = kranke + 1 - jj
        if (.not.(j.lt.n)) go to 390
        rb = ws(j)*w(j,j)
        if (rb.ne.zero) rb = one/rb
        jp1 = j + 1
        do 350 i=jp1,n
          w(i,j) = sdot(n-j,w(i,jp1),mdw,w(j,jp1),mdw)*rb
  350   continue
        gam = sdot(n-j,w(jp1,j),1,w(j,jp1),mdw)*rb
        gam = half*gam
        call saxpy(n-j, gam, w(j,jp1), mdw, w(jp1,j), 1)
        do 370 i=jp1,n
          do 360 k=i,n
            w(i,k) = w(i,k) + w(j,i)*w(k,j) + w(i,j)*w(j,k)
            w(k,i) = w(i,k)
  360     continue
  370   continue
        uj = ws(j)
        vj = gam*uj
        w(j,j) = uj*vj + uj*vj
        do 380 i=jp1,n
          w(j,i) = uj*w(i,j) + vj*w(j,i)
  380   continue
        call scopy(n-j, w(j,jp1), mdw, w(jp1,j), 1)
  390   continue
  400 continue
  410 continue
c
c     apply the scaling to the covariance matrix.
  420 if (.not.(cov)) go to 440
      do 430 i=1,n
        l = n1 + i
        call sscal(n, ws(l-1), w(i,1), mdw)
        call sscal(n, ws(l-1), w(1,i), 1)
  430 continue
  440 continue
  450 continue
c
c     rescale soln. vector.
      if (.not.(mode.le.1)) go to 470
      do 460 j=1,n
        l = n1 + j
        x(j) = x(j)*ws(l-1)
  460 continue
  470 ip(1) = kranke
      ip(3) = ip(3) + 2*kranke + n
      return
  480 continue
c     to process-option-vector
c
c     the nominal tolerance used in the code
c     for the equality constraint equations.
      tau = sqrt(srelpr)
c
c     the nominal column scaling used in the code is
c     the identity scaling.
      ws(n1) = one
      call scopy(n, ws(n1), 0, ws(n1), 1)
c
c     no covariance matrix is nominally computed.
      cov = .false.
c
c     define bound for number of options to change.
      nopt = 1000
      ntimes = 0
c
c     define bound for positive values of link.
      nlink = 100000
      last = 1
      link = prgopt(1)
      if (.not.(link.le.0 .or. link.gt.nlink)) go to 490
      nerr = 3
      iopt = 1
      call xerror(38hlsei( ) the option vector is undefined, 38, nerr,
     * iopt)
      mode = 4
      return
  490 if (.not.(link.gt.1)) go to 540
      ntimes = ntimes + 1
      if (.not.(ntimes.gt.nopt)) go to 500
      nerr = 3
      iopt = 1
      call xerror(
     * 52hlsei( ). the links in the option vector are cycling., 52,
     * nerr, iopt)
      mode = 4
      return
  500 key = prgopt(last+1)
      if (key.eq.1) cov = prgopt(last+2).ne.zero
      if (.not.(key.eq.2 .and. prgopt(last+2).ne.zero)) go to 520
      do 510 j=1,n
        t = snrm2(m,w(1,j),1)
        if (t.ne.zero) t = one/t
        l = n1 + j
        ws(l-1) = t
  510 continue
  520 if (key.eq.3) call scopy(n, prgopt(last+2), 1, ws(n1), 1)
      if (key.eq.4) tau = amax1(srelpr,prgopt(last+2))
      next = prgopt(link)
      if (.not.(next.le.0 .or. next.gt.nlink)) go to 530
      nerr = 3
      iopt = 1
      call xerror(38hlsei( ) the option vector is undefined, 38, nerr,
     * iopt)
      mode = 4
      return
  530 last = link
      link = next
      go to 490
  540 do 550 j=1,n
        l = n1 + j
        call sscal(m, ws(l-1), w(1,j), 1)
  550 continue
      go to 560
  560 go to igo990, (90)
      end
      subroutine lsi(w, mdw, ma, mg, n, prgopt, x, rnorm, mode, ws, ip)
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start editing at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/sasum/dasum/,/sdot/ddot/,
c     / sqrt/ dsqrt/,/amax1/dmax1/,/sswap/dswap/,
c     /scopy/dcopy/,/sscal/dscal/,/saxpy/daxpy/,/e0/d0/,/srelpr/drelpr/
c
c     this is a companion subprogram to lsei( ).
c     the documentation for lsei( ) has more complete
c     usage instructions.
c     written by r. j. hanson, sla.
c
c     solve..
c              ax = b,  a  ma by n  (least squares equations)
c     subject to..
c
c              gx.ge.h, g  mg by n  (inequality constraints)
c
c     input..
c
c      w(*,*) contains  (a b) in rows 1,...,ma+mg, cols 1,...,n+1.
c                       (g h)
c
c     mdw,ma,mg,n
c              contain (resp) var. dimension of w(*,*),
c              and matrix dimensions.
c
c     prgopt(*),
c              program option vector.
c
c     output..
c
c      x(*),rnorm
c
c              solution vector(unless mode=2), length of ax-b.
c
c      mode
c              =0   inequality constraints are compatible.
c              =2   inequality constraints contradictory.
c
c      ws(*),
c              working storage of dimension k+n+(mg+2)*(n+7),
c              where k=max(ma+mg,n).
c      ip(mg+2*n+1)
c              integer working storage
c      revised oct. 1, 1981.
c
c     subroutines called
c
c     lpdp          this subprogram minimizes a sum of squares
c                   of unknowns subject to linear inequality
c                   constraints.  part of this package.
c
c++
c     sdot,sscal    subroutines from the blas package.
c     saxpy,sasum,  see trans. math. soft., vol. 5, no. 3, p. 308.
c     scopy,sswap
c
c     hfti          solves an unconstrained linear least squares
c                   problem.  part of this package.
c
c     h12           subroutine to construct and apply a householder
c                   transformation.
c
c     subroutine lsi(w,mdw,ma,mg,n,prgopt,x,rnorm,mode,ws,ip)
c

C 20060607 tbt - Changed (1) to (*)
C     real             w(mdw,1), prgopt(1), rnorm, ws(1), x(1)
C     integer ip(1)

      real             w(mdw,*), prgopt(*), rnorm, ws(*), x(*)
      integer ip(*)

      real             anorm, srelpr, fac, gam, half, one, rb, tau, tol,
     * xnorm, zero
      real             sasum, sdot, sqrt, amax1
      logical cov
c
      data zero /0.e0/, srelpr /0.e0/, one /1.e0/, half /.5e0/
c     data zero /0.e0/, srelpr /1.0e-12/, one /1.e0/, half /.5e0/
c
C------------------------------------------------------------

c     compute machine precision=srelpr only when necessary.
      if (.not.(srelpr.eq.zero)) go to 30
      srelpr = one
   10 if (one+srelpr.eq.one) go to 20
      srelpr = srelpr*half
      go to 10
   20 srelpr = srelpr + srelpr
   30 mode = 0
      rnorm = zero
      m = ma + mg
      np1 = n + 1
      krank = 0
      if (n.le.0 .or. m.le.0) go to 70
      assign 40 to igo994
      go to 500
c
c     process-option-vector
c
c     compute matrix norm of least squares equas.
   40 anorm = zero
      do 50 j=1,n
        anorm = amax1(anorm,sasum(ma,w(1,j),1))
   50 continue
c
c     set tol for hfti( ) rank test.
      tau = tol*anorm
c
c     compute householder orthogonal decomp of matrix.
      if (n.gt.0) ws(1) = zero
      call scopy(n, ws, 0, ws, 1)
      call scopy(ma, w(1,np1), 1, ws, 1)
      k = max0(m,n)
      minman = min0(ma,n)
      n1 = k + 1
      n2 = n1 + n
      call hfti(w, mdw, ma, n, ws, 1, 1, tau, krank, rnorm, ws(n2),
     * ws(n1), ip)
      fac = one
      gam=ma-krank
      if (krank.lt.ma) fac = rnorm**2/gam
      assign 60 to igo990
      go to 80
c
c     reduce-to-lpdp-and-solve
   60 continue
   70 ip(1) = krank
      ip(2) = n + max0(m,n) + (mg+2)*(n+7)
      return
   80 continue
c
c     to reduce-to-lpdp-and-solve
      map1 = ma + 1
c
c     compute ineq. rt-hand side for lpdp.
      if (.not.(ma.lt.m)) go to 260
      if (.not.(minman.gt.0)) go to 160
      do 90 i=map1,m
        w(i,np1) = w(i,np1) - sdot(n,w(i,1),mdw,ws,1)
   90 continue
      do 100 i=1,minman
        j = ip(i)
c
c     apply permutations to cols of ineq. constraint matrix.
        call sswap(mg, w(map1,i), 1, w(map1,j), 1)
  100 continue
c
c     apply householder transformations to constraint matrix.
      if (.not.(0.lt.krank .and. krank.lt.n)) go to 120
      do 110 ii=1,krank
        i = krank + 1 - ii
        l = n1 + i
        call h12(2, i, krank+1, n, w(i,1), mdw, ws(l-1), w(map1,1),
     *   mdw, 1, mg)
  110 continue
c
c     compute permuted ineq. constr. matrix times r-inverse.
  120 do 150 i=map1,m
        if (.not.(0.lt.krank)) go to 140
        do 130 j=1,krank
          w(i,j) = (w(i,j)-sdot(j-1,w(1,j),1,w(i,1),mdw))/w(j,j)
  130   continue
  140   continue
  150 continue
c
c     solve the reduced problem with lpdp algorithm,
c     the least projected distance problem.
  160 call lpdp(w(map1,1), mdw, mg, krank, n-krank, prgopt, x, xnorm,
     * mdlpdp, ws(n2), ip(n+1))
      if (.not.(mdlpdp.eq.1)) go to 240
      if (.not.(krank.gt.0)) go to 180
c
c     compute soln in original coordinates.
      do 170 ii=1,krank
        i = krank + 1 - ii
        x(i) = (x(i)-sdot(ii-1,w(i,i+1),mdw,x(i+1),1))/w(i,i)
  170 continue
c
c     apply householder trans. to soln vector.
  180 if (.not.(0.lt.krank .and. krank.lt.n)) go to 200
      do 190 i=1,krank
        l = n1 + i
        call h12(2, i, krank+1, n, w(i,1), mdw, ws(l-1), x, 1, 1, 1)
  190 continue
  200 if (.not.(minman.gt.0)) go to 230
c
c     repermute variables to their input order.
      do 210 ii=1,minman
        i = minman + 1 - ii
        j = ip(i)
        call sswap(1, x(i), 1, x(j), 1)
  210 continue
c
c     variables are now in orig. coordinates.
c     add soln of unsconstrained prob.
      do 220 i=1,n
        x(i) = x(i) + ws(i)
  220 continue
c
c     compute the residual vector norm.
      rnorm = sqrt(rnorm**2+xnorm**2)
  230 go to 250
  240 mode = 2
  250 go to 270
  260 call scopy(n, ws, 1, x, 1)
  270 if (.not.(cov .and. krank.gt.0)) go to 490
c
c     compute covariance matrix based on the orthogonal decomp.
c     from hfti( ).
c
      krm1 = krank - 1
      krp1 = krank + 1
c
c     copy diag. terms to working array.
      call scopy(krank, w, mdw+1, ws(n2), 1)
c
c     reciprocate diag. terms.
      do 280 j=1,krank
        w(j,j) = one/w(j,j)
  280 continue
      if (.not.(krank.gt.1)) go to 310
c
c     invert the upper triangular qr factor on itself.
      do 300 i=1,krm1
        ip1 = i + 1
        do 290 j=ip1,krank
          w(i,j) = -sdot(j-i,w(i,i),mdw,w(i,j),1)*w(j,j)
  290   continue
  300 continue
c
c     compute the inverted factor times its transpose.
  310 do 330 i=1,krank
        do 320 j=i,krank
          w(i,j) = sdot(krank+1-j,w(i,j),mdw,w(j,j),mdw)
  320   continue
  330 continue
      if (.not.(krank.lt.n)) go to 450
c
c     zero out lower trapezoidal part.
c     copy upper tri. to lower tri. part.
      do 340 j=1,krank
        call scopy(j, w(1,j), 1, w(j,1), mdw)
  340 continue
      do 350 i=krp1,n
        w(i,1) = zero
        call scopy(i, w(i,1), 0, w(i,1), mdw)
  350 continue
c
c     apply right side transformations to lower tri.
      n3 = n2 + krp1
      do 430 i=1,krank
        l = n1 + i
        k = n2 + i
        rb = ws(l-1)*ws(k-1)
        if (.not.(rb.lt.zero)) go to 420
c
c     if rb.ge.zero, transformation can be regarded as zero.
        rb = one/rb
c
c     store unscaled rank-one householder update in work array.
        ws(n3) = zero
        call scopy(n, ws(n3), 0, ws(n3), 1)
        l = n1 + i
        k = n3 + i
        ws(k-1) = ws(l-1)
        do 360 j=krp1,n
          k = n3 + j
          ws(k-1) = w(i,j)
  360   continue
        do 370 j=1,n
          l = n3 + i
          k = n3 + j
          ws(j) = sdot(j-i,w(j,i),mdw,ws(l-1),1) + sdot(n-j+1,w(j,j),1,
     *     ws(k-1),1)
          ws(j) = ws(j)*rb
  370   continue
        l = n3 + i
        gam = sdot(n-i+1,ws(l-1),1,ws(i),1)*rb
        gam = gam*half
        call saxpy(n-i+1, gam, ws(l-1), 1, ws(i), 1)
        do 410 j=i,n
          if (.not.(i.gt.1)) go to 390
          im1 = i - 1
          k = n3 + j
          do 380 l=1,im1
            w(j,l) = w(j,l) + ws(k-1)*ws(l)
  380     continue
  390     k = n3 + j
          do 400 l=i,j
            il = n3 + l
            w(j,l) = w(j,l) + ws(j)*ws(il-1) + ws(l)*ws(k-1)
  400     continue
  410   continue
  420   continue
  430 continue
c
c     copy lower tri. to upper tri. to symmetrize the covariance matrix.
      do 440 i=1,n
        call scopy(i, w(i,1), mdw, w(1,i), 1)
  440 continue
c
c     repermute rows and cols.
  450 do 470 ii=1,minman
        i = minman + 1 - ii
        k = ip(i)
        if (.not.(i.ne.k)) go to 460
        call sswap(1, w(i,i), 1, w(k,k), 1)
        call sswap(i-1, w(1,i), 1, w(1,k), 1)
        call sswap(k-i-1, w(i,i+1), mdw, w(i+1,k), 1)
        call sswap(n-k, w(i,k+1), mdw, w(k,k+1), mdw)
  460   continue
  470 continue
c
c     put in normalized residual sum of squares scale factor
c     and symmetrize the resulting covariance marix.
      do 480 j=1,n
        call sscal(j, fac, w(1,j), 1)
        call scopy(j, w(1,j), 1, w(j,1), mdw)
  480 continue
  490 go to 540
  500 continue
c
c     to process-option-vector
c
c     the nominal tolerance used in the code,
      tol = sqrt(srelpr)
      cov = .false.
      last = 1
      link = prgopt(1)
  510 if (.not.(link.gt.1)) go to 520
      key = prgopt(last+1)
      if (key.eq.1) cov = prgopt(last+2).ne.zero
      if (key.eq.5) tol = amax1(srelpr,prgopt(last+2))
      next = prgopt(link)
      last = link
      link = next
      go to 510
  520 go to 530
  530 go to igo994, (40)
  540 go to igo990, (60)
      end

      subroutine lpdp(a, mda, m, n1, n2, prgopt, x, wnorm, mode, ws, is)

c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start editing at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/snrm2/dnrm2/,/sdot/ddot/,
c     /scopy/dcopy/,/sscal/dscal/,/abs(/dabs(/, abs/, dabs/,/e0/d0/
c
c     dimension a(mda,n+1),prgopt(*),x(n),ws((m+2)*(n+7)),is(m+n+1),
c     where n=n1+n2.  this is a slight overestimate for ws(*).
c
c     written by r. j. hanson and k. h. haskell, sandia labs
c     revised oct. 1, 1981.
c
c     determine an n1-vector w, and
c               an n2-vector z
c     which minimizes the euclidean length of w
c     subject to g*w+h*z .ge. y.
c     this is the least projected distance problem, lpdp.
c     the matrices g and h are of respective
c     dimensions m by n1 and m by n2.
c
c     called by subprogram lsi( ).
c
c     the matrix
c                (g h y)
c
c     occupies rows 1,...,m and cols 1,...,n1+n2+1 of a(*,*).
c
c     the solution (w) is returned in x(*).
c                  (z)
c
c     the value of mode indicates the status of
c     the computation after returning to the user.
c
c          mode=1  the solution was successfully obtained.
c
c          mode=2  the inequalities are inconsistent.
c
c     subroutines called
c
c     wnnls         solves a nonnegatively constrained linear least
c                   squares problem with linear equality constraints.
c                   part of this package.
c
c++
c     sdot,         subroutines from the blas package.
c     sscal,snrm2,  see trans. math. soft., vol. 5, no. 3, p. 308.
c     scopy

C----------------------------------------------------------------------------

C 20060607 tbt Changed to (*)
C     real             a(mda,1), prgopt(1), ws(1), wnorm, x(1)
C     integer is(1)
 
      real             a(mda,*), prgopt(*), ws(*), wnorm, x(*)
      integer is(*)
 
      real             fac, one, rnorm, sc, ynorm, zero
      real             sdot, snrm2, abs
      data zero, one /0.e0,1.e0/, fac /0.1e0/

C----------------------------------------------------------------------------

      n = n1 + n2
      mode = 1
      if (.not.(m.le.0)) go to 20
      if (.not.(n.gt.0)) go to 10
      x(1) = zero
      call scopy(n, x, 0, x, 1)
   10 wnorm = zero
      return
   20 np1 = n + 1
c
c     scale nonzero rows of inequality matrix to have length one.
      do 40 i=1,m
        sc = snrm2(n,a(i,1),mda)
        if (.not.(sc.ne.zero)) go to 30
        sc = one/sc
        call sscal(np1, sc, a(i,1), mda)
   30   continue
   40 continue
c
c     scale rt.-side vector to have length one (or zero).
      ynorm = snrm2(m,a(1,np1),1)
      if (.not.(ynorm.ne.zero)) go to 50
      sc = one/ynorm
      call sscal(m, sc, a(1,np1), 1)
c
c     scale cols of matrix h.
   50 j = n1 + 1
   60 if (.not.(j.le.n)) go to 70
      sc = snrm2(m,a(1,j),1)
      if (sc.ne.zero) sc = one/sc
      call sscal(m, sc, a(1,j), 1)
      x(j) = sc
      j = j + 1
      go to 60
   70 if (.not.(n1.gt.0)) go to 130
c
c     copy transpose of (h g y) to work array ws(*).
      iw = 0
      do 80 i=1,m
c
c     move col of transpose of h into work array.
        call scopy(n2, a(i,n1+1), mda, ws(iw+1), 1)
        iw = iw + n2
c
c     move col of transpose of g into work array.
        call scopy(n1, a(i,1), mda, ws(iw+1), 1)
        iw = iw + n1
c
c     move component of vector y into work array.
        ws(iw+1) = a(i,np1)
        iw = iw + 1
   80 continue
      ws(iw+1) = zero
      call scopy(n, ws(iw+1), 0, ws(iw+1), 1)
      iw = iw + n
      ws(iw+1) = one
      iw = iw + 1
c
c     solve eu=f subject to (transpose of h)u=0, u.ge.0.  the
c     matrix e = transpose of (g y), and the (n+1)-vector
c     f = transpose of (0,...,0,1).
      ix = iw + 1
      iw = iw + m
c
c     do not check lengths of work arrays in this usage of wnnls( ).
      is(1) = 0
      is(2) = 0
      call wnnls(ws, np1, n2, np1-n2, m, 0, prgopt, ws(ix), rnorm,
     * modew, is, ws(iw+1))
c
c     compute the components of the soln denoted above by w.
      sc = one - sdot(m,a(1,np1),1,ws(ix),1)
      if (.not.(one+fac*abs(sc).ne.one .and. rnorm.gt.zero)) go to 110
      sc = one/sc
      do 90 j=1,n1
        x(j) = sc*sdot(m,a(1,j),1,ws(ix),1)
   90 continue
c
c     compute the vector q=y-gw.  overwrite y with this vector.
      do 100 i=1,m
        a(i,np1) = a(i,np1) - sdot(n1,a(i,1),mda,x,1)
  100 continue
      go to 120
  110 mode = 2
      return
  120 continue
  130 if (.not.(n2.gt.0)) go to 180
c
c     copy transpose of (h q) to work array ws(*).
      iw = 0
      do 140 i=1,m
        call scopy(n2, a(i,n1+1), mda, ws(iw+1), 1)
        iw = iw + n2
        ws(iw+1) = a(i,np1)
        iw = iw + 1
  140 continue
      ws(iw+1) = zero
      call scopy(n2, ws(iw+1), 0, ws(iw+1), 1)
      iw = iw + n2
      ws(iw+1) = one
      iw = iw + 1
      ix = iw + 1
      iw = iw + m
c
c     solve rv=s subject to v.ge.0.  the matrix r =(transpose
c     of (h q)), where q=y-gw.  the (n2+1)-vector s =(transpose
c     of (0,...,0,1)).
c
c     do not check lengths of work arrays in this usage of wnnls( ).
      is(1) = 0
      is(2) = 0
      call wnnls(ws, n2+1, 0, n2+1, m, 0, prgopt, ws(ix), rnorm, modew,
     * is, ws(iw+1))
c
c     compute the components of the soln denoted above by z.
      sc = one - sdot(m,a(1,np1),1,ws(ix),1)
      if (.not.(one+fac*abs(sc).ne.one .and. rnorm.gt.zero)) go to 160
      sc = one/sc
      do 150 j=1,n2
        l = n1 + j
        x(l) = sc*sdot(m,a(1,l),1,ws(ix),1)*x(l)
  150 continue
      go to 170
  160 mode = 2
      return
  170 continue
c
c     account for scaling of rt.-side vector in solution.
  180 call sscal(n, ynorm, x, 1)
      wnorm = snrm2(n1,x,1)
      return
      end

C-------------------------------------------------------

      subroutine wnnls(w, mdw, me, ma, n, l, prgopt, x, rnorm, mode,
     * iwork, work)

c
c     dimension w(mdw,n+1),prgopt(*),x(n),iwork(m+n),work(m+5*n)
c
c     abstract
c
c     this subprogram solves a linearly constrained least squares
c     problem.  suppose there are given matrices e and a of
c     respective dimensions me by n and ma by n, and vectors f
c     and b of respective lengths me and ma.  this subroutine
c     solves the problem
c
c               ex = f, (equations to be exactly satisfied)
c
c               ax = b, (equations to be approximately satisfied,
c                        in the least squares sense)
c
c               subject to components l+1,...,n nonnegative
c
c     any values me.ge.0, ma.ge.0 and 0.le. l .le.n are permitted.
c
c     the problem is reposed as problem wnnls
c
c               (wt*e)x = (wt*f)
c               (   a)    (   b), (least squares)
c               subject to components l+1,...,n nonnegative.
c
c     the subprogram chooses the heavy weight (or penalty parameter) wt.
c
c     the parameters for wnnls are
c
c     input..
c
c     w(*,*),mdw,  the array w(*,*) is double subscripted with first
c     me,ma,n,l    dimensioning parameter equal to mdw.  for this
c                  discussion let us call m = me + ma.  then mdw
c                  must satisfy mdw.ge.m.  the condition mdw.lt.m
c                  is an error.
c
c                  the array w(*,*) contains the matrices and vectors
c
c                       (e  f)
c                       (a  b)
c
c                  in rows and columns 1,...,m and 1,...,n+1
c                  respectively.  columns 1,...,l correspond to
c                  unconstrained variables x(1),...,x(l).  the
c                  remaining variables are constrained to be
c                  nonnegative.  the condition l.lt.0 .or. l.gt.n is
c                  an error.
c
c     prgopt(*)    this array is the option vector.
c                  if the user is satisfied with the nominal
c                  subprogram features set
c
c                  prgopt(1)=1 (or prgopt(1)=1.0)
c
c                  otherwise prgopt(*) is a linked list consisting of
c                  groups of data of the following form
c
c                  link
c                  key
c                  data set
c
c                  the parameters link and key are each one word.
c                  the data set can be comprised of several words.
c                  the number of items depends on the value of key.
c                  the value of link points to the first
c                  entry of the next group of data within
c                  prgopt(*).  the exception is when there are
c                  no more options to change.  in that
c                  case link=1 and the values key and data set
c                  are not referenced. the general layout of
c                  prgopt(*) is as follows.
c
c               ...prgopt(1)=link1 (link to first entry of next group)
c               .  prgopt(2)=key1 (key to the option change)
c               .  prgopt(3)=data value (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...prgopt(link1)=link2 (link to the first entry of
c               .                       next group)
c               .  prgopt(link1+1)=key2 (key to the option change)
c               .  prgopt(link1+2)=data value
c               ...     .
c               .       .
c               .       .
c               ...prgopt(link)=1 (no more options to change)
c
c                  values of link that are nonpositive are errors.
c                  a value of link.gt.nlink=100000 is also an error.
c                  this helps prevent using invalid but positive
c                  values of link that will probably extend
c                  beyond the program limits of prgopt(*).
c                  unrecognized values of key are ignored.  the
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  to prevent cycling in the
c                  processing of the option array a count of the
c                  number of options changed is maintained.
c                  whenever this count exceeds nopt=1000 an error
c                  message is printed and the subprogram returns.
c
c                  options..
c
c                  key=6
c                         scale the nonzero columns of the
c                  entire data matrix
c                  (e)
c                  (a)
c                  to have length one.  the data set for
c                  this option is a single value.  it must
c                  be nonzero if unit length column scaling is
c                  desired.
c
c                  key=7
c                         scale columns of the entire data matrix
c                  (e)
c                  (a)
c                  with a user-provided diagonal matrix.
c                  the data set for this option consists
c                  of the n diagonal scaling factors, one for
c                  each matrix column.
c
c                  key=8
c                         change the rank determination tolerance from
c                  the nominal value of sqrt(eps).  this quantity can
c                  be no smaller than eps, the arithmetic-
c                  storage precision.  the quantity used
c                  here is internally restricted to be at
c                  least eps.  the data set for this option
c                  is the new tolerance.
c
c                  key=9
c                         change the blow-up parameter from the
c                  nominal value of sqrt(eps).  the reciprocal of
c                  this parameter is used in rejecting solution
c                  components as too large when a variable is
c                  first brought into the active set.  too large
c                  means that the proposed component times the
c                  reciprocal of the parameteris not less than
c                  the ratio of the norms of the right-side
c                  vector and the data matrix.
c                  this parameter can be no smaller than eps,
c                  the arithmetic-storage precision.
c
c                  for example, suppose we want to provide
c                  a diagonal matrix to scale the problem
c                  matrix and change the tolerance used for
c                  determining linear dependence of dropped col
c                  vectors.  for these options the dimensions of
c                  prgopt(*) must be at least n+6.  the fortran
c                  statements defining these options would
c                  be as follows.
c
c                  prgopt(1)=n+3 (link to entry n+3 in prgopt(*))
c                  prgopt(2)=7 (user-provided scaling key)
c
c                  call scopy(n,d,1,prgopt(3),1) (copy the n
c                  scaling factors from a user array called d(*)
c                  into prgopt(3)-prgopt(n+2))
c
c                  prgopt(n+3)=n+6 (link to entry n+6 of prgopt(*))
c                  prgopt(n+4)=8 (linear dependence tolerance key)
c                  prgopt(n+5)=... (new value of the tolerance)
c
c                  prgopt(n+6)=1 (no more options to change)
c
c     iwork(1),    the amounts of working storage actually allocated
c     iwork(2)     for the working arrays work(*) and iwork(*),
c                  respectively.  these quantities are compared with
c                  the actual amounts of storage needed for wnnls( ).
c                  insufficient storage allocated for either work(*)
c                  or iwork(*) is considered an error.  this feature
c                  was included in wnnls( ) because miscalculating
c                  the storage formulas for work(*) and iwork(*)
c                  might very well lead to subtle and hard-to-find
c                  execution errors.
c
c                  the length of work(*) must be at least
c
c                  lw = me+ma+5*n
c                  this test will not be made if iwork(1).le.0.
c
c                  the length of iwork(*) must be at least
c
c                  liw = me+ma+n
c                  this test will not be made if iwork(2).le.0.
c
c     output..
c
c     x(*)         an array dimensioned at least n, which will
c                  contain the n components of the solution vector
c                  on output.
c
c     rnorm        the residual norm of the solution.  the value of
c                  rnorm contains the residual vector length of the
c                  equality constraints and least squares equations.
c
c     mode         the value of mode indicates the success or failure
c                  of the subprogram.
c
c                  mode = 0  subprogram completed successfully.
c
c                       = 1  max. number of iterations (equal to
c                            3*(n-l)) exceeded. nearly all problems
c                            should complete in fewer than this
c                            number of iterations. an approximate
c                            solution and its corresponding residual
c                            vector length are in x(*) and rnorm.
c
c                       = 2  usage error occurred.  the offending
c                            condition is noted with the error
c                            processing subprogram, xerror( ).
c
c     user-designated
c     working arrays..
c
c     work(*)      a working array of length at least
c                  m + 5*n.
c
c     iwork(*)     an integer-valued working array of length at least
c                  m+n.
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/, dummy/,sngl(dummy)/
c
c     written by karen h. haskell, sandia laboratories,
c     and r.j. hanson, sandia laboratories.
c     revised feb.25, 1982.
c
c     subroutines called by wnnls( )
c
c++
c     wnlsm         companion subroutine to wnnls( ), where
c                   most of the computation takes place.
c
c     xerror,xerrwv from slatec error processing package.
c                   this is documented in sandia tech. rept.,
c                   sand78-1189.
c
c     references
c
c     1. solving least squares problems, by c.l. lawson
c        and r.j. hanson.  prentice-hall, inc. (1974).
c
c     2. basic linear algebra subprograms for fortran usage, by
c        c.l. lawson, r.j. hanson, d.r. kincaid, and f.t. krogh.
c        toms, v. 5, no. 3, p. 308.  also available as
c        sandia technical report no. sand77-0898.
c
c     3. an algorithm for linear least squares with equality
c        and nonnegativity constraints, by k.h. haskell and
c        r.j. hanson.  available as sandia technical report no.
c        sand77-0552, and math. programming, vol. 21, (1981), p. 98-118.
c
c     4. slatec common math. library error handling
c        package.  by r. e. jones.  available as sandia
c        technical report sand78-1189.
c

C 20060607 tbt - Changed (1) to (*)
C     real              dummy, w(mdw,1), prgopt(1), x(1), work(1), rnorm
C     integer iwork(1)

      real              dummy, w(mdw,*), prgopt(*), x(*), work(*), rnorm
      integer iwork(*)

c
c
      mode = 0
      if (ma+me.le.0 .or. n.le.0) return
      if (.not.(iwork(1).gt.0)) go to 20
      lw = me + ma + 5*n
      if (.not.(iwork(1).lt.lw)) go to 10
      nerr = 2
      iopt = 1
      call xerrwv(70hwnnls( ), insufficient storage allocated for work(*
     *), need lw=i1 below, 70, nerr, iopt, 1, lw, 0, 0, dummy, dummy)
      mode = 2
      return
   10 continue
   20 if (.not.(iwork(2).gt.0)) go to 40
      liw = me + ma + n
      if (.not.(iwork(2).lt.liw)) go to 30
      nerr = 2
      iopt = 1
      call xerrwv(72hwnnls( ), insufficient storage allocated for iwork(
     **), need liw=i1 below, 72, nerr, iopt, 1, liw, 0, 0, dummy, dummy)
      mode = 2
      return
   30 continue
   40 if (.not.(mdw.lt.me+ma)) go to 50

      nerr = 1
      iopt = 1
      call xerror(44hwnnls( ), the value mdw.lt.me+ma is an error, 44,
     * nerr, iopt)
      mode = 2
      return
   50 if (0.le.l .and. l.le.n) go to 60
      nerr = 2
      iopt = 1
      call xerror(39hwnnls( ), l.le.0.and.l.le.n is required, 39, nerr,
     * iopt)
      mode = 2
      return
c
c     the purpose of this subroutine is to break up the arrays
c     work(*) and iwork(*) into separate work arrays
c     required by the main subroutine wnlsm( ).
c
   60 l1 = n + 1
      l2 = l1 + n
      l3 = l2 + me + ma
      l4 = l3 + n
      l5 = l4 + n
c
      call wnlsm(w, mdw, me, ma, n, l, prgopt, x, rnorm, mode, iwork,
     * iwork(l1), work(1), work(l1), work(l2), work(l3), work(l4),
     * work(l5))
      return
      end

C-----------------------------------------------

      subroutine wnlsm(w, mdw, mme, ma, n, l, prgopt, x, rnorm, mode,
     * ipivot, itype, wd, h, scale, z, temp, d)
c
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start changes at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/sasum/dasum/,/srotmg/drotmg/,
c     /snrm2/dnrm2/,/ sqrt/ dsqrt/,/srotm/drotm/,/amax1/dmax1/,
c     /scopy/dcopy/,/sscal/dscal/,/saxpy/daxpy/,/e0/d0/,/sswap/dswap/,
c     /isamax/idamax/,/srelpr/drelpr/
c
c     this is a companion subprogram to wnnls( ).
c     the documentation for wnnls( ) has more complete
c     usage instructions.
c
c     written by karen h. haskell, sandia laboratories,
c     with the help of r.j. hanson, sandia laboratories,
c     december 1976 - january 1978.
c     revised mar. 4, 1982.
c
c     in addition to the parameters discussed in the prologue to
c     subroutine wnnls, the following work arrays are used in
c     subroutine wnlsm  (they are passed through the calling
c     sequence from wnnls for purposes of variable dimensioning).
c     their contents will in general be of no interest to the user.
c
c         ipivot(*)
c            an array of length n.  upon completion it contains the
c         pivoting information for the cols of w(*,*).
c
c         itype(*)
c            an array of length m which is used to keep track
c         of the classification of the equations.  itype(i)=0
c         denotes equation i as an equality constraint.
c         itype(i)=1 denotes equation i as a least squares
c         equation.
c
c         wd(*)
c            an array of length n.  upon completion it contains the
c         dual solution vector.
c
c         h(*)
c            an array of length n.  upon completion it contains the
c         pivot scalars of the householder transformations performed
c         in the case krank.lt.l.
c
c         scale(*)
c            an array of length m which is used by the subroutine
c         to store the diagonal matrix of weights.
c         these are used to apply the modified givens
c         transformations.
c
c         z(*),temp(*)
c            working arrays of length n.
c
c         d(*)
c            an array of length n that contains the
c         column scaling for the matrix (e).
c                                       (a)
c
c     subroutine wnlsm (w,mdw,mme,ma,n,l,prgopt,x,rnorm,mode,
c    1                  ipivot,itype,wd,h,scale,z,temp,d)
c++

C 20060607 tbt (1) to (*)
c     real             w(mdw,1), x(1), wd(1), h(1), scale(1), dope(4)
c     real             z(1), temp(1), prgopt(1), d(1), sparam(5)

      real             w(mdw,*), x(*), wd(*), h(*), scale(*), dope(4)
      real             z(*), temp(*), prgopt(*), d(*), sparam(5)


      real             alamda, alpha, alsq, amax, bnorm, eanorm
      real             srelpr, fac, one, blowup
      real             rnorm, sm, t, tau, two, wmax, zero, zz, z2
      real             amax1, sqrt, snrm2, sasum

C     integer ipivot(1), itype(1), isamax, idope(8)
      integer ipivot(*), itype(*), isamax, idope(8)

      logical hitcon, feasbl, done, pos
c     data zero /0.e0/, one /1.e0/, two /2.e0/, srelpr /1.0e-12/
      data zero /0.e0/, one /1.e0/, two /2.e0/, srelpr /0.e0/
c
c     initialize-variables
      assign 10 to igo998
      go to 180
c
c     perform initial triangularization in the submatrix
c     corresponding to the unconstrained variables using
c     the procedure initially-triangularize.
   10 assign 20 to igo995
      go to 280
c
c     perform wnnls algorithm using the following steps.
c
c     until(done)
c
c        compute-search-direction-and-feasible-point
c
c        when (hitcon) add-constraints
c
c        else perform-multiplier-test-and-drop-a-constraint
c
c        fin
c
c     compute-final-solution
c
   20 if (done) go to 80
c
      assign 30 to igo991
      go to 300
c
c     compute-search-direction-and-feasible-point
c
   30 if (.not.(hitcon)) go to 50
      assign 40 to igo986
      go to 370
   40 go to 70
c
c     when (hitcon) add-constraints
c
   50 assign 60 to igo983
      go to 640
   60 continue
c
c     else perform-multiplier-test-and-drop-a-constraint
c
   70 go to 20
c
   80 assign 90 to igo980
      go to 1000
c
c     compute-final-solution
c
   90 return
  100 continue
c
c     to process-option-vector
      fac = 1.e-4
c
c     the nominal tolerance used in the code,
      tau = sqrt(srelpr)
c
c     the nominal blow-up factor used in the code.
      blowup = tau
c
c     the nominal column scaling used in the code is
c     the identity scaling.
      d(1) = one
      call scopy(n, d, 0, d, 1)
c
c     define bound for number of options to change.
      nopt = 1000
c
c     define bound for positive value of link.
      nlink = 100000
      ntimes = 0
      last = 1
      link = prgopt(1)
      if (.not.(link.le.0 .or. link.gt.nlink)) go to 110
      nerr = 3
      iopt = 1
      call xerror(39hwnnls( ) the option vector is undefined, 39, nerr,
     * iopt)
      mode = 2
      return
  110 if (.not.(link.gt.1)) go to 160
      ntimes = ntimes + 1
      if (.not.(ntimes.gt.nopt)) go to 120
      nerr = 3
      iopt = 1
      call xerror(
     * 53hwnnls( ). the links in the option vector are cycling., 53,
     * nerr, iopt)
      mode = 2
      return
  120 key = prgopt(last+1)
      if (.not.(key.eq.6 .and. prgopt(last+2).ne.zero)) go to 140
      do 130 j=1,n
        t = snrm2(m,w(1,j),1)
        if (t.ne.zero) t = one/t
        d(j) = t
  130 continue
  140 if (key.eq.7) call scopy(n, prgopt(last+2), 1, d, 1)
      if (key.eq.8) tau = amax1(srelpr,prgopt(last+2))
      if (key.eq.9) blowup = amax1(srelpr,prgopt(last+2))
      next = prgopt(link)
      if (.not.(next.le.0 .or. next.gt.nlink)) go to 150
      nerr = 3
      iopt = 1
      call xerror(39hwnnls( ) the option vector is undefined, 39, nerr,
     * iopt)
      mode = 2
      return
  150 last = link
      link = next
      go to 110
  160 do 170 j=1,n
        call sscal(m, d(j), w(1,j), 1)
  170 continue
      go to 1260
  180 continue
c
c     to initialize-variables
c
c     srelpr is the precision for the particular machine
c     being used.  this logic avoids recomputing it every entry.
      if (.not.(srelpr.eq.zero)) go to 210
      srelpr = one
  190 if (one+srelpr.eq.one) go to 200
      srelpr = srelpr/two
      go to 190
  200 srelpr = srelpr*two
  210 m = ma + mme
      me = mme
      mep1 = me + 1
      assign 220 to igo977
      go to 100
c
c     process-option-vector
  220 done = .false.
      iter = 0
      itmax = 3*(n-l)
      mode = 0
      lp1 = l + 1
      nsoln = l
      nsp1 = nsoln + 1
      np1 = n + 1
      nm1 = n - 1
      l1 = min0(m,l)
c
c     compute scale factor to apply to equal. constraint equas.
      do 230 j=1,n
        wd(j) = sasum(m,w(1,j),1)
  230 continue
      imax = isamax(n,wd,1)
      eanorm = wd(imax)
      bnorm = sasum(m,w(1,np1),1)
      alamda = eanorm/(srelpr*fac)
c
c     define scaling diag matrix for mod givens usage and
c     classify equation types.
      if(alamda .gt. 1.0e10)then
        alsq = 1.0e20
        alamda = 1.0e10
      else 
         alsq = alamda * alamda
      endif
      do 260 i=1,m
c
c     when equ i is heavily weighted itype(i)=0, else itype(i)=1.
        if (.not.(i.le.me)) go to 240
        t = alsq
        itemp = 0
        go to 250
  240   t = one
        itemp = 1
  250   scale(i) = t
        itype(i) = itemp
  260 continue
c
c     set the soln vector x(*) to zero and the col interchange
c     matrix to the identity.
      x(1) = zero
      call scopy(n, x, 0, x, 1)
      do 270 i=1,n
        ipivot(i) = i
  270 continue
      go to 1230
  280 continue
c
c     to initially-triangularize
c
c     set first l comps. of dual vector to zero because
c     these correspond to the unconstrained variables.
      if (.not.(l.gt.0)) go to 290
      wd(1) = zero
      call scopy(l, wd, 0, wd, 1)
c
c     the arrays idope(*) and dope(*) are used to pass
c     information to wnlit().  this was done to avoid
c     a long calling sequence or the use of common.
  290 idope(1) = me
      idope(2) = mep1
      idope(3) = 0
      idope(4) = 1
      idope(5) = nsoln
      idope(6) = 0
      idope(7) = 1
      idope(8) = l1
c
      dope(1) = alsq
      dope(2) = eanorm
      dope(3) = fac
      dope(4) = tau
      call wnlit(w, mdw, m, n, l, ipivot, itype, h, scale, rnorm,
     * idope, dope, done)
      me = idope(1)
      mep1 = idope(2)
      krank = idope(3)
      krp1 = idope(4)
      nsoln = idope(5)
      niv = idope(6)
      niv1 = idope(7)
      l1 = idope(8)
      go to 1240
  300 continue
c
c     to compute-search-direction-and-feasible-point
c
c     solve the triangular system of currently non-active
c     variables and store the solution in z(*).
c
c     solve-system
      assign 310 to igo958
      go to 1110
c
c     increment iteration counter and check against max. number
c     of iterations.
  310 iter = iter + 1
      if (.not.(iter.gt.itmax)) go to 320
      mode = 1
      done = .true.
c
c     check to see if any constraints have become active.
c     if so, calculate an interpolation factor so that all
c     active constraints are removed from the basis.
  320 alpha = two
      hitcon = .false.
      if (.not.(l.lt.nsoln)) go to 360
      do 350 j=lp1,nsoln
        zz = z(j)
        if (.not.(zz.le.zero)) go to 340
        t = x(j)/(x(j)-zz)
        if (.not.(t.lt.alpha)) go to 330
        alpha = t
        jcon = j
  330   hitcon = .true.
  340   continue
  350 continue
  360 go to 1220
  370 continue
c
c     to add-constraints
c
c     use computed alpha to interpolate between last
c     feasible solution x(*) and current unconstrained
c     (and infeasible) solution z(*).
      if (.not.(lp1.le.nsoln)) go to 390
      do 380 j=lp1,nsoln
        x(j) = x(j) + alpha*(z(j)-x(j))
  380 continue
  390 feasbl = .false.
      go to 410
  400 if (feasbl) go to 610
c
c     remove col jcon and shift cols jcon+1 through n to the
c     left. swap col jcon into the n-th position.  this achieves
c     upper hessenberg form for the nonactive constraints and
c     leaves an upper hessenberg matrix to retriangularize.
  410 do 420 i=1,m
        t = w(i,jcon)
        call scopy(n-jcon, w(i,jcon+1), mdw, w(i,jcon), mdw)
        w(i,n) = t
  420 continue
c
c     update permuted index vector to reflect this shift and swap.
      itemp = ipivot(jcon)
      if (.not.(jcon.lt.n)) go to 440
      do 430 i=jcon,nm1
        ipivot(i) = ipivot(i+1)
  430 continue
  440 ipivot(n) = itemp
c
c     similarly repermute x(*) vector.
      call scopy(n-jcon, x(jcon+1), 1, x(jcon), 1)
      x(n) = zero
      nsp1 = nsoln
      nsoln = nsoln - 1
      niv1 = niv
      niv = niv - 1
c
c     retriangularize upper hessenberg matrix after adding constraints.
      j = jcon
      i = krank + jcon - l
  450 if (.not.(j.le.nsoln)) go to 570
      if (.not.(itype(i).eq.0 .and. itype(i+1).eq.0)) go to 470
      assign 460 to igo938
      go to 620
c
c     (itype(i).eq.0 .and. itype(i+1).eq.0) zero-ip1-to-i-in-col-j
  460 go to 560
  470 if (.not.(itype(i).eq.1 .and. itype(i+1).eq.1)) go to 490
      assign 480 to igo938
      go to 620
c
c     (itype(i).eq.1 .and. itype(i+1).eq.1) zero-ip1-to-i-in-col-j
  480 go to 560
  490 if (.not.(itype(i).eq.1 .and. itype(i+1).eq.0)) go to 510
      call sswap(np1, w(i,1), mdw, w(i+1,1), mdw)
      call sswap(1, scale(i), 1, scale(i+1), 1)
      itemp = itype(i+1)
      itype(i+1) = itype(i)
      itype(i) = itemp
c
c     swapped row was formerly a pivot elt., so it will
c     be large enough to perform elim.
      assign 500 to igo938
      go to 620
c
c     zero-ip1-to-i-in-col-j
  500 go to 560
  510 if (.not.(itype(i).eq.0 .and. itype(i+1).eq.1)) go to 550
      t = scale(i)*w(i,j)**2/alsq
      if (.not.(t.gt.tau**2*eanorm**2)) go to 530
      assign 520 to igo938
      go to 620
  520 go to 540
  530 call sswap(np1, w(i,1), mdw, w(i+1,1), mdw)
      call sswap(1, scale(i), 1, scale(i+1), 1)
      itemp = itype(i+1)
      itype(i+1) = itype(i)
      itype(i) = itemp
      w(i+1,j) = zero
  540 continue
  550 continue
  560 i = i + 1
      j = j + 1
      go to 450
c
c     see if the remaining coeffs in the soln set are feasible.  they
c     should be because of the way alpha was determined.  if any are
c     infeasible it is due to roundoff error.  any that are non-
c     positive will be set to zero and removed from the soln set.
  570 if (.not.(lp1.le.nsoln)) go to 590
      do 580 jcon=lp1,nsoln
        if (x(jcon).le.zero) go to 600
  580 continue
  590 feasbl = .true.
  600 continue
      go to 400
  610 go to 1200
  620 continue
c
c     to zero-ip1-to-i-in-col-j
      if (.not.(w(i+1,j).ne.zero)) go to 630
      call srotmg(scale(i), scale(i+1), w(i,j), w(i+1,j), sparam)
      w(i+1,j) = zero
      call srotm(np1-j, w(i,j+1), mdw, w(i+1,j+1), mdw, sparam)
  630 go to 1290
  640 continue
c
c     to perform-multiplier-test-and-drop-a-constraint
      call scopy(nsoln, z, 1, x, 1)
      if (.not.(nsoln.lt.n)) go to 650
      x(nsp1) = zero
      call scopy(n-nsoln, x(nsp1), 0, x(nsp1), 1)
  650 i = niv1
  660 if (.not.(i.le.me)) go to 690
c
c     reclassify least squares eqations as equalities as
c     necessary.
      if (.not.(itype(i).eq.0)) go to 670
      i = i + 1
      go to 680
  670 call sswap(np1, w(i,1), mdw, w(me,1), mdw)
      call sswap(1, scale(i), 1, scale(me), 1)
      itemp = itype(i)
      itype(i) = itype(me)
      itype(me) = itemp
      mep1 = me
      me = me - 1
  680 go to 660
c
c     form inner product vector wd(*) of dual coeffs.
  690 if (.not.(nsp1.le.n)) go to 730
      do 720 j=nsp1,n
        sm = zero
        if (.not.(nsoln.lt.m)) go to 710
        do 700 i=nsp1,m
          sm = sm + scale(i)*w(i,j)*w(i,np1)
  700   continue
  710   wd(j) = sm
  720 continue
  730 go to 750
  740 if (pos .or. done) go to 970
c
c     find j such that wd(j)=wmax is maximum.  this determines
c     that the incoming col j will reduce the residual vector
c     and be positive.
  750 wmax = zero
      iwmax = nsp1
      if (.not.(nsp1.le.n)) go to 780
      do 770 j=nsp1,n
        if (.not.(wd(j).gt.wmax)) go to 760
        wmax = wd(j)
        iwmax = j
  760   continue
  770 continue
  780 if (.not.(wmax.le.zero)) go to 790
      done = .true.
      go to 960
c
c     set dual coeff to zero for incoming col.
  790 wd(iwmax) = zero
c
c     wmax .gt. zero, so okay to move col iwmax to soln set.
c     perform transformation to retriangularize, and test
c     for near linear dependence.
c     swap col iwmax into nsoln-th position to maintain upper
c     hessenberg form of adjacent cols, and add new col to
c     triangular decomposition.
      nsoln = nsp1
      nsp1 = nsoln + 1
      niv = niv1
      niv1 = niv + 1
      if (.not.(nsoln.ne.iwmax)) go to 800
      call sswap(m, w(1,nsoln), 1, w(1,iwmax), 1)
      wd(iwmax) = wd(nsoln)
      wd(nsoln) = zero
      itemp = ipivot(nsoln)
      ipivot(nsoln) = ipivot(iwmax)
      ipivot(iwmax) = itemp
c
c     reduce col nsoln so that the matrix of nonactive
c     constraints variables is triangular.
  800 j = m
  810 if (.not.(j.gt.niv)) go to 870
      jm1 = j - 1
      jp = jm1
c
c     when operating near the me line, test to see if the pivot elt.
c     is near zero.  if so, use the largest elt. above it as the pivot.
c     this is to maintain the sharp interface between weighted and
c     non-weighted rows in all cases.
      if (.not.(j.eq.mep1)) go to 850
      imax = me
      amax = scale(me)*w(me,nsoln)**2
  820 if (.not.(jp.ge.niv)) go to 840
      t = scale(jp)*w(jp,nsoln)**2
      if (.not.(t.gt.amax)) go to 830
      imax = jp
      amax = t
  830 jp = jp - 1
      go to 820
  840 jp = imax
  850 if (.not.(w(j,nsoln).ne.zero)) go to 860
      call srotmg(scale(jp), scale(j), w(jp,nsoln), w(j,nsoln), sparam)
      w(j,nsoln) = zero
      call srotm(np1-nsoln, w(jp,nsp1), mdw, w(j,nsp1), mdw, sparam)
  860 j = jm1
      go to 810
c
c     solve for z(nsoln)=proposed new value for x(nsoln).
c     test if this is nonpositive or too large.
c     if this was true or if the pivot term was zero reject
c     the col as dependent.
  870 if (.not.(w(niv,nsoln).ne.zero)) go to 890
      isol = niv
      assign 880 to igo897
      go to 980
c
c     test-proposed-new-component
  880 go to 940
  890 if (.not.(niv.le.me .and. w(mep1,nsoln).ne.zero)) go to 920
c
c     try to add row mep1 as an additional equality constraint.
c     check size of proposed new soln component.
c     reject it if it is too large.
      isol = mep1
      assign 900 to igo897
      go to 980
c
c     test-proposed-new-component
  900 if (.not.(pos)) go to 910
c
c     swap rows mep1 and niv, and scale factors for these rows.
      call sswap(np1, w(mep1,1), mdw, w(niv,1), mdw)
      call sswap(1, scale(mep1), 1, scale(niv), 1)
      itemp = itype(mep1)
      itype(mep1) = itype(niv)
      itype(niv) = itemp
      me = mep1
      mep1 = me + 1
  910 go to 930
  920 pos = .false.
  930 continue
  940 if (pos) go to 950
      nsp1 = nsoln
      nsoln = nsoln - 1
      niv1 = niv
      niv = niv - 1
  950 continue
  960 go to 740
  970 go to 1250
  980 continue
c
c     to test-proposed-new-component
      z2 = w(isol,np1)/w(isol,nsoln)
      z(nsoln) = z2
      pos = z2.gt.zero
      if (.not.(z2*eanorm.ge.bnorm .and. pos)) go to 990
      pos = .not.(blowup*z2*eanorm.ge.bnorm)
  990 go to 1280
 1000 continue
c     to compute-final-solution
c
c     solve system, store results in x(*).
c
      assign 1010 to igo958
      go to 1110
c     solve-system
 1010 call scopy(nsoln, z, 1, x, 1)
c
c     apply householder transformations to x(*) if krank.lt.l
      if (.not.(0.lt.krank .and. krank.lt.l)) go to 1030
      do 1020 i=1,krank
        call h12(2, i, krp1, l, w(i,1), mdw, h(i), x, 1, 1, 1)
 1020 continue
c
c     fill in trailing zeroes for constrained variables not in soln.
 1030 if (.not.(nsoln.lt.n)) go to 1040
      x(nsp1) = zero
      call scopy(n-nsoln, x(nsp1), 0, x(nsp1), 1)
c
c     repermute soln vector to natural order.
 1040 do 1070 i=1,n
        j = i
 1050   if (ipivot(j).eq.i) go to 1060
        j = j + 1
        go to 1050
 1060   ipivot(j) = ipivot(i)
        ipivot(i) = j
        call sswap(1, x(j), 1, x(i), 1)
 1070 continue
c
c     rescale the soln using the col scaling.
      do 1080 j=1,n
        x(j) = x(j)*d(j)
 1080 continue
      if (.not.(nsoln.lt.m)) go to 1100
      do 1090 i=nsp1,m
        t = w(i,np1)
        if (i.le.me) t = t/alamda
        t = (scale(i)*t)*t
        rnorm = rnorm + t
 1090 continue
 1100 rnorm = sqrt(rnorm)
      go to 1210
c
c     to solve-system
c
 1110 continue
      if (.not.(done)) go to 1120
      isol = 1
      go to 1130
 1120 isol = lp1
 1130 if (.not.(nsoln.ge.isol)) go to 1190
c
c     copy rt. hand side into temp vector to use overwriting method.
      call scopy(niv, w(1,np1), 1, temp, 1)
      do 1180 jj=isol,nsoln
        j = nsoln - jj + isol
        if (.not.(j.gt.krank)) go to 1140
        i = niv - jj + isol
        go to 1150
 1140   i = j
 1150   if (.not.(j.gt.krank .and. j.le.l)) go to 1160
        z(j) = zero
        go to 1170
 1160   z(j) = temp(i)/w(i,j)
        call saxpy(i-1, -z(j), w(1,j), 1, temp, 1)
 1170   continue
 1180 continue
 1190 go to 1270
 1200 go to igo986, (40)
 1210 go to igo980, (90)
 1220 go to igo991, (30)
 1230 go to igo998, (10)
 1240 go to igo995, (20)
 1250 go to igo983, (60)
 1260 go to igo977, (220)
 1270 go to igo958, (310, 1010)
 1280 go to igo897, (880, 900)
 1290 go to igo938, (460, 480, 500, 520)
      end
      subroutine wnlit(w, mdw, m, n, l, ipivot, itype, h, scale, rnorm,
     * idope, dope, done)
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (begin changes at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/scopy/dcopy/,/srotm/drotm/,
c     /sscal/dscal/,
c     /sswap/dswap/,/amax1/dmax1/,/isamax/idamax/,/.e-/.d-/,/e0/d0/
c
c     this is a companion subprogram to wnnls( ).
c     the documentation for wnnls( ) has more complete
c     usage instructions.
c
c     note  the m by (n+1) matrix w( , ) contains the rt. hand side
c           b as the (n+1)st col.
c
c
c     triangularize l1 by l1 subsystem, where l1=min(m,l), with
c     col interchanges.
c     revised march 4, 1982
c
c++
C     real             w(mdw,1), h(1), scale(1), dope(4), sparam(5)
      real             w(mdw,*), h(*), scale(*), dope(4), sparam(5)

      real             alsq, amax, eanorm, fac, factor, hbar, one, rn
      real             rnorm, sn, t, tau, tenm3, zero
      real             amax1

C     integer itype(1), ipivot(1), idope(8)
      integer itype(*), ipivot(*), idope(8)

      integer isamax
      logical indep, done, recalc
      data tenm3 /1.e-3/, zero /0.e0/, one /1.e0/
c
      me = idope(1)
      mep1 = idope(2)
      krank = idope(3)
      krp1 = idope(4)
      nsoln = idope(5)
      niv = idope(6)
      niv1 = idope(7)
      l1 = idope(8)
c
      alsq = dope(1)
      eanorm = dope(2)
      fac = dope(3)
      tau = dope(4)
      np1 = n + 1
      lb = min0(m-1,l)
      recalc = .true.
      rnorm = zero
      krank = 0
c     we set factor=1.e0 so that the heavy weight alamda will be
c     included in the test for col independence.
      factor = 1.e0
      i = 1
      ip1 = 2
      lend = l
   10 if (.not.(i.le.lb)) go to 150
c
c     set ir to point to the i-th row.
      ir = i
      mend = m
      assign 20 to igo996
      go to 460
c
c     update-col-ss-and-find-pivot-col
   20 assign 30 to igo993
      go to 560
c
c     perform-col-interchange
c
c     set ic to point to i-th col.
   30 ic = i
      assign 40 to igo990
      go to 520
c
c     test-indep-of-incoming-col
   40 if (.not.(indep)) go to 110
c
c     eliminate i-th col below diag. using mod. givens transformations
c     applied to (a b).
      j = m
      do 100 jj=ip1,m
        jm1 = j - 1
        jp = jm1
c     when operating near the me line, use the largest elt.
c     above it as the pivot.
        if (.not.(j.eq.mep1)) go to 80
        imax = me
        amax = scale(me)*w(me,i)**2
   50   if (.not.(jp.ge.i)) go to 70
        t = scale(jp)*w(jp,i)**2
        if (.not.(t.gt.amax)) go to 60
        imax = jp
        amax = t
   60   jp = jp - 1
        go to 50
   70   jp = imax
   80   if (.not.(w(j,i).ne.zero)) go to 90
        call srotmg(scale(jp), scale(j), w(jp,i), w(j,i), sparam)
        w(j,i) = zero
        call srotm(np1-i, w(jp,ip1), mdw, w(j,ip1), mdw, sparam)
   90   j = jm1
  100 continue
      go to 140
  110 continue
      if (.not.(lend.gt.i)) go to 130
c
c     col i is dependent. swap with col lend.
      max = lend
c
c     perform-col-interchange
      assign 120 to igo993
      go to 560
  120 continue
      lend = lend - 1
c
c     find col in remaining set with largest ss.
      max = isamax(lend-i+1,h(i),1) + i - 1
      hbar = h(max)
      go to 30
  130 continue
      krank = i - 1
      go to 160
  140 i = ip1
      ip1 = ip1 + 1
      go to 10
  150 krank = l1
  160 continue
      krp1 = krank + 1
      if (.not.(krank.lt.me)) go to 290
      factor = alsq
      do 170 i=krp1,me
        if (l.gt.0) w(i,1) = zero
        call scopy(l, w(i,1), 0, w(i,1), mdw)
  170 continue
c
c     determine the rank of the remaining equality constraint
c     equations by eliminating within the block of constrained
c     variables.  remove any redundant constraints.
      lp1 = l + 1
      recalc = .true.
      lb = min0(l+me-krank,n)
      i = lp1
      ip1 = i + 1
  180 if (.not.(i.le.lb)) go to 280
      ir = krank + i - l
      lend = n
      mend = me
      assign 190 to igo996
      go to 460
c
c     update-col-ss-and-find-pivot-col
  190 assign 200 to igo993
      go to 560
c
c     perform-col-interchange
c
c     eliminate elements in the i-th col.
  200 j = me
  210 if (.not.(j.gt.ir)) go to 230
      jm1 = j - 1
      if (.not.(w(j,i).ne.zero)) go to 220
      call srotmg(scale(jm1), scale(j), w(jm1,i), w(j,i), sparam)
      w(j,i) = zero
      call srotm(np1-i, w(jm1,ip1), mdw, w(j,ip1), mdw, sparam)
  220 j = jm1
      go to 210
c
c     set ic=i=col being eliminated
  230 ic = i
      assign 240 to igo990
      go to 520
c
c     test-indep-of-incoming-col
  240 if (indep) go to 270
c
c     remove any redundant or dependent equality constraints.
      jj = ir
  250 if (.not.(ir.le.me)) go to 260
      w(ir,1) = zero
      call scopy(n, w(ir,1), 0, w(ir,1), mdw)
      rnorm = rnorm + (scale(ir)*w(ir,np1)/alsq)*w(ir,np1)
      w(ir,np1) = zero
      scale(ir) = one
c     reclassify the zeroed row as a least squares equation.
      itype(ir) = 1
      ir = ir + 1
      go to 250
c
c     reduce me to reflect any discovered dependent equality
c     constraints.
  260 continue
      me = jj - 1
      mep1 = me + 1
      go to 300
  270 i = ip1
      ip1 = ip1 + 1
      go to 180
  280 continue
  290 continue
  300 continue
      if (.not.(krank.lt.l1)) go to 420
c
c     try to determine the variables krank+1 through l1 from the
c     least squares equations.  continue the triangularization with
c     pivot element w(mep1,i).
c
      recalc = .true.
c
c     set factor=alsq to remove effect of heavy weight from
c     test for col independence.
      factor = alsq
      kk = krp1
      i = kk
      ip1 = i + 1
  310 if (.not.(i.le.l1)) go to 410
c
c     set ir to point to the mep1-st row.
      ir = mep1
      lend = l
      mend = m
      assign 320 to igo996
      go to 460
c
c     update-col-ss-and-find-pivot-col
  320 assign 330 to igo993
      go to 560
c
c     perform-col-interchange
c
c     eliminate i-th col below the ir-th element.
  330 irp1 = ir + 1
      j = m
      do 350 jj=irp1,m
        jm1 = j - 1
        if (.not.(w(j,i).ne.zero)) go to 340
        call srotmg(scale(jm1), scale(j), w(jm1,i), w(j,i), sparam)
        w(j,i) = zero
        call srotm(np1-i, w(jm1,ip1), mdw, w(j,ip1), mdw, sparam)
  340   j = jm1
  350 continue
c
c     test if new pivot element is near zero. if so, the col is
c     dependent.
      t = scale(ir)*w(ir,i)**2
      indep = t.gt.tau**2*eanorm**2
      if (.not.indep) go to 380
c
c     col test passed. now must pass row norm test to be classified
c     as independent.
      rn = zero
      do 370 i1=ir,m
        do 360 j1=ip1,n
          rn = amax1(rn,scale(i1)*w(i1,j1)**2)
  360   continue
  370 continue
      indep = t.gt.tau**2*rn
c
c     if independent, swap the ir-th and krp1-st rows to maintain the
c     triangular form.  update the rank indicator krank and the
c     equality constraint pointer me.
  380 if (.not.(indep)) go to 390
      call sswap(np1, w(krp1,1), mdw, w(ir,1), mdw)
      call sswap(1, scale(krp1), 1, scale(ir), 1)
c     reclassify the least sq. equation as an equality constraint and
c     rescale it.
      itype(ir) = 0
      t = sqrt(scale(krp1))
      call sscal(np1, t, w(krp1,1), mdw)
      scale(krp1) = alsq
      me = mep1
      mep1 = me + 1
      krank = krp1
      krp1 = krank + 1
      go to 400
  390 go to 430
  400 i = ip1
      ip1 = ip1 + 1
      go to 310
  410 continue
  420 continue
  430 continue
c
c     if pseudorank is less than l, apply householder trans.
c     from right.
      if (.not.(krank.lt.l)) go to 450
      do 440 i=1,krank
        j = krp1 - i
        call h12(1, j, krp1, l, w(j,1), mdw, h(j), w, mdw, 1, j-1)
  440 continue
  450 niv = krank + nsoln - l
      niv1 = niv + 1
      if (l.eq.n) done = .true.
c
c  end of initial triangularization.
      idope(1) = me
      idope(2) = mep1
      idope(3) = krank
      idope(4) = krp1
      idope(5) = nsoln
      idope(6) = niv
      idope(7) = niv1
      idope(8) = l1
      return
  460 continue
c
c     to update-col-ss-and-find-pivot-col
c
c     the col ss vector will be updated at each step. when
c     numerically necessary, these values will be recomputed.
c
      if (.not.(ir.ne.1 .and. (.not.recalc))) go to 480
c     update col ss =sum of squares.
      do 470 j=i,lend
        h(j) = h(j) - scale(ir-1)*w(ir-1,j)**2
  470 continue
c
c     test for numerical accuracy.
      max = isamax(lend-i+1,h(i),1) + i - 1
      recalc = hbar + tenm3*h(max).eq.hbar
c
c     if required, recalculate col ss, using rows ir through mend.
  480 if (.not.(recalc)) go to 510
      do 500 j=i,lend
        h(j) = zero
        do 490 k=ir,mend
          h(j) = h(j) + scale(k)*w(k,j)**2
  490   continue
  500 continue
c
c     find col with largest ss.
      max = isamax(lend-i+1,h(i),1) + i - 1
      hbar = h(max)
  510 go to 600
  520 continue
c
c     to test-indep-of-incoming-col
c
c     test the col ic to determine if it is linearly independent
c     of the cols already in the basis.  in the init tri
c     step, we usually want the heavy weight alamda to
c     be included in the test for independence.  in this case the
c     value of factor will have been set to 1.e0 before this
c     procedure is invoked.  in the potentially rank deficient
c     problem, the value of factor will have been
c     set to alsq=alamda**2 to remove the effect of the heavy weight
c     from the test for independence.
c
c     write new col as partitioned vector
c             (a1)  number of components in soln so far = niv
c             (a2)  m-niv components
c     and compute  sn = inverse weighted length of a1
c                  rn = inverse weighted length of a2
c     call the col independent when rn .gt. tau*sn
      sn = zero
      rn = zero
      do 550 j=1,mend
        t = scale(j)
        if (j.le.me) t = t/factor
        t = t*w(j,ic)**2
        if (.not.(j.lt.ir)) go to 530
        sn = sn + t
        go to 540
  530   rn = rn + t
  540   continue
  550 continue
      indep = rn.gt.tau**2*sn
      go to 590
  560 continue
c
c     to perform-col-interchange
c
      if (.not.(max.ne.i)) go to 570
c     exchange elements of permuted index vector and perform col
c     interchanges.
      itemp = ipivot(i)
      ipivot(i) = ipivot(max)
      ipivot(max) = itemp
      call sswap(m, w(1,max), 1, w(1,i), 1)
      t = h(max)
      h(max) = h(i)
      h(i) = t
  570 go to 580
  580 go to igo993, (30, 200, 330, 120)
  590 go to igo990, (40, 240)
  600 go to igo996, (20, 190, 320)
      end
      subroutine hfti(a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h,
     * g, ip)
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (begin changes at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/ sqrt/ dsqrt/,
c     /, abs/, dabs/,/abs(/dabs(/,/e0/d0/
c
c     dimension a(mda,n),(b(mdb,nb) or b(m)),rnorm(nb),h(n),g(n),ip(n)
c
c     written by c. l. lawson and r. j. hanson.  from the book solving
c     least squares problems, prentice-hall, inc. (1974). for further
c     algorithmic details see algorithm hfti in chapter 14.
c
c     abstract
c
c     this subroutine solves a linear least squares problem or a set of
c     linear least squares problems having the same matrix but different
c     right-side vectors.  the problem data consists of an m by n matrix
c     a, an m by nb matrix b, and an absolute tolerance parameter tau
c     whose usage is described below.  the nb column vectors of b
c     represent right-side vectors for nb distinct linear least squares
c     problems.
c
c     this set of problems can also be written as the matrix least
c     squares problem
c
c                       ax = b,
c
c     where x is the n by nb solution matrix.
c
c     note that if b is the m by m identity matrix, then x will be the
c     pseudo-inverse of a.
c
c     this subroutine first transforms the augmented matrix (a b) to a
c     matrix (r c) using premultiplying householder transformations with
c     column interchanges.  all subdiagonal elements in the matrix r are
c     zero and its diagonal elements satisfy
c
c                       abs(r(i,i)).ge.abs(r(i+1,i+1)),
c
c                       i = 1,...,l-1, where
c
c                       l = min(m,n).
c
c     the subroutine will compute an integer, krank, equal to the number
c     of diagonal terms of r that exceed tau in magnitude.  then a
c     solution of minimum euclidean length is computed using the first
c     krank rows of (r c).
c
c     to be specific we suggest that the user consider an easily
c     computable matrix norm, such as, the maximum of all column sums of
c     magnitudes.
c
c     now if the relative uncertainty of b is eps, (norm of uncertainty/
c     norm of b), it is suggested that tau be set approximately equal to
c     eps*(norm of a).
c
c     the user must dimension all arrays appearing in the call list..
c     a(mda,n),(b(mdb,nb) or b(m)),rnorm(nb),h(n),g(n),ip(n).  this
c     permits the solution of a range of problems in the same array
c     space.
c
c     the entire set of parameters for hfti are
c
c     input..
c
c     a(*,*),mda,m,n    the array a(*,*) initially contains the m by n
c                       matrix a of the least squares problem ax = b.
c                       the first dimensioning parameter of the array
c                       a(*,*) is mda, which must satisfy mda.ge.m
c                       either m.ge.n or m.lt.n is permitted.  there
c                       is no restriction on the rank of a.  the
c                       condition mda.lt.m is considered an error.
c
c     b(*),mdb,nb       if nb = 0 the subroutine will perform the
c                       orthogonal decomposition but will make no
c                       references to the array b(*).  if nb.gt.0
c                       the array b(*) must initially contain the m by
c                       nb matrix b of the least squares problem ax =
c                       b.  if nb.ge.2 the array b(*) must be doubly
c                       subscripted with first dimensioning parameter
c                       mdb.ge.max(m,n).  if nb = 1 the array b(*) may
c                       be either doubly or singly subscripted.  in
c                       the latter case the value of mdb is arbitrary
c                       but it should be set to some valid integer
c                       value such as mdb = m.
c
c                       the condition of nb.gt.1.and.mdb.lt. max(m,n)
c                       is considered an error.
c
c     tau               absolute tolerance parameter provided by user
c                       for pseudorank determination.
c
c     h(*),g(*),ip(*)   arrays of working space used by hfti.
c
c     output..
c
c     a(*,*)            the contents of the array a(*,*) will be
c                       modified by the subroutine.  these contents
c                       are not generally required by the user.
c
c     b(*)              on return the array b(*) will contain the n by
c                       nb solution matrix x.
c
c     krank             set by the subroutine to indicate the
c                       pseudorank of a.
c
c     rnorm(*)          on return, rnorm(j) will contain the euclidean
c                       norm of the residual vector for the problem
c                       defined by the j-th column vector of the array
c                       b(*,*) for j = 1,...,nb.
c
c     h(*),g(*)         on return these arrays respectively contain
c                       elements of the pre- and post-multiplying
c                       householder transformations used to compute
c                       the minimum euclidean length solution.
c
c     ip(*)             array in which the subroutine records indices
c                       describing the permutation of column vectors.
c                       the contents of arrays h(*),g(*) and ip(*)
c                       are not generally required by the user.
c
c++
      real             a(mda,n), b(mdb,1), h(n), g(n), rnorm(nb), tau
      real             factor, hmax, sm1, zero, sm, tmp
      real             diff, sqrt, abs
      integer   ip(n)
      zero = 0.e0
      factor = 0.001e0
c
      k = 0
      ldiag = min0(m,n)
      if (ldiag.le.0) go to 310
      if (.not.mda.lt.m) go to 10
      nerr = 2
      iopt = 2
      call xerror(31hhfti mda.lt.m.. probable error., 31, nerr, iopt)
      return
   10 continue
c
      if (.not.(nb.gt.1 .and. max0(m,n).gt.mdb)) go to 20
      nerr = 2
      iopt = 2
      call xerror(49hhfti mdb.lt.max(m,n).and.nb.gt.1. probable error.,
     * 49, nerr, iopt)
      return
   20 continue
c
      do 100 j=1,ldiag
        if (j.eq.1) go to 40
c
c     update squared column lengths and find lmax
c    ..
        lmax = j
        do 30 l=j,n
          h(l) = h(l) - a(j-1,l)**2
          if (h(l).gt.h(lmax)) lmax = l
   30   continue
        if (diff(hmax+factor*h(lmax),hmax)) 40, 40, 70
c
c     compute squared column lengths and find lmax
c    ..
   40   lmax = j
        do 60 l=j,n
          h(l) = zero
          do 50 i=j,m
            h(l) = h(l) + a(i,l)**2
   50     continue
          if (h(l).gt.h(lmax)) lmax = l
   60   continue
        hmax = h(lmax)
c    ..
c     lmax has been determined
c
c     do column interchanges if needed.
c    ..
   70   continue
        ip(j) = lmax
        if (ip(j).eq.j) go to 90
        do 80 i=1,m
          tmp = a(i,j)
          a(i,j) = a(i,lmax)
          a(i,lmax) = tmp
   80   continue
        h(lmax) = h(j)
   90 jcol = min0(j+1,n)
c
c     compute the j-th transformation and apply it to a and b.
c    ..
        call h12(1, j, j+1, m, a(1,j), 1, h(j), a(1,jcol), 1, mda, n-j)
        call h12(2, j, j+1, m, a(1,j), 1, h(j), b, 1, mdb, nb)
  100 continue
c
c     determine the pseudorank, k, using the tolerance, tau.
c    ..
      do 110 j=1,ldiag
        if (abs(a(j,j)).le.tau) go to 120
  110 continue
      k = ldiag
      go to 130
  120 k = j - 1
  130 kp1 = k + 1
c
c     compute the norms of the residual vectors.
c
      if (nb.le.0) go to 170
      do 160 jb=1,nb
        tmp = zero
        if (kp1.gt.m) go to 150
        do 140 i=kp1,m
          tmp = tmp + b(i,jb)**2
  140   continue
  150   rnorm(jb) = sqrt(tmp)
  160 continue
  170 continue
c                                           special for pseudorank = 0
      if (k.gt.0) go to 200
      if (nb.le.0) go to 310
      do 190 jb=1,nb
        do 180 i=1,n
          b(i,jb) = zero
  180   continue
  190 continue
      go to 310
c
c     if the pseudorank is less than n compute householder
c     decomposition of first k rows.
c    ..
  200 if (k.eq.n) go to 220
      do 210 ii=1,k
        i = kp1 - ii
        call h12(1, i, kp1, n, a(i,1), mda, g(i), a, mda, 1, i-1)
  210 continue
  220 continue
c
c
      if (nb.le.0) go to 310
      do 300 jb=1,nb
c
c     solve the k by k triangular system.
c    ..
        do 250 l=1,k
          sm = zero
          i = kp1 - l
          if (i.eq.k) go to 240
          ip1 = i + 1
          do 230 j=ip1,k
            sm = sm + a(i,j)*b(j,jb)
  230     continue
  240     sm1 = sm
          b(i,jb) = (b(i,jb)-sm1)/a(i,i)
  250   continue
c
c     complete computation of solution vector.
c    ..
        if (k.eq.n) go to 280
        do 260 j=kp1,n
          b(j,jb) = zero
  260   continue
        do 270 i=1,k
          call h12(2, i, kp1, n, a(i,1), mda, g(i), b(1,jb), 1, mdb, 1)
  270   continue
c
c      re-order the solution vector to compensate for the
c      column interchanges.
c    ..
  280   do 290 jj=1,ldiag
          j = ldiag + 1 - jj
          if (ip(j).eq.j) go to 290
          l = ip(j)
          tmp = b(l,jb)
          b(l,jb) = b(j,jb)
          b(j,jb) = tmp
  290   continue
  300 continue
c    ..
c     the solution vectors, x, are now
c     in the first  n  rows of the array b(,).
c
  310 krank = k
      return
      end

C------------------------------------------------------------------

      subroutine h12 (mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)


c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (start changes at line with c++ in cols. 1-3.)
c     /real (12 blanks)/double precision/,/sdot/ddot/,/abs,/dabs,/,
c     /sswap/dswap/,/sqrt/dsqrt/,/abs(/ dabs(/,/amax1/dmax1/,
c     /saxpy/daxpy/,/e0/d0/
c
c
c     c.l.lawson and r.j.hanson, jet propulsion laboratory, 1973 jun 12
c     to appear in 'solving least squares problems', prentice-hall, 1974
c
c     modified at sandia labs., may 1977, to --
c
c     1)  remove double precision accumulation, and
c     2)  include usage of the basic linear algebra package for
c         vectors longer than a particular threshold.
c
c     construction and/or application of a single
c     householder transformation..     q = i + u*(u**t)/b
c
c     mode    = 1 or 2   to select algorithm  h1  or  h2 .
c     lpivot is the index of the pivot element.
c     l1,m   if l1 .le. m   the transformation will be constructed to
c            zero elements indexed from l1 through m.   if l1 gt. m
c            the subroutine does an identity transformation.
c     u(),iue,up    on entry to h1 u() contains the pivot vector.
c                   iue is the storage increment between elements.
c                                       on exit from h1 u() and up
c                   contain quantities defining the vector u of the
c                   householder transformation.   on entry to h2 u()
c                   and up should contain quantities previously computed
c                   by h1.  these will not be modified by h2.
c     c()    on entry to h1 or h2 c() contains a matrix which will be
c            regarded as a set of vectors to which the householder
c            transformation is to be applied.  on exit c() contains the
c            set of transformed vectors.
c     ice    storage increment between elements of vectors in c().
c     icv    storage increment between vectors in c().
c     ncv    number of vectors in c() to be transformed. if ncv .le. 0
c            no operations will be done on c().
c
c     subroutine h12 (mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
c++


C     real             u(iue,m), c(1), up
      real             u(iue,m), c(*), up

      real             b, cl, clinv, one, sm, ul1m1
      real             abs, amax1, sqrt, sdot
      one=1.e0
c
      if (0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
      cl=abs(u(1,lpivot))
      if (mode.eq.2) go to 60
c                            ****** construct the transformation. ******
          do 10 j=l1,m
   10     cl=amax1(abs(u(1,j)),cl)
      if (cl) 130,130,20
   20 clinv=one/cl
      sm=(u(1,lpivot)*clinv)**2
          do 30 j=l1,m
   30     sm=sm+(u(1,j)*clinv)**2
      cl=cl*sqrt(sm)
      if (u(1,lpivot)) 50,50,40
   40 cl=-cl
   50 up=u(1,lpivot)-cl
      u(1,lpivot)=cl
      go to 70
c            ****** apply the transformation  i+u*(u**t)/b  to c. ******
c
   60 if (cl) 130,130,70
   70 if (ncv.le.0) return
      b=up*u(1,lpivot)
c                       b  must be nonpositive here.  if b = 0., return.
c
      if (b) 80,130,130
   80 b=one/b
      mml1p2=m-l1+2
      if (mml1p2.gt.20) go to 140
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
          do 120 j=1,ncv
          i2=i2+icv
          i3=i2+incr
          i4=i3
          sm=c(i2)*up
              do 90 i=l1,m
              sm=sm+c(i3)*u(1,i)
   90         i3=i3+ice
          if (sm) 100,120,100
  100     sm=sm*b
          c(i2)=c(i2)+sm*up
              do 110 i=l1,m
              c(i4)=c(i4)+sm*u(1,i)
  110         i4=i4+ice
  120     continue
  130 return
  140 continue
      l1m1=l1-1
      kl1=1+(l1m1-1)*ice
      kl2=kl1
      klp=1+(lpivot-1)*ice
      ul1m1=u(1,l1m1)
      u(1,l1m1)=up
      if (lpivot.eq.l1m1) go to 150
      call sswap(ncv,c(kl1),icv,c(klp),icv)
  150 continue
          do 160 j=1,ncv
          sm=sdot(mml1p2,u(1,l1m1),iue,c(kl1),ice)
          sm=sm*b
          call saxpy (mml1p2,sm,u(1,l1m1),iue,c(kl1),ice)
          kl1=kl1+icv
  160 continue
      u(1,l1m1)=ul1m1
      if (lpivot.eq.l1m1) return
      kl1=kl2
      call sswap(ncv,c(kl1),icv,c(klp),icv)
      return
      end

      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
