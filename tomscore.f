
C 20060607 tbt  output was undefined. Set
C 20060607 tbt Changed       data gam,     gamsq,   rgam,     rgamsq 


      real             function diff(x,y)
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c     use an editing command (change) /string-1/(to)string-2/.
c     (apply changes to entire program unit.)
c     /real (12 blanks)/double precision/
c
c     c.l.lawson and r.j.hanson, jet propulsion laboratory, 1973 june 7
c     to appear in 'solving least squares problems', prentice-hall, 1974
      real              x, y
      diff=x-y
      return
      end



      subroutine srotm (n,sx,incx,sy,incy,sparam)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c
c     (sx(1)     sx(n))
c     (      ...      )
c     (sy(1)     sy(n))
c
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c
      dimension sx(1),sy(1),sparam(5)
      data zero,two /0.e0,2.e0/
c
      sflag=sparam(1)
      if(n .le. 0 .or.(sflag+two.eq.zero)) go to 140
          if(.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if(sflag) 50,10,30
   10          continue
               sh12=sparam(4)
               sh21=sparam(3)
                    do 20 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w+z*sh12
                    sy(i)=w*sh21+z
   20               continue
               go to 140
   30          continue
               sh11=sparam(2)
               sh22=sparam(5)
                    do 40 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z
                    sy(i)=-w+sh22*z
   40               continue
               go to 140
   50          continue
               sh11=sparam(2)
               sh12=sparam(4)
               sh21=sparam(3)
               sh22=sparam(5)
                    do 60 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z*sh12
                    sy(i)=w*sh21+z*sh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if(incx .lt. 0) kx=1+(1-n)*incx
          if(incy .lt. 0) ky=1+(1-n)*incy
c
          if(sflag)120,80,100
   80     continue
          sh12=sparam(4)
          sh21=sparam(3)
               do 90 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w+z*sh12
               sy(ky)=w*sh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          sh11=sparam(2)
          sh22=sparam(5)
               do 110 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z
               sy(ky)=-w+sh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          sh11=sparam(2)
          sh12=sparam(4)
          sh21=sparam(3)
          sh22=sparam(5)
               do 130 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z*sh12
               sy(ky)=w*sh21+z*sh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
          end
      subroutine srotmg (sd1,sd2,sx1,sy1,sparam)
c
c     construct the modified givens transformation matrix h which zeros
c     the second component of the 2-vector  (sqrt(sd1)*sx1,sqrt(sd2)*
c     sy2)**t.
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c
c end of abstract
      dimension sparam(5)
c
      data zero,one,two /0.e0,1.e0,2.e0/,iflag/1/

      data gam,     gamsq,   rgam,     rgamsq 
     1    /4096.E00, 1.678E07, 2.441E-04, 5.960E-08/

C    1    /4096.e0, 1.678e7, 2.441e-4, 5.960e-8/
c
      if(.not. sd1 .lt. zero) go to 10
c       go zero-h-d-and-sx1..
          go to 60
   10 continue
c     case-sd1-nonnegative
      sp2=sd2*sy1
      if(.not. sp2 .eq. zero) go to 20
          sflag=-two
          go to 260
c     regular-case..
   20 continue
      sp1=sd1*sx1
      sq2=sp2*sy1
      sq1=sp1*sx1
c
      if(.not. abs(sq1) .gt. abs(sq2)) go to 40
          sh21=-sy1/sx1
          sh12=sp2/sp1
c
          su=one-sh12*sh21
c
          if(.not. su .le. zero) go to 30
c         go zero-h-d-and-sx1..
               go to 60
   30     continue
               sflag=zero
               sd1=sd1/su
               sd2=sd2/su
               sx1=sx1*su
c         go scale-check..
               go to 100
   40 continue
          if(.not. sq2 .lt. zero) go to 50
c         go zero-h-d-and-sx1..
               go to 60
   50     continue
               sflag=one
               sh11=sp1/sp2
               sh22=sx1/sy1
               su=one+sh11*sh22
               stemp=sd2/su
               sd2=sd1/su
               sd1=stemp
               sx1=sy1*su
c         go scale-check
               go to 100
c     procedure..zero-h-d-and-sx1..
   60 continue
          sflag=-one
          sh11=zero
          sh12=zero
          sh21=zero
          sh22=zero
c
          sd1=zero
          sd2=zero
          sx1=zero
c         return..
          go to 220
c     procedure..fix-h..
   70 continue
      if(.not. sflag .ge. zero) go to 90
c
          if(.not. sflag .eq. zero) go to 80
          sh11=one
          sh22=one
          sflag=-one
          go to 90
   80     continue
          sh21=-one
          sh12=one
          sflag=-one
   90 continue
      go to igo,(120,150,180,210)
c     procedure..scale-check
  100 continue
               if(.not. iflag.eq.1) go to 105
c
c                   recompute rescaling parameters
c                   more accurately..
c
                    rgam = one/gam
                    gamsq = gam**2
                    rgamsq = rgam**2
                    iflag = 2
  105          continue
  110     continue
          if(.not. sd1 .le. rgamsq) go to 130
               if(sd1 .eq. zero) go to 160
               assign 120 to igo
c              fix-h..
               go to 70
  120          continue
               sd1=sd1*gamsq
               sx1=sx1*rgam
               sh11=sh11*rgam
               sh12=sh12*rgam
          go to 110
  130 continue
  140     continue
          if(.not. sd1 .ge. gamsq) go to 160
               assign 150 to igo
c              fix-h..
               go to 70
  150          continue
               sd1=sd1*rgamsq
               sx1=sx1*gam
               sh11=sh11*gam
               sh12=sh12*gam
          go to 140
  160 continue
  170     continue
          if(.not. abs(sd2) .le. rgamsq) go to 190
               if(sd2 .eq. zero) go to 220
               assign 180 to igo
c              fix-h..
               go to 70
  180          continue
               sd2=sd2*gamsq
               sh21=sh21*rgam
               sh22=sh22*rgam
          go to 170
  190 continue
  200     continue
          if(.not. abs(sd2) .ge. gamsq) go to 220
               assign 210 to igo
c              fix-h..
               go to 70
  210          continue
               sd2=sd2*rgamsq
               sh21=sh21*gam
               sh22=sh22*gam
          go to 200
  220 continue
          if(sflag)250,230,240
  230     continue
               sparam(3)=sh21
               sparam(4)=sh12
               go to 260
  240     continue
               sparam(2)=sh11
               sparam(5)=sh22
               go to 260
  250     continue
               sparam(2)=sh11
               sparam(3)=sh21
               sparam(4)=sh12
               sparam(5)=sh22
  260 continue
          sparam(1)=sflag
          return
      end
      subroutine fdump
c     abstract
c        ***note*** machine dependent routine
c        fdump is intended to be replaced by a locally written
c        version which produces a symbolic dump.  failing this,
c        it should be replaced by a version which prints the
c        subprogram nesting list.  note that this dump must be
c        printed on each of up to five files, as indicated by the
c        xgetua routine.  see xsetua and xgetua for details.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  23 may 1979
c
      return
      end
      subroutine xerabt(messg,nmessg)
c
c     abstract
c        ***note*** machine dependent routine
c        xerabt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg and nmessg are as in xerror, except that nmessg may
c        be zero, in which case no message is being supplied.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      dimension messg(nmessg)
      if (.true.) stop
      return
      end



      function j4save(iwhich,ivalue,iset)
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                 = 1 refers to current error number.
c                 = 2 refers to current error control flag.
c                 = 3 refers to current unit number to which error
c                     messages are to be sent.  (0 means use standard.)
c                 = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                 = 5 refers to the total number of units to which
c                     each error message is to be written.
c                 = 6 refers to the 2nd unit for error messages
c                 = 7 refers to the 3rd unit for error messages
c                 = 8 refers to the 4th unit for error messages
c                 = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c     written by ron jones, with slatec common math library subcommittee
c     adapted from bell laboratories port library error handler
c end of abstract
c     latest revision ---  23 may 1979
c
      logical iset
      integer iparam(9)
 
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/

      j4save  = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue

      return
      end




      function numxer(nerr)
c
c     abstract
c        numxer returns the most recent error number,
c        in both numxer and the parameter nerr.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      nerr = j4save(1,0,.false.)
      numxer = nerr
      return
      end
      subroutine s88fmt(n,ivalue,ifmt)
c
c     abstract
c        s88fmt replaces ifmt(1), ... ,ifmt(n) with the
c        characters corresponding to the n least significant
c        digits of ivalue.
c
c     taken from the bell laboratories port library error handler
c end of abstract
c     latest revision ---  7 june 1978
c
      dimension ifmt(n),idigit(10)
      data idigit(1),idigit(2),idigit(3),idigit(4),idigit(5),
     1     idigit(6),idigit(7),idigit(8),idigit(9),idigit(10)
     2     /1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
      nt = n
      it = ivalue
   10    if (nt .eq. 0) return
         index = mod(it,10)
         ifmt(nt) = idigit(index+1)
         it = it/10
         nt = nt - 1
         go to 10
      end
      subroutine xerclr
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary to do in order to determine that
c        a certain error has occurred again since the last time
c        numxer was referenced.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      junk = j4save(1,0,.true.)
      return
      end
      subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made) a call is made to xerctl.
c        if the user has provided his own version of xerctl, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        messg1 - the first word (only) of the error message.
c        nmessg - same as in the call to xerror or xerrwv.
c        nerr   - same as in the call to xerror or xerrwv.
c        level  - same as in the call to xerror or xerrwv.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c end of abstract
c
      return
      end
      subroutine xerdmp
c
c     abstract
c        xerdmp prints an error table showing all errors which
c        have occurred during the current execution, or since xerdmp
c        was last called.  after printing, the error table is cleared,
c        and if program execution is continued accumulation of the
c        error table begins at zero.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      call xersav(1h ,0,0,0,kount)
      return
      end
      subroutine xermax(max)
c
c     abstract
c        xermax sets the maximum number of times any message
c        is to be printed.  that is, non-fatal messages are
c        not to be printed after they have occured max times.
c        such non-fatal messages may be printed less than
c        max times even if they occur max times, if error
c        suppression mode (kontrl=0) is ever in effect.
c
c        the default value for max is 10.
c
c     description of parameter
c      --input--
c        max - the maximum number of times any one message
c              is to be printed.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      junk = j4save(4,max,.true.)
      return
      end
      subroutine xerprt(messg,nmessg)
c
c     abstract
c        print the hollerith message in messg, of length messg,
c        on each file indicated by xgetua.
c        this version prints exactly the right number of characters,
c        not a number of words, and thus should work on machines
c        which do not blank fill the last word of the hollerith.
c
c     ron jones, june 1980
c end of abstract
c
      integer f(10),g(14),lun(5)
      dimension messg(nmessg)
      data f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10)
     1   / 1h( ,1h1 ,1hx ,1h, ,1h  ,1h  ,1ha ,1h  ,1h  ,1h) /
      data g(1),g(2),g(3),g(4),g(5),g(6),g(7),g(8),g(9),g(10)
     1   / 1h( ,1h1 ,1hx ,1h  ,1h  ,1h  ,1h  ,1h  ,1h  ,1h  /
      data g(11),g(12),g(13),g(14)
     1   / 1h   ,1h   ,1h   ,1h)  /
      data la/1ha/,lcom/1h,/,lblank/1h /
c     prepare format for whole lines
      nchar = i1mach(6)
      nfield = 72/nchar
      call s88fmt(2,nfield,f(5))
      call s88fmt(2,nchar,f(8))
c     prepare format for last, partial line, if needed
      ncharl = nfield*nchar
      nlines = nmessg/ncharl
      nword  = nlines*nfield
      nchrem = nmessg - nlines*ncharl
      if (nchrem.le.0) go to 40
         do 10 i=4,13
10          g(i) = lblank
         nfield = nchrem/nchar
         if (nfield.le.0) go to 20
c        prepare whole word fields
            g(4) = lcom
            call s88fmt(2,nfield,g(5))
            g(7) = la
            call s88fmt(2,nchar,g(8))
20       continue
         nchlst = mod(nchrem,nchar)
         if (nchlst.le.0) go to 30
c        prepare partial word field
            g(10) = lcom
            g(11) = la
            call s88fmt(2,nchlst,g(12))
30       continue
40    continue
c     print the message
      nword1 = nword+1
      nword2 = (nmessg+nchar-1)/nchar
      call xgetua(lun,nunit)
      do 50 kunit = 1,nunit
         iunit = lun(kunit)
         if (iunit.eq.0) iunit = i1mach(4)
         if (nword.gt.0) write (iunit,f) (messg(i),i=1,nword)
         if (nchrem.gt.0) write (iunit,g) (messg(i),i=nword1,nword2)
50    continue
      return
      end
      subroutine xerror(messg,nmessg,nerr,level)
c
c     abstract
c        xerror processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed, containing
c                no more than 72 characters.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c
c     examples
c        call xerror(23hsmooth -- num was zero.,23,1,2)
c        call xerror(43hinteg  -- less than full accuracy achieved.,
c                    43,2,1)
c        call xerror(65hrooter -- actual zero of f found before interval
c    1 fully collapsed.,65,3,0)
c        call xerror(39hexp    -- underflows being set to zero.,39,1,-1)
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     revised by k haskell to check input args, 2/18/80
c
      dimension messg(nmessg)
c     check for valid input
      lkntrl = j4save (2,0,.false.)
      if (nmessg.gt.0) go to 10
      if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (33hxerror -- nmessg must be positive,33)
      if (lkntrl.gt.0) call fdump
      if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1 29)
      if (lkntrl.gt.0) call xersav (1h ,0,0,0,kdummy)
      call xerabt (23hxerror -- invalid input,23)
      return
   10 continue
      if (nerr.ne.0) go to 15
      if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (28hxerror -- nerr=0 is an error,28)
      if (lkntrl.gt.0) call fdump
      if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1 29)
      if (lkntrl.gt.0) call xersav (1h ,0,0,0,kdummy)
      call xerabt (23hxerror -- invalid input,23)
      return
   15 continue
      if ((level.ge.(-1)).and.(level.le.2)) go to 20
      if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (32hxerror -- invalid value of level,32)
      if (lkntrl.gt.0) call fdump
      if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1 29)
      if (lkntrl.gt.0) call xersav (1h ,0,0,0,kdummy)
      call xerabt (23hxerror -- invalid input,23)
      return
   20 continue
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)
      return
      end
      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
c
c     abstract
c        xerrwv processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c        in addition, up to two integer values and two real
c        values may be printed along with the message.
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c        ni    - number of integer values to be printed. (o to 2)
c        i1    - first integer value.
c        i2    - second integer value.
c        nr    - number of real values to be printed. (0 to 2)
c        r1    - first real value.
c        r2    - second real value.
c
c     examples
c        call xerrwv(29hsmooth -- num (=i1) was zero.,29,1,2,
c    1   1,num,0,0,0.,0.)
c        call xerrwv(54hquadxy -- requested error (r1) less than minimum
c    1 (r2).,54,77,1,0,0,0,2,errreq,errmin)
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  19 mar 1980
c     revised by k haskell to check input args, 2/18/80
c
      dimension messg(nmessg),lun(5)
c     get flags
      lkntrl = j4save(2,0,.false.)
      maxmes = j4save(4,0,.false.)
c     check for valid input
      if (nmessg.gt.0) go to 2
         if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (33hxerrwv -- nmessg must be positive,33)
      if (lkntrl.gt.0) call fdump
      if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1 29)
      if (lkntrl.gt.0) call xersav(1h ,0,0,0,kdummy)
      call xerabt (23hxerrwv -- invalid input,23)
      return
    2 continue
      if (nerr.ne.0) go to 4
      if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (28hxerrwv -- nerr=0 is an error,28)
      if (lkntrl.gt.0) call fdump
      if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1 29)
      if (lkntrl.gt.0) call xersav(1h ,0,0,0,kdummy)
      call xerabt (23hxerrwv -- invalid input,23)
      return
    4 continue
      if ((level.ge.(-1)).and.(level.le.2)) go to 10
      if (lkntrl.gt.0) call xerprt(17hfatal error in...,17)
      call xerprt (32hxerrwv -- invalid value of level,32)
         if (lkntrl.gt.0) call fdump
         if (lkntrl.gt.0) call xerprt(29hjob abort due to fatal error.,
     1   29)
         if (lkntrl.gt.0) call xersav(1h ,0,0,0,kdummy)
         call xerabt(23hxerror -- invalid input,23)
         return
   10 continue
c     record message
      junk = j4save(1,nerr,.true.)
      call xersav(messg,nmessg,nerr,level,kount)
c     let user override
      lfirst = messg(1)
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
c     reset to original values
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max0(-2,min0(2,lkntrl))
      mkntrl = iabs(lkntrl)
c     decide whether to print message
      if ((llevel.lt.2).and.(lkntrl.eq.0)) go to 100
      if (((llevel.eq.(-1)).and.(kount.gt.min0(1,maxmes)))
     1.or.((llevel.eq.0)   .and.(kount.gt.maxmes))
     2.or.((llevel.eq.1)   .and.(kount.gt.maxmes).and.(mkntrl.eq.1))
     3.or.((llevel.eq.2)   .and.(kount.gt.max0(1,maxmes)))) go to 100
         if (lkntrl.le.0) go to 20
            call xerprt(1h ,1)
c           introduction
            if (llevel.eq.(-1)) call xerprt
     1(57hwarning message...this message will only be printed once.,57)
            if (llevel.eq.0) call xerprt(13hwarning in...,13)
            if (llevel.eq.1) call xerprt
     1      (23hrecoverable error in...,23)
            if (llevel.eq.2) call xerprt(17hfatal error in...,17)
   20    continue
c        message
         call xerprt(messg,lmessg)
         call xgetua(lun,nunit)
         do 50 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
            if (ni.ge.1) write (iunit,22) i1
            if (ni.ge.2) write (iunit,23) i2
            if (nr.ge.1) write (iunit,24) r1
            if (nr.ge.2) write (iunit,25) r2
   22       format (11x,21hin above message, i1=,i10)
   23       format (11x,21hin above message, i2=,i10)
   24       format (11x,21hin above message, r1=,e20.10)
   25       format (11x,21hin above message, r2=,e20.10)
            if (lkntrl.le.0) go to 40
c              error number
               write (iunit,30) lerr
   30          format (15h error number =,i10)
   40       continue
   50    continue
c        trace-back
         call fdump
  100 continue
      ifatal = 0
      if ((llevel.eq.2).or.((llevel.eq.1).and.(mkntrl.eq.2)))
     1ifatal = 1
c     quit here if message is not fatal
      if (ifatal.le.0) return
      if (lkntrl.le.0) go to 120
c        print reason for abort
         if (llevel.eq.1) call xerprt
     1   (35hjob abort due to unrecovered error.,35)
         if (llevel.eq.2) call xerprt
     1   (29hjob abort due to fatal error.,29)
c        print error summary
         call xersav(1h ,0,0,0,kdummy)
  120 continue
c     abort
      if ((llevel.eq.2).and.(kount.gt.max0(1,maxmes))) lmessg = 0
      call xerabt(messg,lmessg)
      return
      end
      subroutine xersav(messg,nmessg,nerr,level,icount)
c
c     abstract
c        record that this error occurred.
c
c     description of parameters
c     --input--
c       messg, nmessg, nerr, level are as in xerror,
c       except that when nmessg=0 the tables will be
c       dumped and cleared, and when nmessg is less than zero the
c       tables will be dumped and not cleared.
c     --output--
c       icount will be the number of times this message has
c       been seen, or zero if the table has overflowed and
c       does not contain this message specifically.
c       when nmessg=0, icount will not be altered.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c
      integer f(17),lun(5)
      dimension messg(1)
      dimension mestab(10),nertab(10),levtab(10),kount(10)
c     next three data statements are needed merely to satisfy
c     certain conventions for compilers which dynamically
c     allocate storage.
      data mestab(1),mestab(2),mestab(3),mestab(4),mestab(5),
     1     mestab(6),mestab(7),mestab(8),mestab(9),mestab(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data nertab(1),nertab(2),nertab(3),nertab(4),nertab(5),
     1     nertab(6),nertab(7),nertab(8),nertab(9),nertab(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data levtab(1),levtab(2),levtab(3),levtab(4),levtab(5),
     1     levtab(6),levtab(7),levtab(8),levtab(9),levtab(10)
     2     /0,0,0,0,0,0,0,0,0,0/
c     next two data statements are necessary to provide a blank
c     error table initially
      data kount(1),kount(2),kount(3),kount(4),kount(5),
     1     kount(6),kount(7),kount(8),kount(9),kount(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data kountx/0/
c     next data statement sets up output format
      data f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),
     1     f(11),f(12),f(13),f(14),f(15),f(16),f(17)
     2     /1h( ,1h1 ,1hx ,1h, ,1ha ,1h  ,1h  ,1h, ,1hi ,1h  ,
     3      1h  ,1h, ,1h2 ,1hi ,1h1 ,1h0 ,1h) /
      if (nmessg.gt.0) go to 80
c     dump the table
         if (kount(1).eq.0) return
c        prepare format
         nchar = i1mach(6)
         call s88fmt(2,nchar,f(6))
         ncol = 20 - nchar
         call s88fmt(2,ncol,f(10))
c        print to each unit
         call xgetua(lun,nunit)
         do 60 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c           print table header
            write (iunit,10)
   10       format (32h0          error message summary/
     1              41h first word      nerr     level     count)
c           print body of table
            do 20 i=1,10
               if (kount(i).eq.0) go to 30
               write (iunit,f) mestab(i),nertab(i),levtab(i),kount(i)
   20       continue
   30       continue
c           print number of other errors
            if (kountx.ne.0) write (iunit,40) kountx
   40       format (41h0other errors not individually tabulated=,i10)
            write (iunit,50)
   50       format (1x)
   60    continue
         if (nmessg.lt.0) return
c        clear the error tables
         do 70 i=1,10
   70       kount(i) = 0
         kountx = 0
         return
   80 continue
c     process a message...
c     search for this messg, or else an empty slot for this messg,
c     or else determine that the error table is full.
      do 90 i=1,10
         ii = i
         if (kount(i).eq.0) go to 110
         if (messg(1).ne.mestab(i)) go to 90
         if (nerr.ne.nertab(i)) go to 90
         if (level.ne.levtab(i)) go to 90
         go to 100
   90 continue
c     three possible cases...
c     table is full
         kountx = kountx+1
         icount = 1
         return
c     message found in table
  100    kount(ii) = kount(ii) + 1
         icount = kount(ii)
         return
c     empty slot found for new message
  110    mestab(ii) = messg(1)
         nertab(ii) = nerr
         levtab(ii) = level
         kount(ii)  = 1
         icount = 1
         return
      end
      subroutine xgetf(kontrl)
c
c     abstract
c        xgetf returns the current value of the error control flag
c        in kontrl.  see subroutine xsetf for flag value meanings.
c        (kontrl is an output parameter only.)
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      kontrl = j4save(2,0,.false.)
      return
      end
      subroutine xgetua(iunit,n)
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n.lt.5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c
      dimension iunit(5)
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunit(i) = j4save(index,0,.false.)
   30 continue
      return
      end
      subroutine xgetun(iunit)
c
c     abstract
c        xgetun gets the (first) output file to which error messages
c        are being sent.  to find out if more than one file is being
c        used, one must use the xgetua routine.
c
c     description of parameter
c      --output--
c        iunit - the logical unit number of the  (first) unit to
c                which error messages are being sent.
c                a value of zero means that the default file, as
c                defined by the i1mach routine, is being used.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision --- 23 may 1979
c
      iunit = j4save(3,0,.false.)
      return
      end
      subroutine xsetf(kontrl)
c
c     abstract
c        xsetf sets the error control flag value to kontrl.
c        (kontrl is an input parameter only.)
c        the following table shows how each message is treated,
c        depending on the values of kontrl and level.  (see xerror
c        for description of level.)
c
c        if kontrl is zero or negative, no information other than the
c        message itself (including numeric values, if any) will be
c        printed.  if kontrl is positive, introductory messages,
c        trace-backs, etc., will be printed in addition to the message.
c
c              iabs(kontrl)
c        level        0              1              2
c        value
c          2        fatal          fatal          fatal
c
c          1     not printed      printed         fatal
c
c          0     not printed      printed        printed
c
c         -1     not printed      printed        printed
c                                  only           only
c                                  once           once
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  23 may 1979
c
      if ((kontrl.ge.(-2)).and.(kontrl.le.2)) go to 10
         call xerrwv(39hxsetf  -- invalid value of kontrl (i1).,33,1,2,
     1   1,kontrl,0,0,0.,0.)
         return
   10 junk = j4save(2,kontrl,.true.)
      return
      end
      subroutine xsetua(iunit,n)
c
c     abstract
c        xsetua may be called to declare a list of up to five
c        logical units, each of which is to receive a copy of
c        each error message processed by this package.
c        the purpose of xsetua is to allow simultaneous printing
c        of each error message on, say, a main output file,
c        an interactive terminal, and other files such as graphics
c        communication files.
c
c     description of parameters
c      --input--
c        iunit - an array of up to five unit numbers.
c                normally these numbers should all be different.
c                (but duplicates are not prohibited.)
c        n     - the number of unit numbers provided in iunit.
c                must have 1 .le. n .le. 5.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision --- 23 may 1979
c
      dimension iunit(5)
      if ((n.ge.1).and.(n.le.5)) go to 10
         call xerrwv(34hxsetua -- invalid value of n (i1).,34,1,2,
     1   1,n,0,0,0.,0.)
         return
   10 continue
      do 20 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         junk = j4save(index,iunit(i),.true.)
   20 continue
      junk = j4save(5,n,.true.)
      return
      end
      subroutine xsetun(iunit)
c
c     abstract
c        xsetun sets the output file to which error messages are to
c        be sent.  only one file will be used.  see xsetua for
c        how to declare more than one file.
c
c     description of parameter
c      --input--
c        iunit - an input parameter giving the logical unit number
c                to which error messages are to be sent.
c
c     written by ron jones, with slatec common math library subcommittee
c end of abstract
c     latest revision ---  7 june 1978
c
      junk = j4save(3,iunit,.true.)
      junk = j4save(5,1,.true.)
      return
      end
      integer function i1mach(i)
c***begin prologue  i1mach
c***revision date  811015   (yymmdd)
c***category no.  q
c***keywords  machine constants,integer
c***date written   1975
c***author fox p.a.,hall a.d.,schryer n.l. (bell labs)
c***purpose
c returns integer machine dependent constants
c***description
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   these machine constant routines must be activated for
c   a particular environment.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     i1mach can be used to obtain machine-dependent parameters
c     for the local machine environment.  it is a function
c     subroutine with one (input) argument, and can be called
c     as follows, for example
c
c          k = i1mach(i)
c
c     where i=1,...,16.  the (output) value of k above is
c     determined by the (input) value of i.  the results for
c     various values of i are discussed below.
c
c  i/o unit numbers.
c    i1mach( 1) = the standard input unit.
c    i1mach( 2) = the standard output unit.
c    i1mach( 3) = the standard punch unit.
c    i1mach( 4) = the standard error message unit.
c
c  words.
c    i1mach( 5) = the number of bits per integer storage unit.
c    i1mach( 6) = the number of characters per integer storage unit.
c
c  integers.
c    assume integers are represented in the s-digit, base-a form
c
c               sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c               where 0 .le. x(i) .lt. a for i=0,...,s-1.
c    i1mach( 7) = a, the base.
c    i1mach( 8) = s, the number of base-a digits.
c    i1mach( 9) = a**s - 1, the largest magnitude.
c
c  floating-point numbers.
c    assume floating-point numbers are represented in the t-digit,
c    base-b form
c               sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c               where 0 .le. x(i) .lt. b for i=1,...,t,
c               0 .lt. x(1), and emin .le. e .le. emax.
c    i1mach(10) = b, the base.
c
c  single-precision
c    i1mach(11) = t, the number of base-b digits.
c    i1mach(12) = emin, the smallest exponent e.
c    i1mach(13) = emax, the largest exponent e.
c
c  double-precision
c    i1mach(14) = t, the number of base-b digits.
c    i1mach(15) = emin, the smallest exponent e.
c    i1mach(16) = emax, the largest exponent e.
      integer imach(16),output
      data imach( 1) /   5/
      data imach( 2) /   6/
      data imach( 3) /   6/
      data imach( 4) /   7/
      data imach( 5) /  31/
      data imach( 6) /  4/
      data imach(7) /    2 /
      data imach(8) /   31 /
      data imach(9) / 47483647 /
      data imach(10)/    2 /
      data imach(11)/   31 /
      data imach(12)/ -37 /
      data imach(13)/  38 /
      data imach(14)/   31 /
      data imach(15)/ -37 /
      data imach(16)/  38 /
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.  also, the values of
c  i1mach(1) - i1mach(4) should be checked for consistency
c  with the local operating system.
c
c***references
c  fox p.a., hall a.d., schryer n.l.,*framework for a portable library*,
c  acm transaction on mathematical software, vol. 4, no. 2,
c  june 1978, pp. 177-188.
c***routines called  xerror
c***end prologue  i1mach
c
c
c     equivalence (imach(4),output)
c
c     machine constants for the burroughs 1700 system.
c
c     data imach( 1) /    7 /
c     data imach( 2) /    2 /
c     data imach( 3) /    2 /
c     data imach( 4) /    2 /
c     data imach( 5) /   36 /
c     data imach( 6) /    4 /
c     data imach( 7) /    2 /
c     data imach( 8) /   33 /
c     data imach( 9) / z1ffffffff /
c     data imach(10) /    2 /
c     data imach(11) /   24 /
c     data imach(12) / -256 /
c     data imach(13) /  255 /
c     data imach(14) /   60 /
c     data imach(15) / -256 /
c     data imach(16) /  255 /
c
c     machine constants for the burroughs 5700 system.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  48 /
c     data imach( 6) /   6 /
c     data imach( 7) /   2 /
c     data imach( 8) /  39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /   8 /
c     data imach(11) /  13 /
c     data imach(12) / -50 /
c     data imach(13) /  76 /
c     data imach(14) /  26 /
c     data imach(15) / -50 /
c     data imach(16) /  76 /
c
c     machine constants for the burroughs 6700/7700 systems.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  48 /
c     data imach( 6) /   6 /
c     data imach( 7) /   2 /
c     data imach( 8) /  39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /   8 /
c     data imach(11) /  13 /
c     data imach(12) / -50 /
c     data imach(13) /  76 /
c     data imach(14) /  26 /
c     data imach(15) / -32754 /
c     data imach(16) /  32780 /
c
c     machine constants for the cdc 6000/7000 series.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    7 /
c     data imach( 4) /6loutput/
c     data imach( 5) /   60 /
c     data imach( 6) /   10 /
c     data imach( 7) /    2 /
c     data imach( 8) /   48 /
c     data imach( 9) / 00007777777777777777b /
c     data imach(10) /    2 /
c     data imach(11) /   48 /
c     data imach(12) / -974 /
c     data imach(13) / 1070 /
c     data imach(14) /   96 /
c     data imach(15) / -927 /
c     data imach(16) / 1070 /
c
c     machine constants for the cray 1
c
c     data imach( 1) /   100 /
c     data imach( 2) /   101 /
c     data imach( 3) /   102 /
c     data imach( 4) /   101 /
c     data imach( 5) /    64 /
c     data imach( 6) /     8 /
c     data imach( 7) /     2 /
c     data imach( 8) /    63 /
c     data imach( 9) /  777777777777777777777b /
c     data imach(10) /     2 /
c     data imach(11) /    48 /
c     data imach(12) / -8192 /
c     data imach(13) /  8191 /
c     data imach(14) /    96 /
c     data imach(15) / -8192 /
c     data imach(16) /  8191 /
c
c     machine constants for the data general eclipse s/200
c
c     data imach( 1) /   11 /
c     data imach( 2) /   12 /
c     data imach( 3) /    8 /
c     data imach( 4) /   10 /
c     data imach( 5) /   16 /
c     data imach( 6) /    2 /
c     data imach( 7) /    2 /
c     data imach( 8) /   15 /
c     data imach( 9) /32767 /
c     data imach(10) /   16 /
c     data imach(11) /    6 /
c     data imach(12) /  -64 /
c     data imach(13) /   63 /
c     data imach(14) /   14 /
c     data imach(15) /  -64 /
c     data imach(16) /   63 /
c
c     machine constants for the harris 220
c
c     data imach( 1) /       5 /
c     data imach( 2) /       6 /
c     data imach( 3) /       0 /
c     data imach( 4) /       6 /
c     data imach( 5) /      24 /
c     data imach( 6) /       3 /
c     data imach( 7) /       2 /
c     data imach( 8) /      23 /
c     data imach( 9) / 8388607 /
c     data imach(10) /       2 /
c     data imach(11) /      23 /
c     data imach(12) /    -127 /
c     data imach(13) /     127 /
c     data imach(14) /      38 /
c     data imach(15) /    -127 /
c     data imach(16) /     127 /
c
c     machine constants for the honeywell 600/6000 series.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /   43 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    6 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -127 /
c     data imach(13) /  127 /
c     data imach(14) /   63 /
c     data imach(15) / -127 /
c     data imach(16) /  127 /
c
c     machine constants for the hp 2100
c     3 word double precision option with ftn4
c
c     data imach(1) /      5/
c     data imach(2) /      6 /
c     data imach(3) /      4 /
c     data imach(4) /      1 /
c     data imach(5) /     16 /
c     data imach(6) /      2 /
c     data imach(7) /      2 /
c     data imach(8) /     15 /
c     data imach(9) /  32767 /
c     data imach(10)/      2 /
c     data imach(11)/     23 /
c     data imach(12)/   -128 /
c     data imach(13)/    127 /
c     data imach(14)/     39 /
c     data imach(15)/   -128 /
c     data imach(16)/    127 /
c
c     machine constants for the hp 2100
c     4 word double precision option with ftn4
c
c     data imach(1) /      5 /
c     data imach(2) /      6 /
c     data imach(3) /      4 /
c     data imach(4) /      1 /
c     data imach(5) /     16 /
c     data imach(6) /      2 /
c     data imach(7) /      2 /
c     data imach(8) /     15 /
c     data imach(9) /  32767 /
c     data imach(10)/      2 /
c     data imach(11)/     23 /
c     data imach(12)/   -128 /
c     data imach(13)/    127 /
c     data imach(14)/     55 /
c     data imach(15)/   -128 /
c     data imach(16)/    127 /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9, the sel systems 85/86, and
c     the perkin elmer (interdata) 7/32.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  32 /
c     data imach( 6) /   4 /
c     data imach( 7) /   2 /
c     data imach( 8) /  31 /
c     data imach( 9) / z7fffffff /
c     data imach(10) /  16 /
c     data imach(11) /   6 /
c     data imach(12) / -64 /
c     data imach(13) /  63 /
c     data imach(14) /  14 /
c     data imach(15) / -64 /
c     data imach(16) /  63 /
c
c     machine constants for the pdp-10 (ka processor).
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    5 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   54 /
c     data imach(15) / -101 /
c     data imach(16) /  127 /
c
c     machine constants for the pdp-10 (ki processor).
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    5 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   62 /
c     data imach(15) / -128 /
c     data imach(16) /  127 /
c
c     machine constants for pdp-11 fortran[s supporting
c     32-bit integer arithmetic.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   32 /
c     data imach( 6) /    4 /
c     data imach( 7) /    2 /
c     data imach( 8) /   31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /    2 /
c     data imach(11) /   24 /
c     data imach(12) / -127 /
c     data imach(13) /  127 /
c     data imach(14) /   56 /
c     data imach(15) / -127 /
c     data imach(16) /  127 /
c
c     machine constants for pdp-11 fortran[s supporting
c     16-bit integer arithmetic.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   16 /
c     data imach( 6) /    2 /
c     data imach( 7) /    2 /
c     data imach( 8) /   15 /
c     data imach( 9) / 32767 /
c     data imach(10) /    2 /
c     data imach(11) /   24 /
c     data imach(12) / -127 /
c     data imach(13) /  127 /
c     data imach(14) /   56 /
c     data imach(15) / -127 /
c     data imach(16) /  127 /
c
c     machine constants for the univac 1100 series. ftn compiler
c
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    1 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    4 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   60 /
c     data imach(15) /-1024 /
c     data imach(16) / 1023 /
c
c     machine constants for the univac 1100 series. for compiler
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    7 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    6 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   60 /
c     data imach(15) /-1024/
c     data imach(16) / 1023 /
c
c
c     machine constants for the vax 11/780
c
c     data imach(1) /    5 /
c     data imach(2) /    6 /
c     data imach(3) /    5 /
c     data imach(4) /    6 /
c     data imach(5) /   32 /
c     data imach(6) /    4 /
c     data imach(7) /    2 /
c     data imach(8) /   31 /
c     data imach(9) /2147483647 /
c     data imach(10)/    2 /
c     data imach(11)/   24 /
c     data imach(12)/ -127 /
c     data imach(13)/  127 /
c     data imach(14)/   56 /
c     data imach(15)/ -127 /
c     data imach(16)/  127 /
c
c***first executable statement  i1mach
c
      if (i .lt. 1  .or.  i .gt. 16) go to 10
c
      i1mach=imach(i)
      return
c
   10 continue

C 20060607 tbt  output was undefined. See above
      output = imach(4)

      write(output,9000)
9000  format(39h1error    1 in i1mach - i out of bounds )
c
c     call fdump
c
c
      stop
      end
