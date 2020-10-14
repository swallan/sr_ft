      DOUBLE PRECISION FUNCTION PPND(P,IER)
C
C       ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
C
C       PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA = P.
C
C	See also AS 241 which contains alternative routines accurate to
C	about 7 and 16 decimal digits.
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DATA SPLIT /0.42D0/
        DATA A0,A1,A2,A3/2.50662823884D0,-18.61500062529D0,
     1  41.39119773534D0,-25.44106049637D0/, B1,B2,B3,B4/
     2  -8.47351093090D0,23.08336743743D0,-21.06224101826D0,
     3  3.13082909833D0/, C0,C1,C2,C3/-2.78718931138D0,-2.29796479134D0,
     4  4.85014127135D0,2.32121276858D0/, D1,D2/3.54388924762D0,
     5  1.63706781897D0/
      DATA ZERO/0.D0/, ONE/1.D0/, HALF/0.5D0/
C
      IER = 0
      Q = P-HALF
      IF (ABS(Q).GT.SPLIT) GO TO 10
C
C       0.08 < P < 0.92
C
      R = Q*Q
      PPND = Q*(((A3*R + A2)*R + A1)*R + A0)/((((B4*R + B3)*R + B2)*R
     1      + B1)*R + ONE)
      RETURN
C
C       P < 0.08 OR P > 0.92, SET R = MIN(P,1-P)
C
   10 R = P
      IF (Q.GT.ZERO) R = ONE-P
      IF (R.LE.ZERO) GO TO 20
      R = SQRT(-LOG(R))
      PPND = (((C3*R + C2)*R + C1)*R + C0)/((D2*R + D1)*R + ONE)
      IF (Q.LT.ZERO) PPND = -PPND
      RETURN
   20   IER = 1
      PPND = ZERO
      RETURN
      END
c     This file includes the Applied Statistics algorithm AS 66 for calculating
c     the tail area under the normal curve, and two alternative routines which
c     give higher accuracy.   The latter have been contributed by Alan Miller of
c     CSIRO Division of Mathematics & Statistics, Clayton, Victoria.   Notice
c     that each function or routine has different call arguments.
c
c     Allan Miller suggests that NPROB (fortran adaptation of Adam's routine) is more accurate than Hill's ALNORM.  This is because in the process of restructuring Hill's code, Miller inadvertently introduced typographical errors into two of Hill's constants.
c
c     Miller gives his variable q the value of
c     0.39990348504d0
c     which corresponds to Hill's A2 with value
c      0.399903438504
c     This typo reduces the accuracy of Miller's ALNORM by up to 4 decimal digits for |x|<1.28.  The omission of a digit would have been visually obvious in Hill's original layout.
c  
c     Miller gives his variable b1 the value of
c       -29.8213557807d0
c     but in Hill's original variable A4, the final digit was 8, not 7.
c
c     With these two corrected constants, the double precision versions of ALNORM and NPROB give 
C     identical results (as Hill said they should), except for values beyond Hill's large tail cutoff 
C     value of LTONE (or beyond the small tail cutoff of 12.7 that Miller introduced into Adam's routine).
C     Hill said the cutoff values should be altered for double precision.  I suggest LTONE=8.3 and
C     UTZERO=37.52, which like Hill's original cutoffs merely avoid useless work without limiting the working range.
c 
c     Jerry W. Lewis, PhD,  JWLewis53@verizon.net
c
      double precision function alnorm(x,upper)
c
c         Algorithm AS66 Applied Statistics (1973) vol22 no.3
c
c       Evaluates the tail area of the standardised normal curve
c       from x to infinity if upper is .true. or
c       from minus infinity to x if upper is .false.
c
      double precision zero,one,half
      double precision con,z,y,x
      double precision p,q,r,a1,a2,a3,b1,b2,c1,c2,c3,c4,c5,c6
      double precision d1,d2,d3,d4,d5
      logical upper,up
c*** machine dependent constants
      double precision ltone,utzero
      data zero/0.0d0/, one/1.0d0/, half/0.5d0/
      data ltone/7.0d0/,utzero/18.66d0/
      data con/1.28d0/
      data p/0.398942280444d0/,q/0.39990348504d0/,r/0.398942280385d0/   
      data a1/5.75885480458d0/,a2/2.62433121679d0/,a3/5.92885724438d0/  
      data b1/-29.8213557807d0/,b2/48.6959930692d0/
      data c1/-3.8052d-8/,c2/3.98064794d-4/,c3/-0.151679116635d0/
      data c4/4.8385912808d0/,c5/0.742380924027d0/,c6/3.99019417011d0/  
      data d1/1.00000615302d0/,d2/1.98615381364d0/,d3/5.29330324926d0/  
      data d4/-15.1508972451d0/,d5/30.789933034d0/
c
      up=upper
      z=x
      if(z.ge.zero)goto 10
      up=.not.up
      z=-z
   10 if(z.le.ltone.or.up.and.z.le.utzero)goto 20
      alnorm=zero
      goto 40
   20 y=half*z*z
      if(z.gt.con) goto 30
c
      alnorm=half-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
      goto 40
   30 alnorm=r*dexp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/   
     2   (z+c6))))))
   40 if(.not.up)alnorm=one-alnorm
      return
      end
C
C    
      SUBROUTINE NPROB(Z,P,Q,PDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       P, Q = PROBABILITIES TO THE LEFT AND RIGHT OF Z
C       FOR THE STANDARD NORMAL DISTRIBUTION.
C       PDF  = THE PROBABILITY DENSITY FUNCTION
C
C       REFERENCE: ADAMS,A.G. AREAS UNDER THE NORMAL CURVE,
C       ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969.
C
C       LATEST REVISION - 23 JANUARY 1981
C
C          **************************************************************
C
      DATA A0,A1,A2,A3,A4,A5,A6,A7/0.5D0,0.398942280444D0,
     1 0.399903438504D0, 5.75885480458D0, 29.8213557808D0,
     2 2.62433121679D0, 48.6959930692D0, 5.92885724438D0/,
     3 B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11/0.398942280385D0,
     4 3.8052D-8, 1.00000615302D0, 3.98064794D-4, 1.98615381364D0,
     5 0.151679116635D0, 5.29330324926D0, 4.8385912808D0,
     6 15.1508972451D0, 0.742380924027D0, 30.789933034D0,
     7 3.99019417011D0/
C
      ZABS = ABS(Z)
      IF(ZABS.GT.12.7D0) GO TO 20
      Y = A0*Z*Z
      PDF = EXP(-Y)*B0
      IF(ZABS.GT.1.28D0) GO TO 10
C
C       Z BETWEEN -1.28 AND +1.28
C
      Q = A0-ZABS*(A1-A2*Y/(Y+A3-A4/(Y+A5+A6/(Y+A7))))
      IF(Z.LT.0.D0) GO TO 30
      P = 1.D0-Q
      RETURN
C
C       ZABS BETWEEN 1.28 AND 12.7
C
   10 Q = PDF/(ZABS-B1+B2/(ZABS+B3+B4/(ZABS-B5+B6/(ZABS+B7-B8/
     1  (ZABS+B9+B10/(ZABS+B11))))))
      IF(Z.LT.0.D0) GO TO 30
      P = 1.D0-Q
      RETURN
C
C       Z FAR OUT IN TAIL
C
   20 Q = 0.D0
      PDF = 0.D0
      IF (Z.LT.0.D0) GO TO 30
      P = 1.D0
      RETURN
C
C       NEGATIVE Z, INTERCHANGE P AND Q
C
   30 P = Q
      Q = 1.D0-P
      RETURN
      END
c
c
c
C
      double precision function prtrng(q, v, r, ifault)
      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
c	 Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)
c
c        Evaluates the probability from 0 to q for a studentized
c        range having v degrees of freedom and r samples.
c
c        Uses subroutine ALNORM = algorithm AS66.
c
c        Arrays vw and qw store transient values used in the
c        quadrature summation.  Node spacing is controlled by
c        step.  pcutj and pcutk control truncation.
c        Minimum and maximum number of steps are controlled by
c        jmin, jmax, kmin and kmax.  Accuracy can be increased
c        by use of a finer grid - Increase sizes of arrays vw
c        and qw, and jmin, jmax, kmin, kmax and 1/step proportionally.
c
      double precision q,v,r,vw(30), qw(30), pcutj, pcutk, step,
     #    vmax,zero,fifth,half,one,two,cv1,cv2,cvmax,cv(4)
      double precision g, gmid, r1, c, h, v2, gstep, pk1, pk2, gk, pk
      double precision w0, pz, x, hj, ehj, pj
      data pcutj, pcutk, step, vmax /0.00003d0, 0.0001d0, 0.45d0,
     #   120.0d0/, zero, fifth, half, one, two /0.0d0, 0.2d0, 0.5d0,
     #   1.0d0, 2.0d0/, cv1, cv2, cvmax /0.193064705d0, 0.293525326d0,
     #   0.39894228d0/, cv(1), cv(2), cv(3), cv(4) /0.318309886d0,
     #   -0.268132716d-2, 0.347222222d-2, 0.833333333d-1/
      data jmin, jmax, kmin, kmax/3, 15, 7, 15/
c
c        Check initial values
c
      prtrng = zero
      ifault = 0
      if (v .lt. one .or. r .lt. two) ifault = 1
      if (q .le. zero .or. ifault .eq. 1) goto 99
c
c        Computing constants, local midpoint, adjusting steps.
c
      g = step * r ** (-fifth)
      gmid = half * log(r)
      r1 = r - one
      c = log(r * g * cvmax)
      if(v.gt.vmax) goto 20
c
      h = step * v ** (-half)
      v2 = v * half
      if (v .eq. one) c = cv1
      if (v .eq. two) c = cv2
      if (.not. (v .eq. one .or. v .eq. two)) c = sqrt(v2)
     #    * cv(1) / (one + ((cv(2) / v2 + cv(3)) / v2 + cv(4)) / v2)
      c = log(c * r * g * h)
c
c        Computing integral
c        Given a row k, the procedure starts at the midpoint and works
c        outward (index j) in calculating the probability at nodes
c        symmetric about the midpoint.  The rows (index k) are also
c        processed outwards symmetrically about the midpoint.  The
c        centre row is unpaired.
c
   20 gstep = g
      qw(1) = -one
      qw(jmax + 1) = -one
      pk1 = one
      pk2 = one
      do 28 k = 1, kmax
      gstep = gstep - g
   21   gstep = -gstep
      gk = gmid + gstep
      pk = zero
      if (pk2 .le. pcutk .and. k .gt. kmin) goto 26
      w0 = c - gk * gk * half
        pz = alnorm(gk, .true.)
        x = alnorm(gk - q, .true.) - pz
      if (x .gt. zero) pk = exp(w0 + r1 * log(x))
      if (v .gt. vmax) goto 26
c
      jump = -jmax
   22   jump = jump + jmax
      do 24 j = 1, jmax
      jj = j + jump
      if (qw(jj) .gt. zero) goto 23
      hj = h * j
      if (j .lt. jmax) qw(jj + 1) = -one
      ehj = exp(hj)
      qw(jj) = q * ehj
      vw(jj) = v * (hj + half - ehj * ehj * half)
c
   23 pj = zero
       x = alnorm(gk - qw(jj), .true.) - pz
      if (x .gt. zero) pj = exp(w0 + vw(jj) + r1 * log(x))
      pk = pk + pj
      if (pj .gt. pcutj) goto 24
      if (jj .gt. jmin .or. k .gt. kmin) goto 25
   24   continue
   25   h = -h
      if (h .lt. zero) goto 22
c
   26   prtrng = prtrng + pk
      if (k .gt. kmin .and. pk .le. pcutk .and. pk1 .le. pcutk)goto 99
      pk2 = pk1
      pk1 = pk
      if (gstep .gt. zero) goto 21
   28 continue
c
   99 return
      end
c
c
      double precision function qtrng(p, v, r, ifault)
      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190.1  Appl. Statist. (1983) Vol.32, No. 2
c
c        Approximates the quantile p for a studentized range
c        distribution having v degrees of freedom and r samples
c        for probability p, p.ge.0.90 .and. p.le.0.99.
c
c        Uses functions  alnorm, ppnd, prtrng and qtrng0 -
c        Algorithms AS 66, AS 111, AS 190 and AS 190.2
c
      double precision p, v, r, pcut, p75, p80, p90, p99, p995
      double precision p175, one, two, five
      double precision q1, p1, q2, p2, e1, e2
      double precision eps
      data jmax, pcut, p75, p80, p90, p99, p995, p175, one, two, five
     #  /8, 0.001d0, 0.75d0, 0.80d0, 0.90d0, 0.99d0, 0.995d0, 1.75d0,
     #  1.0d0, 2.0d0, 5.0d0/
      data eps/1.0d-4/
c
c        Check input parameters
c
      ifault = 0
      nfault = 0
      if (v .lt. one .or. r.lt. two) ifault = 1
      if (p .lt. p90 .or. p .gt. p99) ifault = 2
      if (ifault .ne. 0) goto 99
c
c        Obtain initial values
c
      q1 = qtrng0(p, v, r, nfault)
      if (nfault .ne. 0) goto 99
      p1 = prtrng(q1, v, r, nfault)
      if (nfault .ne. 0) goto 99
      qtrng = q1
      if (abs(p1-p) .lt. pcut) goto 99
      if (p1 .gt. p) p1 = p175 * p - p75 * p1
      if (p1 .lt. p) p2 = p + (p - p1) * (one - p) / (one - p1) * p75
      if (p2 .lt. p80) p2 = p80
      if (p2 .gt. p995) p2 = p995
      q2 = qtrng0(p2, v, r, nfault)
      if (nfault .ne. 0) goto 99
c
c        Refine approximation
c
      do 14 j = 2, jmax
      p2 = prtrng(q2, v, r, nfault)
      if (nfault .ne. 0) goto 99
      e1 = p1 - p
      e2 = p2 - p
      qtrng = (q1 + q2) / two
      d = e2 - e1
      if (abs(d) .gt. eps) qtrng = (e2 * q1 - e1 * q2) / d
      if(abs(e1) .lt. abs(e2)) goto 12
      q1 = q2
      p1 = p2
   12   if (abs(p1 - p) .lt. pcut * five) goto 99
      q2 = qtrng
   14 continue
c
   99 if (nfault .ne. 0) ifault = 9
      return
      end
c
c
      double precision function qtrng0(p, v, r, ifault)
      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
c
c        Calculates an initial quantile p for a studentized range
c        distribution having v degrees of freedom and r samples
c        for probability p, p.gt.0.80 .and. p.lt.0.995.
c
c        Uses function ppnd - Algorithm AS 111
c
      double precision p, v, r, q, t, vmax, half, one, four, c1, c2, c3
      double precision c4, c5
      data vmax, half, one, four, c1, c2, c3, c4, c5 / 120.0d0, 0.5d0,
     #  1.0d0, 4.0d0, 0.8843d0, 0.2368d0, 1.214d0, 1.208d0, 1.4142d0/
c
      t=ppnd(half + half * p,ifault)
      if (v .lt. vmax) t = t + (t * t* t + t) / v / four
      q = c1 - c2 * t
      if (v .lt. vmax) q = q - c3 / v + c4 * t / v
      qtrng0 = t * (q * log(r - one) + c5)
      return
      end
      
