c   See Copenhaver, Margaret DiPonzio & Holland, Burt S.: 
c    Multiple Comparisons of Simple Effects in the Two-Way Analysis

c    of Variance with Fixed Effects.  Journal of Statistical Computation
c    and Simulation, vol.30, pp.1-15, 1988.
c
c   FUNCTION SUBPROGRAM CV                                              
c        
c   USES SECANT METHOD TO FIND CRITICAL VALUES                          
c        
c        
c   P = CONFIDENCE LEVEL (1 - ALPHA)                                    
c   RR = NO. OF ROWS OR GROUPS                                          
c   CC = NO. OF COLUMNS OR TREATMENTS                                   
c   DF = DEGREES OF FREEDOM OF ERROR TERM                               
c        
c   IR(1) = ERROR FLAG = 1 IF WPROB PROBABILITY > 1                     
c   IR(2) = ERROR FLAG = 1 IF QPROB PROBABILITY > 1                     
c   IR(3) = ERROR FLAG = 1 IF CONVERGENCE NOT REACHED IN 50 ITERATIONS  
c                      = 2 IF DF < 2                                    
c        
c   CV = RETURNED CRITICAL VALUE                                        
c        
c   PROGRAM WILL NOT TERMINATE IF IR(1) OR IR(2) ARE RAISED.            
c   PROGRAM WILL TERMINATE IF IR(3) IS RAISED.                          
c        
c                                                                       
c        
      function cv(p, rr, cc, df, ir)
c                                                                       
c        
c   IF DIFFERENCE BETWEEN SUCCESSIVE ITERATES < EPS, SEARCH IS TERMINATE
cD       
      implicit double precision (o-z, a-h)
c                                                                       
c        
      parameter (eps = 0.0001d0)
c                                                                       
c  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c        
c   DF MUST BE > 1                                                      
c        
 9    dimension ir(3), it(2)
      if (df .lt. 2.0) then
      ir(3) = 2
      goto 900
      end if
      iter = 0
c                                                                       
c        
c  INITIAL VALUE USING USER-WRITTEN FUNCTION                            
c        
      ir(1) = 0
      ir(2) = 0
      ir(3) = 0
c                                                                       
c        
c   FIND PROB(VALUE < X0)                                               
c        
 4    x0 = qinv(p,cc,df)
      valx0 = qprob(x0,rr,cc,df,it) - p
      if (it(1) .eq. 1) ir(1) = 1
c                                                                       
c        
c   FIND SECOND ITERATE AND PROB(VALUE < X1)                            
c        
c   IF FIRST ITERATE HAS PROBABILITY VALUE EXCEEDING P THEN SECOND      
c        
c    ITERATE IS 1 LESS THAN FIRST ITERATE;   ELSE IT IS 1 GREATER.      
c        
      if (it(2) .eq. 1) ir(2) = 1
      if (valx0 .gt. 0.0) then
      x1 = dmax1(0.0d0,x0 - 1.0)
      else
      x1 = x0 + 1.0
c                                                                       
c        
      end if
      valx1 = qprob(x1,rr,cc,df,it) - p
      if (it(1) .eq. 1) ir(1) = 1
c                                                                       
c        
c   FIND NEW ITERATE                                                    
c        
      if (it(2) .eq. 1) ir(2) = 1
   50 cv = x1 - ((valx1 * (x1 - x0)) / (valx1 - valx0))
      valx0 = valx1
c                                                                       
c        
c   NEW ITERATE MUST BE >= 0                                            
c        
      x0 = x1
      if (cv .lt. 0.0) then
      cv = 0.0
      valx1 = - p
c                                                                       
c        
c     FIND PROB(VALUE < NEW ITERATE)                                    
c        
      end if
      valx1 = qprob(cv,rr,cc,df,it) - p
      if (it(1) .eq. 1) ir(1) = 1
c                                                                       
c        
      if (it(2) .eq. 1) ir(2) = 1
      x1 = cv
      iter = iter + 1
      if (iter .eq. 51) goto 100
c   IF THE DIFFERENCE BETWEEN TWO SUCCESSIVE ITERATES < EPSILON, STOP   
c        
      xabs = dabs(x1 - x0)
      if (xabs .lt. eps) goto 900
c                                                                       
c        
      goto 50
  100 ir(3) = 1
  900 return 
c                                                                       
c        
      end
c                                                                       
c        
c   FUNCTION SUBPROGRAM QPROB                                           
c        
c                                                                       
c        
c   Q = VALUE OF STUDENTIZED RANGE                                      
c        
c   RR = NO. OF ROWS OR GROUPS                                          
c        
c   CC = NO. OF COLUMNS OR TREATMENTS                                   
c        
c   DF = DEGREES OF FREEDOM OF ERROR TERM                               
c        
c   IR(1) = ERROR FLAG = 1 IF WPROB PROBABILITY > 1                     
c        
c   IR(2) = ERROR FLAG = 1 IF QPROB PROBABILITY > 1                     
c        
c   QPROB = RETURNED PROBABILITY INTEGRAL FROM (0, Q)                   
c        
c   PROGRAM WILL NOT TERMINATE IF ER(1) OR ER(2) ARE RAISED.            
c        
c                                                                       
c        
c   ALL REFERENCES IN WPROB AND QLGAMA TO ABRAMOWITZ AND STEGUN         
c        
c   ARE FROM THE FOLLOWING REFERENCE:                                   
c        
c    ABRAMOWITZ, MILTON AND STEGUN, IRENE A.  HANDBOOK OF MATHEMATICAL  
c        
c    FUNCTIONS.  NEW YORK:  DOVER PUBLICATIONS, INC. (1970).            
c        
c    ALL CONSTANTS TAKEN FROM THIS TEXT ARE GIVEN TO 25 SIGNIFICANT DIGI
cTS.     
c                                                                       
c        
      function qprob(q, rr, cc, df, ir)
c                                                                       
c        
c   NLEGQ = ORDER OF LEGENDRE QUADRATURE                                
c        
c   IHALFQ = INT ((NLEGQ + 1) / 2)                                      
c        
c   EPS = MAX. ALLOWABLE VALUE OF INTEGRAL                              
c        
c   EPS1 & EPS2 = VALUES BELOW WHICH THERE IS NO CONTRIBUTION TO INTEGRA
cL.      
c   D.F. <= DHAF:   INTEGRAL IS DIVIDED INTO ULEN1 LENGTH INTERVALS.  EL
cSE      
c   D.F. <= DQUAR:  INTEGRAL IS DIVIDED INTO ULEN2 LENGTH INTERVALS.  EL
cSE      
c   D.F. <= DEIGH:  INTEGRAL IS DIVIDED INTO ULEN3 LENGTH INTERVALS.  EL
cSE      
c   D.F. <= DLARG:  INTEGRAL IS DIVIDED INTO ULEN4 LENGTH INTERVALS.    
c        
c   D.F. > DLARG:   THE RANGE IS USED TO CALCULATE INTEGRAL.            
c        
      implicit double precision (o-z, a-h)
c                                                                       
c        
      parameter (nlegq = 16, ihalfq = 8, eps = 1.0d0, eps1 = -30.0d0, 
     &eps2 = 1.0d-14, dhaf = 100.0d0, dquar = 800.0d0, deigh = 5000.0d0
     &, dlarg = 25000.0d0, ulen1 = 1.0d0, ulen2 = 0.5d0, ulen3 = 0.25d0
     &, ulen4 = 0.125d0)
      dimension ir(2)
c                                                                       
c        
c   R2 = LOG(2)                                                         
c        
c   XLEGQ = LEGENDRE 16-POINT NODES                                     
c        
c   ALEGQ = LEGENDRE 16-POINT COEFFICIENTS                              
c        
      dimension xlegq(ihalfq), alegq(ihalfq)
c                                                                       
c        
c                                                                       
c        
c   THE COEFFICIENTS AND NODES FOR THE LEGENDRE QUADRATURE USED IN QPROB
c AND    
c   WPROB WERE CALCULATED USING THE ALGORITHMS FOUND IN:                
c        
c   STROUD, A. H. AND SECREST, D.  GAUSSIAN QUADRATURE FORMULAS.  ENGLEW
cOOD     
c   CLIFFS, NEW JERSEY:  PRENTICE-HALL, INC, 1966.                      
c        
c   ALL VALUES MATCHED THE TABLES (PROVIDED IN SAME REFERENCE) TO 30    
c        
c   SIGNIFICANT DIGITS.                                                 
c        
c                                                                       
c        
c   FORTRAN FUNCTIONS ERF, ERFC, QEXP, LOG, AND QSQRT                
c        
c   HAVE MAXIMUM RELATIVE ERROR:                                        
c        
c   MAX(CALC(X) - TRUE(X)) / TRUE(X))                                   
c        
c   OF 8 * 10 ** -33.                                                   
c        
c                                                                       
c        
c   F(X) = .5 + ERF(X / SQRT(2)) / 2      FOR X > 0                    
c        
c   F(X) = ERFC( -X / SQRT(2)) / 2        FOR X < 0                    
c        
c   WHERE F(X) IS STANDARD NORMAL C. D. F.                              
c        
c                                                                       
c        
c  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cXXXXXXX 
c                                                                       
c        
c                                                                       
c        
c   IF DEGREES OF FREEDOM LARGE, APPROXIMATE INTEGRAL WITH RANGE DISTRIB
cUTION.  
      data r2 / 0.693147180559945309417232121458d0 /
      data xlegq / 0.989400934991649932596154173450d+00, 
     &0.944575023073232576077988415535d+00, 
     &0.865631202387831743880467897712d+00, 
     &0.755404408355003033895101194847d+00, 
     &0.617876244402643748446671764049d+00, 
     &0.458016777657227386342419442984d+00, 
     &0.281603550779258913230460501460d+00, 
     &0.950125098376374401853193354250d-01 /
      data alegq / 0.271524594117540948517805724560d-01, 
     &0.622535239386478928628438369944d-01, 
     &0.951585116824927848099251076022d-01, 
     &0.124628971255533872052476282192d+00, 
     &0.149595988816576732081501730547d+00, 
     &0.169156519395002538189312079030d+00, 
     &0.182603415044923588866763667969d+00, 
     &0.189450610455068496285396723208d+00 /
      ir(1) = 0
      ir(2) = 0
      if (df .gt. dlarg) then
      qprob = wprob(q,rr,cc,it)
      if (it .eq. 1) ir(1) = 1
      goto 900
c                                                                       
c        
c   CALCULATE LEADING CONSTANT                                          
c        
c   QLGAMA IS USER-WRITTEN FUNCTION.                                    
c        
      end if
      f2 = df * 0.5d0
      f2lf = ((f2 * dlog(df)) - (df * r2)) - qlgama(f2)
      f21 = f2 - 1.0
c                                                                       
c        
c   INTEGRAL IS DIVIDED INTO UNIT, HALF-UNIT, QUARTER-UNIT, OR          
c        
c   EIGHTH-UNIT LENGTH INTERVALS DEPENDING ON THE VALUE OF THE          
c        
c   DEGREES OF FREEDOM.                                                 
c        
      ff4 = df * 0.25d0
      if (df .le. dhaf) then
      ulen = ulen1
      else if (df .le. dquar) then
      ulen = ulen2
      else if (df .le. deigh) then
      ulen = ulen3
      else
      ulen = ulen4
      end if
c                                                                       
c        
      f2lf = f2lf + dlog(ulen)
c                                                                       
c        
c   INTEGRATE OVER EACH SUBINTERVAL                                     
c        
      qprob = 0.0
      do 300 i = 1, 50
      otsum = 0.0
c                                                                       
c        
c       LEGENDRE QUADRATURE WITH ORDER = NLEGQ                          
c        
c       NODES (STORED IN XLEGQ) ARE SYMMETRIC AROUND ZERO.              
c        
      twa1 = ((2.d0 * i) - 1.0) * ulen
c                                                                       
c        
      do 200 jj = 1, nlegq
      if (ihalfq .lt. jj) then
      j = jj - ihalfq
      t1 = (f2lf + (f21 * dlog(twa1 + (xlegq(j) * ulen)))) - (((xlegq(j)
     & * ulen) + twa1) * ff4)
      else
      j = jj
      t1 = (f2lf + (f21 * dlog(twa1 - (xlegq(j) * ulen)))) + (((xlegq(j)
     & * ulen) - twa1) * ff4)
c                                                                       
c        
c           IF EXP(T1) < 9E-14, THEN DOESN'T CONTRIBUTE TO INTEGRAL     
c        
      end if
      if (t1 .ge. eps1) then
      if (ihalfq .lt. jj) then
      qsqz = q * dsqrt(((xlegq(j) * ulen) + twa1) * 0.5d0)
      else
      qsqz = q * dsqrt(((- (xlegq(j) * ulen)) + twa1) * 0.5d0)
c                                                                       
c        
c             CALL WPROB TO FIND INTEGRAL OF RANGE PORTION              
c        
      end if
c                                                                       
c        
      wprb = wprob(qsqz,rr,cc,it)
      if (it .eq. 1) ir(1) = 1
      rotsum = (wprb * alegq(j)) * dexp(t1)
      otsum = rotsum + otsum
c                                                                       
c        
      end if
c                                                                       
c        
c       END LEGENDRE INTEGRAL FOR INTERVAL I                            
c        
c                                                                       
c        
c       IF INTEGRAL FOR INTERVAL I < 1E-14, THEN STOP.  HOWEVER,        
c        
c       IN ORDER TO AVOID SMALL AREA UNDER LEFT TAIL, AT LEAST          
c        
c       1 / ULEN INTERVALS ARE CALCULATED.                              
c        
  200 continue
      if (((i * ulen) .ge. 1.0) .and. (otsum .le. eps2)) goto 400
c                                                                       
c        
c   END OF INTERVAL I                                                   
c        
c                                                                       
c        
  300 qprob = qprob + otsum
  400 if (qprob .gt. eps) ir(2) = 1
      if ((qprob .gt. 1.0) .and. (qprob .le. eps)) qprob = 1.0
  900 return 
c                                                                       
c        
      end
c                                                                       
c        
c   THIS FUNCTION SUBPROGRAM CALCULATES INTEGRAL OF HARTLEY'S FORM      
c        
c   OF THE RANGE.                                                       
c        
c                                                                       
c        
c   W = VALUE OF RANGE                                                  
c        
c   RR = NO. OF ROWS OR GROUPS                                          
c        
c   CC = NO. OF COLUMNS OR TREATMENTS                                   
c        
c   IR = ERROR FLAG = 1 IF WPROB PROBABILITY > 1                        
c        
c   WPROB = RETURNED PROBABILITY INTEGRAL FROM (0, W)                   
c        
c   PROGRAM WILL NOT TERMINATE IF IR IS RAISED.                         
c        
c                                                                       
c        
      function wprob(w, rr, cc, ir)
c                                                                       
c        
c   BB = UPPER LIMIT OF LEGENDRE INTEGRATION                            
c        
c   EPS = MAXIMUM ACCEPTABLE VALUE OF INTEGRAL                          
c        
c   NLEG = ORDER OF LEGENDRE QUADRATURE                                 
c        
c   IHALF = INT ((NLEG + 1) / 2)                                        
c        
c   WLAR = VALUE OF RANGE ABOVE WHICH WINCR1 INTERVALS ARE USED TO      
c        
c   CALCULATE SECOND PART OF INTEGRAL,                                  
c        
c   ELSE WINCR2 INTERVALS ARE USED.                                     
c        
c   EPS1, EPS2, EPS3 = VALUES WHICH ARE USED AS CUTOFFS FOR TERMINATING 
c        
c   OR MODIFYING A CALCULATION.                                         
c        
      implicit double precision (o-z, a-h)
c                                                                       
c        
      parameter (bb = 8.0d0, eps = 1.0d0, nleg = 12, ihalf = 6, wlar = 
     &3.0d0, wincr1 = 2.0d0, wincr2 = 3.0d0, eps1 = -30.0d0, eps2 = 
     &-50.0d0, eps3 = 60.0d0)
c                                                                       
c        
c   SQ2PII = 1 / SQRT(2 * PI);  FROM ABRAMOWITZ & STEGUN, P. 3.         
c        
c   QSQR2 = SQRT(2)                                                     
c        
c   XLEG = LEGENDRE 12-POINT NODES                                      
c        
c   ALEG = LEGENDRE 12-POINT COEFFICIENTS                               
c        
      dimension xleg(ihalf), aleg(ihalf)
c                                                                       
c        
c  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cXXXXXXXX
c                                                                       
c        
      data sq2pii / 0.3989422804014326779399461d0 /
      data qsqr2 / 1.41421356237309504880168872421d0 /
      data xleg / 0.981560634246719250690549090149d+00, 
     &0.904117256370474856678465866119d+00, 
     &0.769902674194304687036893833213d+00, 
     &0.587317954286617447296702418941d+00, 
     &0.367831498998180193752691536644d+00, 
     &0.125233408511468915472441369464d+00 /
      data aleg / 0.471753363865118271946159614850d-01, 
     &0.106939325995318430960254718194d+00, 
     &0.160078328543346226334652529543d+00, 
     &0.203167426723065921749064455810d+00, 
     &0.233492536538354808760849898925d+00, 
     &0.249147045813402785000562436043d+00 /
      ir = 0
      qsqz = w * 0.5d0
c   IF W >= 16 THEN INTEGRAL LOWER BOUND (OCCURS FOR C=20) IS 0.99999999
c999995  
c   SO RETURN A VALUE OF 1.                                             
c        
      wprob = 1.0
c                                                                       
c        
c   FIND (F(W/2) - 1) ** CC (FIRST TERM IN INTEGRAL OF HARTLEY'S FORM). 
c        
      if (qsqz .ge. bb) goto 900
c   IF WPROB ** CC < 2E-22 THEN SET WPROB = 0                           
c        
      wprob = derf(qsqz / qsqr2)
      if (wprob .ge. dexp(eps2 / cc)) then
      wprob = wprob ** cc
      else
      wprob = 0.0
c                                                                       
c        
c   IF W IS LARGE THEN SECOND COMPONENT OF INTEGRAL IS SMALL, SO FEWER  
c        
c   INTERVALS ARE NEEDED.                                               
c        
      end if
      if (w .gt. wlar) then
      wincr = wincr1
      else
      wincr = wincr2
c                                                                       
c        
c   FIND INTEGRAL OF SECOND TERM OF HARTLEY'S FORM FOR INTEGRAL         
c        
c   OF THE RANGE FOR EQUAL-LENGTH INTERVALS USING LEGENDRE QUADRATURE.  
c        
c   LIMITS OF INTEGRATION ARE FROM (W/2, 8).                            
c        
c   TWO  OR THREE EQUAL-LENGTH INTERVALS ARE USED.                      
c        
c                                                                       
c        
c   BLB AND BUB ARE LOWER AND UPPER LIMITS OF INTEGRATION.              
c        
      end if
      blb = qsqz
      binc = (bb - qsqz) / wincr
      bub = blb + binc
      einsum = 0.0
c                                                                       
c        
c   INTEGRATE OVER EACH INTERVAL                                        
c        
      cc1 = cc - 1.0
      do 300 wi = 1, wincr
      elsum = 0.0
      a = 0.5d0 * (bub + blb)
c                                                                       
c        
c       LEGENDRE QUADRATURE WITH ORDER = NLEG                           
c        
      b = 0.5d0 * (bub - blb)
c                                                                       
c        
      do 100 jj = 1, nleg
      if (ihalf .lt. jj) then
      j = (nleg - jj) + 1
      xx = xleg(j)
      else
      j = jj
      xx = - xleg(j)
c                                                                       
c        
      end if
      c = b * xx
      ac = a + c
c           IF EXP(-QEXPO/2) < 9E-14, THEN DOESN'T CONTRIBUTE TO INTEGRA
cL       
      qexpo = ac * ac
c                                                                       
c        
      if (qexpo .gt. eps3) goto 200
      if (ac .gt. 0.0) then
      pplus = 1.0 + derf(ac / qsqr2)
      else
      pplus = derfc(- (ac / qsqr2))
c                                                                       
c        
      end if
      if (ac .gt. w) then
      pminus = 1.0 + derf((ac / qsqr2) - (w / qsqr2))
      else
      pminus = derfc((w / qsqr2) - (ac / qsqr2))
c                                                                       
c        
      end if
c                                                                       
c        
c           IF RINSUM ** (CC-1) < 9E-14, THEN DOESN'T CONTRIBUTE TO INTE
cGRAL    
      rinsum = (pplus * 0.5d0) - (pminus * 0.5d0)
      if (rinsum .ge. dexp(eps1 / cc1)) then
      rinsum = (aleg(j) * dexp(- (0.5d0 * qexpo))) * (rinsum ** cc1)
      elsum = elsum + rinsum
c                                                                       
c        
      end if
c       END LEGENDRE QUADRATURE                                         
c        
c                                                                       
c        
  100 continue
  200 elsum = (((2.0 * b) * cc) * sq2pii) * elsum
      einsum = einsum + elsum
      blb = bub
c                                                                       
c        
      bub = bub + binc
c   END INTEGRATION OF SECOND TERM                                      
c        
c                                                                       
c        
  300 continue
c                                                                       
c        
c   IF WPROB ** RR < 9E-14, THEN RETURN 0.0                             
c        
      wprob = einsum + wprob
      if (wprob .le. dexp(eps1 / rr)) then
      wprob = 0.0
      goto 900
c                                                                       
c        
      end if
      wprob = wprob ** rr
      if (wprob .gt. eps) ir = 1
      if ((wprob .gt. 1.0) .and. (wprob .le. eps)) wprob = 1.0
  900 return 
c                                                                       
c        
      end
c                                                                       
c        
c   THIS FUNCTION SUBPROGRAM CALCULATES THE LOG GAMMA OF X.             
c        
c   FOR X SMALL, CALCULATION IS EXACT; ELSE AN APPROXIMATION IS USED.   
c        
c                                                                       
c        
      function qlgama(x)
c                                                                       
c        
c   XX IS BOUNDARY FOR EXACT VS APPROXIMATE CALCULATION.                
c        
      implicit double precision (o-z, a-h)
c                                                                       
c        
c   SQRPI= SQRT(PI) & GAMMA(1/2);   FROM ABRAMOWITZ & STEGUN , P. 3.    
c        
c   GSQR2P = 0.5 * LOG(2 * PI) ;  FROM ABRAMOWITZ & STEGUN, P. 3.       
c        
      parameter (xx = 21.0d0)
c                                                                       
c        
c  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cXXXXXXXX
c                                                                       
c        
c   IF X < XX & X IS AN INTEGER THEN GAMMA(X) = (X - 1) FACTORIAL       
c        
      data sqrpi / 1.772453850905516027298167d0 /
      data gsqr2p / 0.9189385332046727417803297d0 /
      data o12 / 0.833333333333333333333333333333d-01 /
      data o30 / 30.d0 /
      data o35 / 3.5d0 /
      data o13 / 1.33333333333333333333333333333d0 /
      data o71 / 0.707142857142857142857142857143d0 /
      frac = x - int(x)
      if ((x .lt. xx) .and. (frac .eq. 0.0)) then
      qlgama = 1.0
      y = x - 1.0
      do 100 a = 1.0, y
      qlgama = qlgama * a
  100 continue
c                                                                       
c        
c   IF X < XX & FRAC = 0.5 THEN GAMMA(X) IS                             
c        
c   FORMULA 6.1.12 IN ABRAMOWITZ & STEGUN, P. 255.                      
c        
      qlgama = dlog(qlgama)
      else if (x .lt. xx) then
      qlgama = 1.0
      y = (2.0 * x) - 1.0
      do 200 a = 1.0, y, 2.0
      qlgama = (qlgama * a) * 0.5d0
  200 continue
c                                                                       
c        
c   IF X >= XX THEN LOG GAMMA(X) IS FORMULA 6.1.41 IN                   
c        
c    ABRAMOWITZ & STEGUN, P. 257.                                       
c        
c   THE MAXIMUM RELATIVE ERROR ASSOCIATED WITH THIS FORMULA IS 1.3E-19  
c        
c   WHICH OCCURS AT X = 21.  ERROR DECREASES AS X INCREASES.  FOR EXAMPL
cE,      
c   THE ERROR FOR X = 50 IS 2.5E-24.                                    
c        
      qlgama = dlog(qlgama * sqrpi)
      else
      y = 1.0 / (x * x)
      qlgama = (o12 * (1.0 - ((y / o30) * (1.0 - ((y / o35) * (1.0 - ((y
     & / o13) * (1.0 - (y / o71))))))))) / x
      qlgama = ((qlgama + gsqr2p) + ((x - 0.5d0) * dlog(x))) - x
c                                                                       
c        
      end if
  900 return 
c                                                                       
c        
      end
      FUNCTION QINV (P, C, V)
C                                                                               
C   THIS FUNCTION FINDS PERCENTAGE POINT OF THE STUDENTIZED RANGE               
C   WHICH IS USED AS INITIAL ESTIMATE FOR THE SECANT METHOD.                    
C   FUNCTION IS ADAPTED FROM PORTION OF ALGORITHM AS 70                         
C   FROM APPLIED STATISTICS (1974) ,VOL. 23, NO. 1                              
C   BY ODEH, R. E. AND EVANS, J. O.                                             
C                                                                               
C   P = PERCENTAGE POINT                                                        
C   C = NO. OF COLUMNS OR TREATMENTS                                            
C   V = DEGREES OF FREEDOM                                                      
C   QINV = RETURNED INITIAL ESTIMATE                                            
C                                                                               
      implicit double precision (o-z, a-h)

C                                                                               
C   VMAX IS CUTOFF ABOVE WHICH DEGREES OF FREEDOM IS TREATED AS INFINITY.       
      PARAMETER (VMAX = 120.0)                                                  
C                                                                               
      DATA P0 /0.322232421088/, Q0 /0.993484626060E-01/,
     # P1 /-1.0/, Q1 /0.588581570495/,
     # P2 /-0.342242088547/, Q2 /0.531103462366/,
     # P3 /-0.204231210125/, Q3 /0.103537752850/,
     # P4 /-0.453642210148E-04/, Q4 /0.38560700634E-02/
     # C1 /0.8832/, C2 /0.2368/, C3 /1.214/,
     # C4 /1.208/, C5 /1.4142/
C                                                                               
C  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      PS = 0.5 - 0.5 * P
      YI = SQRT (LOG (1.0 / (PS * PS)))
      T = YI + (((( YI * P4 + P3) * YI + P2) * YI + P1) * YI + P0)
     # / (((( YI * Q4 + Q3) * YI + Q2) * YI + Q1) * YI + Q0)
      IF (V .LT. VMAX) T = T + (T * T * T + T) / V / 4.0
      Q = C1 - C2 * T
      IF (V .LT. VMAX) Q = Q - C3 / V + C4 * T / V
      QINV = T * (Q * LOG (C - 1.0) + C5)
  900 RETURN
      END

      SUBROUTINE cv2(p,rr,cc,df,ir,q)
      DOUBLE PRECISION p,rr,cc,df,q,cv
      INTEGER ir(3)
      q = cv(p,rr,cc,df,ir)
      END

      SUBROUTINE qprob2(q, rr, cc, df, ir,p)
      DOUBLE PRECISION p,rr,cc,df,q,qprob
      integer ir(2)
      p = qprob(q,rr,cc,df,ir)
      END
