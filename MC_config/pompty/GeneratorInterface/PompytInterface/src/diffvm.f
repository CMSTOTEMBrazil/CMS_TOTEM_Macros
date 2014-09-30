       SUBROUTINE PRODIFFVM

c      IMPLICIT NONE
*
*  PROTON DISSOCIATION A LA DIFFVM
*
c#include "rglujets.inc"
c#include "rgpdiss.inc"
c#include "rgpart.inc"

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
ccc       INTEGER PYK,PYCHGE,PYCOMP
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)


c      DOUBLE PRECISION DMXP       
c      INTEGER IBP,IDIFP,J,ICHECK,IPROT,I
c      INTEGER IPOM,KPART
c      REAL RFRAC
c      Double Precision draprn
c      EXTERNAL draprn
	
cc      DATA ICHECK /0/


      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2  

      DOUBLE PRECISION CMGP(4)

      integer pdvm
      common/prodvm/pdvm

cy in ep CMS FRAME

* Type of dissociation
      IFRAGP=1

* beam proton
c      DO 5 I=1,2
c         IF(IABS(K(I,2)).EQ.2212) IBP = I
c 5    CONTINUE

      IBP=1
 

c      IF (ICHECK.EQ.1) THEN
c         write(6,*)'PRODIFFVM: beam proton ibp',IBP
c      ENDIF
 
      DO 10 I=3,N
cc         IF(IABS(K(I,2)).EQ.2212) IPROT = I
         IF(IABS(K(I,2)).EQ.990) IPOM = I
 10   CONTINUE

      iprot=3

cc      IF(IPROT.LT.0) RETURN
      IDIFP=IPROT
      KPART = K(IPROT,2)  

      DMXP=P(IDIFP,5)


c      write(6,*) 'prodiffvm MY = ', DMXP, sqrt(MY2),ifragp

* Determine the state
      IF (IFRAGP .EQ. 1 .OR. IFRAGP .EQ. 2) THEN
         IF (DMXP .LT. 1.48) THEN
            K(IDIFP,2) = SIGN (12212,KPART)
            K(IDIFP,1)=11
            npdc1 = idifp
            npdc2 = 0
         ELSE IF (DMXP .LT. 1.60) THEN
            K(IDIFP,2) = SIGN (2124,KPART)
            K(IDIFP,1)=11
            npdc1 = idifp
            npdc2 = 0
         ELSE IF (DMXP .LT. 1.90) THEN
            RFRAC = PYR(0)
            IF (RFRAC .LT. 0.5) THEN
               K(IDIFP,2)   = SIGN (12216, KPART)
               K(IDIFP,1)=11
               npdc1 = idifp
               npdc2 = 0
            ELSE IF (RFRAC .LT. 0.83) THEN
               K(IDIFP,2) = SIGN (22124, KPART)
               K(IDIFP,1)=11
                npdc1 = idifp
                npdc2 = 0
            ELSE
               K(IDIFP,2) = SIGN (42212, KPART)
               K(IDIFP,1)=11
               npdc1 = idifp
               npdc2 = 0
            END IF
         ELSE
* make proton a diffracted p_diffr+ state 
c            K(IDIFP,2) = SIGN (2210, KPART)
            K(IDIFP,2) = SIGN (9902210, KPART)
            pdvm=1
         END IF
      ELSE
         K(IDIFP,2) = SIGN (ABS (IFRAGP), KPART)
         write(6,*) 'PRODIFFVM:   incorrect value of IFRAGP'
      END IF

      DO 20 I=1,4
         CMGP(I) =P(2,I)+P(3,I)
   20 CONTINUE   

* THE FOLLOWING IS FROM GENEVT IN DIFFVM
*
* If necessary, perform fragmentation of diffractive proton system
*
      IF (ABS (K (IDIFP,2)) .EQ. 2212) THEN
         write(6,*) 'PRODIFFVM:Should not be true:Elastic proton!'
      ELSE IF (ABS (K (IDIFP,2)) .EQ. 9902210) THEN
         IF (IFRAGP .EQ. 1) THEN   
c	      write(6,*) ' prodiffvm: no frag. for FRAGPJ '         
            CALL FRAGPJ (IBP,IDIFP)
         ELSE
c	      write(6,*) ' prodiffvm: no frag. for FRAGPX '         
ccc            CALL FRAGPX (IBP,IDIFP)
         END IF
      endif

* do the decay of the resonances later on in pdisdcvm.f


      RETURN
      END
c ---------------------------------------------------------------------
c     FRAGPJ
c ---------------------------------------------------------------------
      SUBROUTINE FRAGPJ (IBP, IDIFP)

c      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Take hadronic mass at p vertex and let it decay
*     --------    
*
*     Input:      IBP:    Number of beam proton in PYJETS
*     ------      IDIFP:  Number of diffractive state in PYJETS
*
*     Called by:  PRODIFFVM
*     ----------
*
*     Author:     Benno List, 17.1.94
*     -------
*
*     Changed by: Paul Thompson
*     -----------
*
*     Calling:    SPLITP
*     --------    
*
************************************************************************
c #include "rglujets.inc"
c #include "rgpdiss.inc"

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
ccc       INTEGER PYK,PYCHGE,PYCOMP
C
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)


********************************************************
      DOUBLE PRECISION DMP,DMPI0
      PARAMETER (DMP = 0.938270D0,DMPI0 = 0.139568D0)
      DOUBLE PRECISION PQQ (5), PQ (5), PPOM (5)
      INTEGER IBP, IDIFP
********************************************************
      DOUBLE PRECISION PPCMS8(5)
c      Double Precision draprn
c      EXTERNAL draprn
c      INTEGER IPOM,IBEAMP,IQ,IQQ,I,IERR,ICHECK
c      REAL R
c      DATA ICHECK /0/

      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2

      DOUBLE PRECISION la_pg(5),la_prem(5)
	DOUBLE PRECISION ptes(5)

       common/pomcut1/ pomcut(20) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




*****
cc      INTEGER IARR12 (2) /1, 2/
      
*****


      Integer I_glu
      I_glu=0
      if(PYR(0).lt.pomcut(10)) I_glu=1

      pomcut(9)=I_glu   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cccc      write(6,*) 'POMCUT(10)=', pomcut(10)

c	Data I_glu/0/  !  old diffvm treatment 
c	Data I_glu/1/  !  new treatment, adding a gluon, done by Leszek Adamczyk



********************************************************

C     PRINT *, '--- FRAGPJ ---'
C     PRINT *, 'IBP   = ', IBP
C     PRINT *, 'IDIFP = ', IDIFP

c      write(6,*) ' fragjp I_glu :',I_glu

      IF ( P(IDIFP,5) .LT. DMP+DMPI0) THEN
         PRINT *, '### FRAGPJ error: not enough energy! '
         PRINT *, 'IBP:   ', IBP
         PRINT *, 'IDIFP: ', IDIFP
         PRINT *, 'P(IDIFP,5): ', P(IDIFP,5)
         PRINT *, 'DMP + DMPI0:       ', DMP+DMPI0
c        CALL ERRLOG (40, 'F: FRAGPJ: Not enough energy!')
      END IF
*
* Split proton into q - qq system: choose q and qq codes
* The pomeron couples with the same strength to u and d quarks.
* ud_1 and uu_1 are spin triplets, ud_0 is a spin singulett.
* Therefore q-qq is 1/7 u-ud_0, 3/7 u-ud_1, and 3/7 d-uu_1.
*
      IBEAMP=K(IBP,2)
      IF (ICHECK.EQ.1) write(6,*) 'SPLITP: ibeamp',ibeamp
      R = 7.0*PYR(0)
      IF (R .LE. 1.0) THEN
         IQ  = SIGN (2, IBEAMP)       !u
         IQQ = SIGN (2101, IBEAMP)    !ud spin 0
      ELSE IF (R .LE. 4.0) THEN
         IQ  = SIGN (2, IBEAMP)       !u
         IQQ = SIGN (2103, IBEAMP)    !ud spin 1
      ELSE
         IQ  = SIGN (1, IBEAMP)       !d
         IQQ = SIGN (2203, IBEAMP)    !uu spin 1
      END IF
      
      IF (ICHECK.EQ.1) THEN
      write (6,*) 'FRAGPJ IQ:  ', IQ,  ', mass: ', PYMASS (IQ)
      write (6,*) 'FRAGPJ IQQ: ', IQQ, ', mass: ', PYMASS (IQQ)
      ENDIF

      DO I=3,N
         IF(IABS(K(I,2)).EQ.990) IPOM = I
      enddo
* Fill 5 vector of the beam proton and pomeron 
      do i=1,5
         PPCMS8(I)= P(IBP,I)
chju         PPOM (I) = P(IPOM,I)
* To avoid numerical problems, calculate pomeron momentum from
* proton side
	   PPom(I) = P(IDIFP,I) - P(IBP,I)
         IF (ICHECK.EQ.1) THEN
         write(6,*)'I,PPOM,P(IPOM),PROTON',I,PPOM(I),P(IPOM,I),P(IBP,I)
         ENDIF
      enddo
      PPOM (5) = P(IPOM,5)

      IF (ICHECK.EQ.1) THEN
         write(6,*) 'FRAGPJ: PPOM(5)',PPOM(5)
         write(6,*) 'FRAGPJ: CALL SPLITP'
      ENDIF

      If(I_glu.eq.0) then 
          CALL SPLITP (PPCMS8, PPOM, PQ, PQQ, IERR)
      else
cnew
          CALL SPLITP (PPCMS8, PPOM, la_pg, la_prem, IERR)
          do i=1,5
           pq(i) = 1./3. * la_prem(i)
           pqq(i) = 2./3. * la_prem(i)
	     ptes(i)=pq(i)+pqq(i)+la_pg(i)
          enddo
c	    write(6,*) 'ptest ',ptes
	endif
cnew
      IF (ICHECK.EQ.1) THEN
         write (6,*) 'N before',N
         write (6,*) 'SPLITP',IERR
      ENDIF

* Proton  decay

 
c      write(*,*) 'pylist 1',Ierr
c      call pylist(1)

c    write(6,*) 'PRODIFFVM: idifp, ipom, N', idifp,ipom,n

      K(IDIFP,1)=11
      N=N+1
      npdc1 = N
      DO I=1,5
         P(N,I)=PQQ (I)
      ENDDO
      K(N,1)=12
cy      K(N,1)=11  ! diquark 
      K(N,2)=IQQ
      K(N,3)=IDIFP

c      IF (ICHECK.EQ.1) write(6,*)'in between N',N
      If(I_glu.eq.1) then
cnew
         N=N+1
         DO I=1,5
            P(N,I)=la_Pg(I)
         ENDDO
         K(N,1)=12
         K(N,2)=21
         K(N,3)=IDIFP
cnew	
      Endif
      N=N+1
      npdc2 = N
      DO I=1,5
         P(N,I)=PQ(I)
      ENDDO
      K(N,1)=11
cy      K(N,1)=11 ! quark
      K(N,2)=IQ
      K(N,3)=IDIFP



      IF (ICHECK.EQ.1) THEN
         write(6,*)'after adding 2'
         write(6,*)'N',N
         write(6,*) 'Fragpj: end pylist'
         CALL pylist (1)
      ENDIF

c      write(*,*) 'FRAGPJ'
c      call pylist(1)
   
      RETURN
      END
c ----------------------------------------
c     SPLITP
c -----------------------------------------
      SUBROUTINE SPLITP (P, Q, PQ, PREMN, IERR)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:  A quark in a hadron with 5-momentum P is struck
*     --------  by a boson with 5-momentum Q;
*               Q is to be negative for spacelike particles.
*               The hadron splits into
*               a remnant with 5-momentum (1-x)*P(5) and
*               a quark with mass P (5) and 3-momentum x*P(1:3)+Q(1:3).
*               x is determined so that 4-vector conservation holds.
*
*     Input:    P:         5-momentum of the hadron
*     ------    Q:         5-momentum of the boson
*
*     Output:   PQ:        5-vector of the quark
*     -------   PREMN:     5-vector of the remnant
*               IERR:      0 if everything is O.K.,
*                         -1 if x<0 or x>1.
*                          1 if no solution was found.
*
*     Author:   Benno List, 20.1.94
*     -------
*
************************************************************************
      Double Precision P (5), Q (5), PQ (5), PREMN (5)
      INTEGER IERR,I
      Double Precision  X, D, P2, Q2
      
      P2= P (5)*ABS (P (5))
      Q2= Q (5)*ABS (Q (5))
      D = (P(4)+Q(4))**2-(P(1)+Q(1))**2-(P(2)+Q(2))**2-(P(3)+Q(3))**2
     +     - P2 - Q2
      IF (D .EQ. 0D0) THEN
         IERR = 1
         PRINT *, '### SPLITP error: D = 0!'
         PRINT *, 'P = (', P(1),',',P(2),',',P(3),',',P(4),'), ',P(5)
         PRINT *, 'Q = (', Q(1),',',Q(2),',',Q(3),',',Q(4),'), ',Q(5)
      ELSE
         X = -Q2 / D
         IF (X .LT. 0D0 .OR. X .GT. 1D0) THEN
            IERR = -1
            write(6,*) 'SPLITP: x<0 or x>1 X,D',X,D
         ELSE
            IERR = 0
         END IF
         DO 10, I = 1, 4
            PREMN (I) = (1D0 - X)*P (I)
 10      CONTINUE
         PREMN (5) = ABS (1D0 - X)*P (5)
         PQ (5) = ABS (X)*P (5)
         PQ (4) = X**2 * P2
         DO 20, I = 1, 3
            PQ (I) = X*P (I) + Q (I)
            PQ (4) = PQ (4) + PQ (I)**2
 20      CONTINUE
         IF (PQ (4) .LT. 0D0) THEN
            IERR = 1
            write(6,*)'SPLITP: No solution found'
         END IF
         PQ (4) = SQRT (ABS (PQ (4)))

c         PRINT *, '   SPLITP: X =             ', X
c         PRINT *, '   SPLITP: (P+Q)**2 =      ', (P (4) + Q (4))**2
c     +        -(P (1)+Q (1))**2-(P (2)+Q (2))**2-(P (3)+Q (3))**2
c         PRINT *, '   SPLITP: (PQ+PREMN)**2 = ',(PQ (4) + PREMN (4))**2
c     +        -(PQ(1)+PREMN(1))**2-(PQ(2)+PREMN(2))**2
c     +        -(PQ(3)+PREMN(3))**2
c         PRINT *, '   SPLITP: PQ (5) =        ', PQ (5)
c         PRINT *, '   SPLITP: PREMN (5) =     ', PREMN (5)
      END IF

      RETURN
      END
c -------------------------------------------------------------------
c
c ---------------------------------------------------------------------
      SUBROUTINE PDISDCVM

c      IMPLICIT NONE
*
*  PROTON DISSOCIATION A LA DIFFVM
*
c#include "rglujets.inc"
c#include "rgherup.inc"

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
ccc       INTEGER PYK,PYCHGE,PYCOMP
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)


c      DOUBLE PRECISION DMXP       
c      INTEGER IFRAGP,IBP,IDIFP,J,ICHECK,IPROT,I
c      INTEGER IPOM,KPART
c      REAL RFRAC
c     Double Precision draprn
c     EXTERNAL draprn
c     DATA ICHECK /1/

      DOUBLE PRECISION SSS,CM,DBCMS
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)


      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2  

cc      DOUBLE PRECISION CMGP(4)
      
cc      DOUBLE PRECISION beta(3)

      double precision my2
      integer pdvm
      common/mykine/ MY2
      common/prodvm/pdvm

cy in ep CMS FRAME

* Type of dissociation
      IFRAGP=1
      pdvm=0
* beam proton

      if(npdc1.eq.0.and.npdc2.eq.0) write(6,*) ' pdisdcvm ',npdc1,npdc2
      if(npdc2.ne.0) then
         K(npdc1,1) = 2
	   if(npdc1+1.ne.npdc2) then
	     K(npdc1+1,1) = 2
	   endif
         K(npdc2,1) = 1

ccc         if(ILHA.eq.0) call pyexec
         call pyexec

      else

* decay the resonances
c        K(npdc1,1) = 1

c        write(*,*) 'pylist in pdiscvm '
c        call pylist(1)
 
        IF (ABS (K (npdc1,2)) .EQ. 12212 .OR.
     +        ABS (K (npdc1,2)) .EQ.  2124 .OR.
     +        ABS (K (npdc1,2)) .EQ. 22212 .OR.
     +        ABS (K (npdc1,2)) .EQ. 32212 .OR.
     +        ABS (K (npdc1,2)) .EQ.  2216 .OR.
     +        ABS (K (npdc1,2)) .EQ. 12216 .OR.
     +        ABS (K (npdc1,2)) .EQ. 22124 .OR.
     +        ABS (K (npdc1,2)) .EQ. 42212 .OR.
     +        ABS (K (npdc1,2)) .EQ. 32124) THEN
c           write(6,*)  'K(IDIFP,2)=', K(IDIFP,2)
c                 write(6,*) 'pdisdcvm: no decnst '
              K(npdc1,1) = 1
              CALL DECNST (npdc1)
         ELSE IF (ABS (K (IDIFP,2)) .EQ. 2214) THEN
              K(npdc1,1) = 1
                 write(6,*) 'pdisdcvm: no decdel  DELTA '
                 stop
c              CALL DECDEL (npdc1)
         ELSE
         END IF
 
      endif

      RETURN
      END
c --------------------------------------------------
c
c ---------------------------------------------------
	subroutine pdismass(ipdg,xm,width)
	Implicit none
	double precision xm,width
	Integer ipdg,i,j,ncomax
	Integer ip
	Real am,ch,tl,wi
	Logical found
      PARAMETER (NCOMAX=55)
      CHARACTER*12 NA(NCOMAX)
      DIMENSION IP(NCOMAX),AM(NCOMAX),CH(NCOMAX),TL(NCOMAX),WI(NCOMAX)
      DATA (IP(J),NA(J),AM(J),CH(J),TL(J),WI(J),J=  1, ncomax)/
     +      22,'gamma       ',0.0000E+00, 0.00, 0.1000E+16, 0.0000E+00,
     +     111,'pi0         ',0.1350E+00, 0.00, 0.8439E-16, 0.7800E-08,
     +     113,'rho0        ',0.7700E+00, 0.00, 0.4368E-23, 0.1507E+00,
     +     211,'pi+         ',0.1396E+00, 1.00, 0.2603E-07, 0.2528E-16,
     +    -211,'pi-         ',0.1396E+00,-1.00, 0.2603E-07, 0.2528E-16,
     +     213,'rho+        ',0.7700E+00, 1.00, 0.4368E-23, 0.1507E+00,
     +    -213,'rho-        ',0.7700E+00,-1.00, 0.4368E-23, 0.1507E+00,
     +     321,'K+          ',0.4937E+00, 1.00, 0.1239E-07, 0.5314E-16,
     +    -321,'K-          ',0.4937E+00,-1.00, 0.1239E-07, 0.5314E-16,
     +    2212,'p+          ',0.9383E+00, 1.00, 0.1000E+16, 0.0000E+00,
     +   -2212,'p-#         ',0.9383E+00,-1.00, 0.1000E+16, 0.0000E+00,
     +    2112,'n0          ',0.9396E+00, 0.00, 0.8866E+03, 0.7424E-27,
     +   -2112,'n0#         ',0.9396E+00, 0.00, 0.8866E+03, 0.7424E-27,
     +    2114,'Delta0      ',0.1232E+01, 0.00, 0.5485E-23, 0.1200E+00,
     +   -2114,'Delta0#     ',0.1232E+01, 0.00, 0.5485E-23, 0.1200E+00,
     +    2214,'Delta+      ',0.1232E+01, 1.00, 0.5485E-23, 0.1200E+00,
     +   -2214,'Delta-#     ',0.1232E+01,-1.00, 0.5485E-23, 0.1200E+00,
     +    2224,'Delta++     ',0.1232E+01, 2.00, 0.5485E-23, 0.1200E+00,
     +   -2224,'Delta--#    ',0.1232E+01, 2.00, 0.5485E-23, 0.1200E+00,
     +    3122,'Lambda0     ',0.1116E+01, 0.00, 0.2632E-09, 0.2501E-14,
     +   -3122,'Lambda0#    ',0.1116E+01, 0.00, 0.2632E-09, 0.2501E-14,
     +   12112,'N(1440)0    ',0.1440E+01, 0.00, 0.1881E-23, 0.3500E+00,
     +  -12112,'N(1440)#    ',0.1440E+01, 0.00, 0.1881E-23, 0.3500E+00,
     +   12212,'N(1440)+    ',0.1440E+01, 1.00, 0.1881E-23, 0.3500E+00,
     +  -12212,'N(1440)-#   ',0.1440E+01,-1.00, 0.1881E-23, 0.3500E+00,
     +    2124,'N(1520)+    ',0.1520E+01, 1.00, 0.5485E-23, 0.1200E+00,
     +   -2124,'N(1520)-#   ',0.1520E+01,-1.00, 0.5485E-23, 0.1200E+00,
     +   22112,'N(1535)0    ',0.1535E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +  -22112,'N(1535)0#   ',0.1535E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +   22212,'N(1535)+    ',0.1535E+01, 1.00, 0.4388E-23, 0.1500E+00,
     +  -22212,'N(1535)-#   ',0.1535E+01,-1.00, 0.4388E-23, 0.1500E+00,
     +   32112,'N(1650)0    ',0.1650E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +  -32112,'N(1650)0#   ',0.1650E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +   32212,'N(1650)+    ',0.1650E+01, 1.00, 0.4388E-23, 0.1500E+00,
     +  -32212,'N(1650)-#   ',0.1650E+01,-1.00, 0.4388E-23, 0.1500E+00,
     +    2116,'N(1675)0    ',0.1675E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +   -2116,'N(1675)0#   ',0.1675E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +    2216,'N(1675)+    ',0.1675E+01, 1.00, 0.4388E-23, 0.1500E+00,
     +   -2216,'N(1675)-#   ',0.1675E+01,-1.00, 0.4388E-23, 0.1500E+00,
     +   12116,'N(1680)0    ',0.1680E+01, 0.00, 0.5063E-23, 0.1300E+00,
     +  -12116,'N(1680)0#   ',0.1680E+01, 0.00, 0.5063E-23, 0.1300E+00,
     +   12216,'N(1680)+    ',0.1680E+01, 1.00, 0.5063E-23, 0.1300E+00,
     +  -12216,'N(1680)-#   ',0.1680E+01,-1.00, 0.5063E-23, 0.1300E+00,
     +   21214,'N(1700)0    ',0.1700E+01, 0.00, 0.6582E-23, 0.1000E+00,
     +  -21214,'N(1700)0#   ',0.1700E+01, 0.00, 0.6582E-23, 0.1000E+00,
     +   22124,'N(1700)+    ',0.1700E+01, 1.00, 0.6582E-23, 0.1000E+00,
     +  -22124,'N(1700)-#   ',0.1700E+01,-1.00, 0.6582E-23, 0.1000E+00,
     +   42112,'N(1710)0    ',0.1710E+01, 0.00, 0.6582E-23, 0.1000E+00,
     +  -42112,'N(1710)0#   ',0.1710E+01, 0.00, 0.6582E-23, 0.1000E+00,
     +   42212,'N(1710)+    ',0.1710E+01, 1.00, 0.6582E-23, 0.1000E+00,
     +  -42212,'N(1710)-#   ',0.1710E+01,-1.00, 0.6582E-23, 0.1000E+00, 
     +   31214,'N(1720)0    ',0.1720E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +  -31214,'N(1720)0#   ',0.1720E+01, 0.00, 0.4388E-23, 0.1500E+00,
     +   32124,'N(1720)+    ',0.1720E+01, 1.00, 0.4388E-23, 0.1500E+00,
     +  -32124,'N(1720)-#   ',0.1720E+01,-1.00, 0.4388E-23, 0.1500E+00/


      found = .false.
      Do I=1,NCOMAX
	if(ip(i).eq.ipdg) then
	  found=.true.
	  xm=am(i)
	  width = wi(i)
	endif
	enddo
	
	if(.not.found) then
	  write(6,*) ' pdismass: no proper particle code found ',ipdg
	  stop
	Endif
	
	Return
	end
c ---------------------------------------------------------------------------
c
c ---------------------------------------------------------------------------
      SUBROUTINE DECNST (INSTAR)

c      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Let a N* decay
*     --------    Possible decay modes (also charge conjugates)
*                 N*+ -> p pi0, n pi+
*                 N*+ -> Delta++ pi-, Delta+ pi0, Delta0 pi+
*                 N*+ -> p (pi+ pi-)_I=0, p (pi0 pi0)_I=0
*                 N*+ -> p rho0, n rho+
*                 N*+ -> p eta
*                 N*+ -> Lambda K+
*                 N*0 -> p pi-, n pi0
*                 N*0 -> Delta+ pi-, Delta0 pi0, Delta- pi+
*                 N*0 -> n (pi+ pi-)_I=0, n (pi0 pi0)_I=0
*                 N*0 -> p rho-, n rho0
*                 N*0 -> n eta
*                 N*0 -> Lambda K0S, Lambda K0L
*                 Branching ratios between lines taken from
*                   Rev. Part. Prop. 1994: Phys. Rev. D50, Pt. I, 1173.
*                 Branching ratios between particles in one line
*                 calculated from CLEBSCH-GORDAN coefficients
*                 N* can be:
*                 name     PDG code      name     PDG code
*                 N(1440)+    12212  N(1440)0    12112
*                 N(1520)+     2124  N(1520)0     1214
*                 N(1535)+    22212  N(1535)0    22112
*                 N(1650)+    32212  N(1650)0    32112
*                 N(1675)+     2216  N(1675)0     2116
*                 N(1680)+    12216  N(1680)0    12116
*                 N(1700)+    22124  N(1700)0    21214
*                 N(1710)+    42212  N(1710)0    42112
*                 N(1720)+    32124  N(1720)0    31214
*
*     Input:      INSTAR:  Number of N* in /PYJETS/
*     ------      
*
*     Output:     
*     -------     
*
*     Called by:  PRODIFFVM
*     ----------
*
*     Author:     Benno List, 1.9.94
*     -------
*
*     Changed by:
*     -----------
*
*     Calling:  
*               draprn         Random number generator
*
*               rgrambo        RAndom Momenta BOoster (phase space gen.)
*               DECPI0       pi0 decay
*               DECRHO       rho decay
*               DECDEL       Delta decay
*               DECETA       eta decay
*
************************************************************************
c #include "rglujets.inc"
c #include "rgludat2.inc"
c #include "rgpdiss.inc"
c #include "rgherup.inc"

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
cc       INTEGER PYK,PYCHGE,PYCOMP
C
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)

       COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
       COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
       COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
       COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)

****************************************************************
c      Double Precision xm1,xw1

c      INTEGER INSTAR, IDEC, NCODES,NDECAY,I,J,NNSTAR,ICHECK
c      INTEGER NCHARG,NDEC1,NDEC2,NSUM,ITER,IPART,IOK,IHELI,IDUMMY
c      REAL R,RANBW
c	Double Precision draprn
c      EXTERNAL draprn

cc      DATA ICHECK/0/
      PARAMETER (NCODES = 9)
      PARAMETER (NDECAY = 6)
      INTEGER M
      DOUBLE PRECISION beta(3)

* Legal N*+ codes
      INTEGER NPLCOD (1:NCODES)
     +        /12212, 2124, 22212, 32212, 2216, 12216, 22124,
     +         42212, 32124/
* Legal N*0 codes
      INTEGER N0COD (1:NCODES)
     +        /12112, 1214, 22112, 32112, 2116, 12116, 21214,
     +         42112, 31214/
* Cummulative branching ratios for the channels
* N* -> N pi, Delta pi, N rho, N (pi pi)_I=0, N eta, Lambda K
* for the different mass states
      REAL DECCHN (1:NDECAY-1, 1:NCODES)
     +        /0.67, 0.92, 0.92, 1.00, 1.00,
     +         0.60, 0.80, 1.00, 1.00, 1.00,
     +         0.50, 0.50, 0.50, 0.50, 1.00,
     +         0.75, 0.80, 0.90, 0.90, 0.90,
     +         0.45, 1.00, 1.00, 1.00, 1.00,
     +         0.67, 0.78, 0.89, 1.00, 1.00,
     +         0.10, 0.45, 0.65, 1.00, 1.00,
     +         0.15, 0.40, 0.55, 0.85, 0.85,
     +         0.15, 0.15, 0.90, 0.90, 0.90/
* Cummulative branching ratios for different particle charges
      REAL CGCOEF (1:2, 1:NDECAY)
     +        /0.33, 1.00,
     +         0.50, 0.83,
     +         0.33, 1.00,
     +         0.67, 1.00,
     +         1.00, 1.00,
     +         1.00, 1.00/
* Decay particles for N*+ decays
      INTEGER  NPLDC1 (1:3, 1:NDECAY), NPLDC2 (1:3, 1:NDECAY),
     +         NPLDC3 (1:3, 1:NDECAY)
      DATA ((NPLDC1(I,J), NPLDC2(I,J), NPLDC3(I,J), I=1, 3), J=1,NDECAY)
* N* -> N pi
     +        /2212,  111,    0,
     +         2112,  211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         2224, -211,    0,
     +         2214,  111,    0,
     +         2114,  211,    0,
* N* -> N rho
     +         2212,  113,    0,
     +         2112,  213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2212,  211, -211,
     +         2212,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2212,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  321,    0,
     +            0,    0,    0,
     +            0,    0,    0/

* Decay particles for N*0 decays
      INTEGER N0DC1 (1:3, 1:NDECAY), N0DC2 (1:3, 1:NDECAY),
     +        N0DC3 (1:3, 1:NDECAY)
      DATA ((N0DC1(I,J), N0DC2(I,J), N0DC3(I,J), I=1, 3), J=1, NDECAY)
* N* -> N pi
     +        /2112,  111,    0,
     +         2212, -211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         1114,  211,    0,
     +         2114,  111,    0,
     +         2214, -211,    0,
* N* -> N rho
     +         2112,  113,    0,
     +         2212, -213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2112,  211, -211,
     +         2112,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2112,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  311,    0,
     +            0,    0,    0,
     +            0,    0,    0/



      REAL*8 XM (100), PMOM (4, 100), WT, SUM, WIDTH, XMMIN
      REAL*8 PDPCMS (1:4), PDPRS (1:4)

      IDEC  = 0

* if particle already decayed --> error 

      IF (K(INSTAR,1) .NE. 1) THEN
        PRINT *, 'DECNST error: ISTAT (', INSTAR, ') = ', K (INSTAR,1)
        RETURN
      END IF

*
* Look up particle code
*
      NNSTAR = 0
      NCHARG = 0
      DO 1 I = 1, NCODES

*** check if part type is N*+ or N*0

        IF (ABS (K(INSTAR,2)) .EQ. NPLCOD (I)) THEN
          NNSTAR = I
          NCHARG = 1
        ELSE IF (ABS (K(INSTAR,2)) .EQ. N0COD (I)) THEN
          NNSTAR = I
          NCHARG = 0
        END IF
 1    CONTINUE

c      PRINT *, 'NNSTAR = ', NNSTAR, ', NCHARG = ', NCHARG
      IF (NNSTAR .EQ. 0) THEN
        PRINT *, '*** DECNST error: Don''t know particle ',
     +           K(INSTAR,2)
c        CALL ERRLOG (201, 'S: DECNST: Wrong particle!')
        RETURN
      END IF

*
* Choose decay channel
*
 7    CONTINUE

      R = PYR(0)
      NDEC1 = 1
 2    CONTINUE
        IF (NDEC1 .EQ. NDECAY .OR. DECCHN (NDEC1, NNSTAR) .GE. R) GOTO 3
        NDEC1 = NDEC1 + 1
      GO TO 2
 3    CONTINUE

c      IF (NDEC1 .LT. NDECAY) THEN
c        PRINT *, 'R = ', R, ', DECCHN (', NDEC1, ', ', NNSTAR, ') = ',
c     +           DECCHN (NDEC1, NNSTAR)
c      ELSE
c        PRINT *, 'R = ', R, ', DECCHN (', NDEC1-1, ', ', NNSTAR, ') = ',
c     +           DECCHN (NDEC1-1, NNSTAR)
c      END IF

      R = PYR(0)
      NDEC2 = 1
 4    CONTINUE
        IF (NDEC2 .EQ. 3 .OR. CGCOEF (NDEC2, NDEC1) .GE. R) GO TO 5
        NDEC2 = NDEC2 + 1
      GO TO 4
 5    CONTINUE

*
* Fill the mass array XM
*

      IDEC  = N + 1
c	write(6,*) ' decnst IDEC',IDEC,NDEC2,N,NSUM

      IF (NCHARG .EQ. 1) THEN
        K(IDEC,2)   = NPLDC1 (NDEC2, NDEC1)
        K(IDEC+1,2) = NPLDC2 (NDEC2, NDEC1)
        K(IDEC+2,2) = NPLDC3 (NDEC2, NDEC1)
      ELSE
        K(IDEC,2)   = N0DC1 (NDEC2, NDEC1)
        K(IDEC+1,2) = N0DC2 (NDEC2, NDEC1)
        K(IDEC+2,2) = N0DC3 (NDEC2, NDEC1)
      END IF
      IF (K(INSTAR,2) .LT. 0) THEN
        K(IDEC,2) = -K(IDEC,2)
        IF (K(IDEC+1,2) .NE. 111 .AND. K(IDEC+1,2) .NE. 221) THEN
          K(IDEC+2,2) = -K(IDEC+2,2)
          K(IDEC+3,2) = -K(IDEC+3,2)
        END IF
      END IF


c Decay to K0: 50% K0S, 50% K0L
      IF (ABS (K(IDEC+1,2)) .EQ. 311) THEN
        IF (PYR(0) .LT. 0.50) THEN
          K(IDEC+1,2) = 130
        ELSE
          K(IDEC+1,2) = 310
        END IF
      END IF

      IF (K(IDEC+2,2) .EQ. 0) THEN
        NSUM = 2
c       PRINT *, ' Decay to ', K(IDEC,2), K(IDEC+1,2)
      ELSE
        NSUM = 3
c       PRINT *, ' Decay to ', K(IDEC,2), K(IDEC+1,2), K(IDEC+2,2)
      END IF

      ITER = 0
c	write(6,*) ' decnst: nsum,instar ',nsum,instar,idec
c	call pylist(1)
 9    CONTINUE
      ITER = ITER + 1
      SUM = 0D0
      DO 10  I = 1, NSUM
cccc        IPART = K(instar+I,2)
        IPART = K(idec-1+I,2)
	  call pdismass(IPART,XM1,XW1)
        WIDTH = XW1
        XM (I) = XM1
c	  write(6,*) ' decnst ', width,xm(i)
c        write(*,*) XM(i)
        XMMIN = XM (I) - 2D0*WIDTH
        IF (ABS (IPART) .EQ. 1114 .OR. ABS (IPART) .EQ. 2114 .OR.
     +      ABS (IPART) .EQ. 2214 .OR. ABS (IPART) .EQ. 2224) THEN
          XMMIN = 1.1D0
        END IF
        IF (WIDTH .GT. 1D-6) THEN
          XM (I) = RANBW (XM (I), WIDTH, XMMIN, XM (I) + 2D0*WIDTH)
        END IF
c	  write(6,*) ' sum,xm ',sum,xm(i),width
        SUM = SUM + XM (I)
 10   CONTINUE


      IF (SUM .GT. P(INSTAR,5) .AND. ITER .LT. 100) GO TO 9
      IF (SUM .GT. P(INSTAR,5)) GO TO 7

*
* Call rgrambo: distribute particles uniformely in phase space
*
         
       do 15 i=1,3
           beta(i)=P(instar,i)/P(instar,4)       
 15     continue

       CALL PYROBO(instar,instar,0.d0,0.d0,-beta(1),-beta(2) 
     &                                ,-beta(3))


       CALL rgrambo (NSUM, P(INSTAR,5), XM, PMOM, WT, 1, IOK)
       IF (IOK .NE. 0) THEN
	    write(6,*) ' decnst: RGRAMBO fails '
c         CALL ERRLOG (202, 'S: DECNST: rgrambo failed!')
       END IF

* Update /PYJETS/
*
 
       do 20 i=1,nsum            
          P(N+i,5)=0
          do 40 j=1,4      
             P(N+i,j)=PMOM(j,i)
c             P(N+i,5)=P(N+i,5)+PMOM(j,i)
             P(N+i,5)=XM(i)
 40       continue
c          P(N+i,4)=PMOM(4,i)
c          P(N+i,5)=P(N+i,4)-P(N+i,5)

 20    continue   
 
c       write(*,*) ' decnst: p(instar) ', P(instar,1),P(instar,2),
c     &              P(instar,3),P(instar,4),P(instar,5)
c       do i=1,nsum
c       write(*,*) ' decnst: p(n+i) ',i,P(N+i,1),P(N+i,2),P(N+i,3),
c     & P(N+i,4),P(N+i,5) 
c       enddo

       CALL PYROBO(instar,instar,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))
       

       M=N

      K(INSTAR,1) = 11
      K(INSTAR,4) = N + 1
      K(INSTAR,5) = N + NSUM
      M=N
      N=N+NSUM
c	write(6,*) ' decnst ',instar,K(INSTAR,4),K(INSTAR,5),N,NSUM
      DO 50, I = 1, NSUM
        P (N + I,5) = XM (I)
        K (M + I,1) = 1
        K (M + I,3) = INSTAR
        K (M + I,4) = 0
c        K (M + I,5) = 0
 50   CONTINUE

*       N = N + NSUM

       CALL PYROBO(M+1,N,0.d0,0.d0,beta(1),beta(2),beta(3))
      
c       write(*,*) 'rgrambo dcnst'
c
c       call pylist(1)


       ijetset=0     ! jetset performs resonanse decays

*
* Perform Delta, pi0, rho, eta decay
*

      IF (ABS (K(IDEC,2)) .EQ. 1114 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2114 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2214 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2224) THEN
 	   K(IDEC,1)=1
c      if(ijetset.eq.1) CALL DECDEL (IDEC)
      END IF

      IF (K(IDEC+1,2) .EQ. 111) THEN
	   K(IDEC+1,1)=1	
c         if(ijetset.eq.1) CALL DECPI0 (IDEC+1)
      ELSE IF (K(IDEC+1,2) .EQ. 113 .OR.
     +         ABS (K(IDEC+1,2)) .EQ. 213) THEN
        IHELI = INT (3.0*PYR(0)) - 1
	   K(IDEC+1,1)=1
c        if(ijetset.eq.1) CALL DECRHO (IDEC+1, IHELI, IDUMMY)
      ELSE IF (K(IDEC+1,2) .EQ. 221) THEN
	   K(IDEC+1,1)=1
c        if(ijetset.eq.1) CALL DECETA (IDEC+1)
      END IF

      IF (K(IDEC+2,2) .EQ. 111) THEN
	   K(IDEC+2,1)=1
c        if(ijetset.eq.1) CALL DECPI0 (IDEC+2)
      END IF
	
cc      if(ijetset.eq.0.and.ILHA.eq.0) call pyexec
      if(ijetset.eq.0) call pyexec
      
      
      K(INSTAR,1) = 11
c      call pylist(2)
      RETURN
      END
c ------------------------------------------------------------------------------
c
c ------------------------------------------------------------------------------
      SUBROUTINE RGRAMBO(N,ET,XM,P,WT,LW,IOK)
c double precision version of RAMBO	
C.      A NEW MONTE-CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT HIGH
C.      ENERGIES - SUBROUTINE RAMBO (RANDOM MOMENTUM BOOSTER).
C.      AUTHORS:  R.KLEISS AND W.J.STIRLING, CERN, GENEVA
C.                S.D.ELLIS, U.WASHINGTON, SEATTLE.
C.
C.                 DOCUMENTATION = CERN-TH.4299/85, OCTOBER 1985
C. published as: KLEISS, R., W.J. STIRLING & S.D. ELLIS (1986):
C.   A new Monte Carlo treatment of multiparticle phase space at high
C.   energies. - Comp. Phys. Comm. 40, 359 - 373, Amsterdam.
C.             ORIGIN= R.KLEISS '86, STATUS= USED, USES IMPLICIT REAL*8.
C------------------------------------------------------
C
C                       RAMBO
C
C             RA(NDOM)  M(OMENTA)  BO(OSTER)
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=N )
C    P  = PARTICLE MOMENTA ( DIM=(4,N) )
C    WT = WEIGHT OF THE EVENT
C    LW = FLAG FOR EVENT WEIGHTING:
C         LW = 0 WEIGHTED EVENTS
C         LW = 1 UNWEIGHTED EVENTS ( FLAT PHASE SPACE )
C
C CHANGES MADE: RN --> draprn 
C BL: Iteration counting (LOOPS); return with a warning if LOOPS > LPMAX
C BL: IOK=0 indicates that all was O.K.; use H1STOP
C BL: Use ERRLOG
C------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
c      Double Precision draprn
c      EXTERNAL draprn
	
      DIMENSION XM(N),P(4,N),Q(4,100),Z(100),R(4),
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(6)
      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/,IWARN/6*0/
      PARAMETER (LPMAX = 10000)
      PARAMETER (TWOPI = 6.283 185 307 179 586 476 D0)
      DATA  PO2LOG /4.515 827 052 93 D-1/
      CHARACTER ERTEXT*132
C
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      PO2LOG = DLOG (TWOPI / 4D0)
      Z(2)=PO2LOG
      DO 101 K=3,100
  101 Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))
      DO 102 K=3,100
  102 Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))
C
C CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      WRITE (ERTEXT, 1001) N
c      CALL ERRLOG (270, ERTEXT)
      STOP
C
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+DABS(XM(I))
      IF(XMT.LE.ET) GOTO 106
      PRINT 1002,XMT,ET
      WRITE (ERTEXT, 1002) XMT,ET
c      CALL ERRLOG (271, ERTEXT)
      STOP
C
C CHECK ON THE WEIGHTING OPTION
  106 IF(LW.EQ.1.OR.LW.EQ.0) GOTO 200
      PRINT 1003,LW
      WRITE (ERTEXT, 1003) LW
c      CALL ERRLOG (272, ERTEXT)
      STOP
C
C THE PARAMETER VALUES ARE NOW ACCEPTED
C
  200 LOOPS = 0
      IOK = 0
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
  201 LOOPS = LOOPS + 1
      DO 202 I=1,N
      C=2.*PYR(1)-1.
      S=DSQRT(1.-C*C)
      F=TWOPI*PYR(2)
      Q(4,I)=-DLOG(PYR(3)*PYR(4))
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*DCOS(F)
  202 Q(1,I)=Q(4,I)*S*DSIN(F)
C
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS
C
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
C
C RETURN FOR UNWEIGHTED MASSLESS MOMENTA
      WT=1.D0
      IF(NM.EQ.0.AND.LW.EQ.1) RETURN
C
C CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*DLOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IOK = 1
      WRITE (ERTEXT, 1004) WT
c      CALL ERRLOG (273, ERTEXT)
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IOK = 2
      WRITE (ERTEXT, 1005) WT
c      CALL ERRLOG (274, ERTEXT)
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
C
C RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=DEXP(WT)
      RETURN
C
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=DSQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=DSQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(DABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX,ACCU
      WRITE (ERTEXT, 1006) ITMAX,ACCU
c      CALL ERRLOG (275, ERTEXT)
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
C
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)
      IF(LW.EQ.1) GOTO 401
C
C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IOK = 3
      IF(IWARN(3).LE.5) PRINT 1004,WT
      WRITE (ERTEXT, 1004) WT
c      CALL ERRLOG (273, ERTEXT)
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IOK = 4
      IF(IWARN(4).LE.5) PRINT 1005,WT
      WRITE (ERTEXT, 1005) WT
c      CALL ERRLOG (274, ERTEXT)
      IWARN(4)=IWARN(4)+1
  310 WT=DEXP(WT)
      RETURN
C
C UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT
  401 WT=DEXP(WTM)
      IF(NM.GT.1) GOTO 402
C
C ONE MASSIVE PARTICLE
      WTMAX=XMAX**(4*N-6)
      GOTO 405
  402 IF(NM.GT.2) GOTO 404
C
C TWO MASSIVE PARTICLES
      SM2=0.
      PM2=0.
      DO 403 I=1,N
      IF(XM(I).EQ.0.D0) GOTO 403
      SM2=SM2+XM2(I)
      PM2=PM2*XM2(I)
  403 CONTINUE
      WTMAX=((1.-SM2/(ET**2))**2-4.*PM2/ET**4)**(N-1.5)
      GOTO 405
C
C MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY
  404 WTMAX=XMAX**(2*N-5+NM)
C
C DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT
  405 W=WT/WTMAX
      IF(W.LE.1.001D0) GOTO 406
      IOK = 5
      IF(IWARN(5).LE.5) PRINT 1007,WTMAX,W
      WRITE (ERTEXT, 1007) WTMAX,W
c      CALL ERRLOG (276, ERTEXT)
      IWARN(5)=IWARN(5)+1
  406 CONTINUE
      IF(W.LT.PYR(5) .AND. LOOPS .LE. LPMAX) GOTO 201
      WT=1.D0
      IF (LOOPS .LE. LPMAX) GO TO 407
      IOK = 6
      IF(IWARN(6).LE.5) PRINT 1008,LPMAX,W
      WRITE (ERTEXT, 1008) LPMAX,W
c      CALL ERRLOG (277, ERTEXT)
      IWARN(6)=IWARN(6)+1
  407 CONTINUE
      RETURN
 1001 FORMAT('F: RGRAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT('F: RGRAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1003 FORMAT('F: RGRAMBO FAILS: LW=',I3,' IS NOT AN ALLOWED OPTION')
 1004 FORMAT('W: RGRAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT('W: RGRAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT('W: RGRAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D15.6)
 1007 FORMAT('W: RGRAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT =',D15.6,
     . '     EXCEEDED BY A FACTOR ',D15.6)
 1008 FORMAT('W: RGRAMBO WARNS: MORE THAN',I7,' ITERATIONS. LAST W=',
     . D15.6)
      END
c ----------------------------------------------------------------------
c
c ---------------------------------------------------------------------
      REAL*8 FUNCTION RANBW (ER, GAMMA, EMIN, EMAX)
      
      IMPLICIT NONE
************************************************************************
*
*     Purpose:   Generate random number with BREIT-WIGNER distribution
*     --------
*
*     Input:     ER:     Maximum of distribution
*     ------     GAMMA:  Width of distribution
*                EMIN:   Minimal value of RANBW
*                EMAX:   Maximal value of RANBW
*
*     Output:    Function value: Random number between EMIN and EMAX
*     -------    with BREIT-WIGNER distribution:
*                1 / ((E - ER)**2 + GAMMA**2/4)
*
*     Author:    B. List,  18.7.93
*     -------
*
************************************************************************
      REAL*8 ER, GAMMA, EMIN, EMAX,A,B,E  !,DRAPRN
      double precision PYR

      IF (GAMMA .LT. 1.0E-3*ER) THEN
         RANBW = SNGL(ER)
      ELSE
         A = ATAN (2.0*(EMAX - ER)/GAMMA)
         B = ATAN (2.0*(EMIN - ER)/GAMMA)
         E = ER + 0.5*GAMMA*TAN (PYR(0)*(A - B) + B)
         IF (E .LE. EMAX) THEN
            RANBW = SNGL(E)
         ELSE
            RANBW = SNGL(EMAX)
         END IF
      END IF

      RETURN
      END
c ----------------------------------------------------
c
c ---------------------------------------------------
