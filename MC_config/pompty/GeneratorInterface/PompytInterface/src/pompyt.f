c -----------------------------------------------------
c
c ----------------------------------------------------
       subroutine POMINI(comenergy)
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
cc       INTEGER PYK,PYCHGE,PYCOMP
C
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
       COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
       COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)

       common/pompyt1/ pomi(20),ipom(20)
       common/pomcut1/ pomcut(20)
       external pomfluxb

       pomi(1)=comenergy

       ipom(19)=1    !!!!  Pomeron
cc       ipom(19)=2    !!!!  Reggeon

       if(pomcut(2).gt.0.21) then
         write(6,*) 'TOO high xi_max=',pomcut(2)
         stop
       endif

       if(ipom(1).ne.1.and.ipom(1).ne.2) then
        write(6,*) 'pomPdiss must be 1 or 2'
        stop
       endif

       if(pomi(20).lt.0.0) then
        write(6,*) 'pdiss_ap<0'
        stop
       endif
       if(pomi(20).gt.0.5) then
        write(6,*) 'pdiss_ap>0.5'
        stop
       endif

       if(pomi(19).lt.0.0) then
        write(6,*) 'pdiss_b0<0'
        stop
       endif

       write(6,*) 'POMINI **************************************'
       do i=1,20
       write(6,*) i,' ipom ',ipom(i), ' pomi ',pomi(i),
     X              ' pomcut ',pomcut(i)
       enddo


       MSTP(171)=1
       MSTP(172)=2
       MSTP(173)=0


c pomeron
       p(1,5)=0
       p(1,1)=0
       p(1,2)=0
       p(1,3)=pomcut(2)*pomi(2)*pomi(1)/2.D0
       p(1,4)=dSQRT(P(1,5)**2+P(1,1)**2+P(1,2)**2+P(1,3)**2)
c proton
       p(2,5)=0
       p(2,1)=0
       p(2,2)=0
       p(2,3)=-pomi(2)*pomi(1)/2.D0
       p(2,4)=dSQRT(P(2,5)**2+P(2,1)**2+P(2,2)**2+P(2,3)**2)



       CALL PYINIT('5MOM','pomeron','p',WIN)


cc       flux integration

cc       write(6,*) 'flux' ,pomfluxb(0.01D0),pomcut(1),pomcut(2)
       eps=0.001
       CALL DGADAP(pomcut(1),pomcut(2),pomfluxb,eps,result)
ccc       CALL DGADAP(0.001D0,0.01D0,test,eps,sum)
       pomi(9)=result

ccc       write(6,*) 'flux', pomi(9)




       return
       end
c ------------------------------------------------
c
c -------------------------------------------------
      subroutine POMPDF(x,Q2,xpq)
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       common/pompyt1/ pomi(20),ipom(20)

       CHARACTER*20 PARM(20)
       DOUBLE PRECISION valu(20)

       dimension xpq(-6:6)
       dimension F2(2),FL(2),C2(2),CL(2)
       data ifirst/0/

       if(ifirst.eq.0) then
        ifirst=1
c  pomeron
        if(ipom(19).eq.1) then
         bb=0.5
         qq=10.
         call QCD_2006(bb,qq,2,XPQ,F2,FL,C2,CL)
        else
c reggeon
            parm(1)='DEFAULT'
ccc            parm(1)='NPTYPE'
            valu(1)=211
            CALL PDFSET(PARM,VALU)
        endif
       endif


ccc H1 pdfs pomeron
      if(ipom(19).eq.1) then  ! pomeron
       call QCD_2006(x,Q2,0,XPQ,F2,FL,C2,CL)
      else ! reggeon
cc         CALL STRUCTF(x,dsqrt(q2),UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU)
         CALL STRUCTM(x,dsqrt(q2),UPV,DNV,USEA,DSEA,STR,CHM,BTM,TOP,GLU)
ccc         VAL = UPV/2.D0
         XPQ(0)= GLU
         XPQ(1)= (0.5D0*DNV)+DSEA
         XPQ(-1)= (0.5D0*DNV)+DSEA
         XPQ(2)= (0.5D0*UPV)+USEA
         XPQ(-2)= (0.5D0*UPV)+USEA
         XPQ(3)= STR
         XPQ(-3)= STR
         XPQ(4)= CHM
         XPQ(-4)= CHM
         XPQ(5)= BOT
         XPQ(-5)= BOT
         XPQ(6)= TOP
         XPQ(-6)= TOP
      endif

       return
       end
c ------------------------------------------------------
c
c ------------------------------------------------------
       subroutine POMPYT
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
cc       INTEGER PYK,PYCHGE,PYCOMP
C
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
       COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

       common/pompyt1/ pomi(20),ipom(20)
       common/pomcut1/ pomcut(20)
       double precision MX,MY,pp(5),pb1(5),pb2(5)



 1111  continue
       if(ipom(19).eq.1) then
        call pomflux(xp,t)
       else
        call regflux(xp,t)
       endif

       pomi(11)=xp
       pomi(12)=t

       MY=pymass(2212)
       if(ipom(1).eq.2) then   ! p-diss
         xxxmin=pomcut(11)
         if(xxxmin.lt.1.17455D0) xxxmin=1.17455D0
         MY=xxxmin*(pomcut(12)/xxxmin)**PYR(0)
       endif
c Y momentum
       PHI=6.283185*PYR(0)
       PT=dsqrt((-t)*(1.d0-xp))
       pp(5)=MY
       pp(1)=PT*cos(phi)
       pp(2)=PT*sin(phi)
       pp(3)=pomi(2)*((1.D0-xp)*pomi(1)/2.D0-
     X        (MY**2+PT**2)/(1.D0-xp)/2.D0/pomi(1))
       pp(4)=dsqrt(pp(5)**2+pp(1)**2+pp(2)**2+pp(3)**2)
c beam1
       pb1(5)=pymass(2212)
       pb1(1)=0
       pb1(2)=0
       pb1(3)=pomi(2)*pomi(1)/2.D0
       pb1(4)=dSQRT(pb1(5)**2+pb1(1)**2+pb1(2)**2+pb1(3)**2)
c beam2
       pb2(5)=pymass(2212)
       pb2(1)=0
       pb2(2)=0
       pb2(3)=-pomi(2)*pomi(1)/2.D0
       pb2(4)=dSQRT(pb2(5)**2+pb2(1)**2+pb2(2)**2+pb2(3)**2)
c pomeron mass
       tt=(pb1(4)-pp(4))**2-(pb1(1)-pp(1))**2
     X   -(pb1(2)-pp(2))**2-(pb1(3)-pp(3))**2
ccc       p(4,5)=-dSQRT(dABS(tt))
C pomeron
       p(1,5)=-dSQRT(dABS(tt))
       p(1,1)=-pp(1)
       p(1,2)=-pp(2)
       p(1,3)=pb1(3)-pp(3)
       p(1,4)=pb1(4)-pp(4)
cc       p(1,4)=dSQRT(P(1,5)**2+P(1,1)**2+P(1,2)**2+P(1,3)**2)
C beam proton
       p(2,5)=pymass(2212)
       p(2,1)=0
       p(2,2)=0
       p(2,3)=-pomi(2)*pomi(1)/2.D0
       p(2,4)=dSQRT(P(2,5)**2+P(2,1)**2+P(2,2)**2+P(2,3)**2)

       MX=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     X   -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2
       MX=dsqrt(MX)

       call pyevnt
       if(msti(61).ne.0) go to 1111

cccccccccccccccccc       write(6,*) 'MXMY',MX,MY

       if((MX+MY).gt.(pomi(1)*0.9D0)) then
        write(6,*) 'TOO BIG MASSES', MX,MY,pomi(1)
        stop
       endif


       CALL POMOVE(1,3)
       K(4,3)=1
       K(5,3)=2
C proton, Y
       if(ipom(1).eq.1) then
       p(3,5)=pymass(2212)
       k(3,1)=1
       k(3,2)=2212
       k(3,3)=1
       elseif(ipom(1).eq.2) then    ! pdiss     
       p(3,5)=MY
       k(3,2)=9902210 
       k(3,1)=21 
       k(3,3)=1
      endif
       p(3,1)=pp(1)
       p(3,2)=pp(2)
       p(3,3)=pp(3)
       p(3,4)=pp(4)
c real beams
       p(1,5)=pymass(2212)
       p(1,1)=0
       p(1,2)=0
       p(1,3)=pomi(2)*pomi(1)/2.D0
       p(1,4)=dSQRT(P(1,5)**2+P(1,1)**2+P(1,2)**2+P(1,3)**2)
       k(1,1)=21
       k(1,2)=2212
       k(1,3)=0 

       p(2,5)=pymass(2212)
       p(2,1)=0
       p(2,2)=0
       p(2,3)=-pomi(2)*pomi(1)/2.D0
       p(2,4)=dSQRT(P(2,5)**2+P(2,1)**2+P(2,2)**2+P(2,3)**2)
       k(2,1)=21
       k(2,2)=2212
       k(2,3)=0 
ccccccccccccccccccccccccccccccccccccccccccccccccc

c proton dissociation
       if(ipom(1).eq.2) then
        call prodiffvm
        call PDISDCVM
cc        call pyexec
       endif



       pomi(13)=p(4,5)   ! pomeron mass
       pomi(14)=p(3,5)   ! MY
       pomi(15)=MX   ! MX


       return
       end
c ---------------------------
c
c -------------------------------
       subroutine pomflux(xp,t)
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       common/pomcut1/ pomcut(20)
       common/pompyt1/ pomi(20),ipom(20)

c fit B pomeron
       data a0 /1.1110/, b0 /5.5/, ap /0.06/
       data ifirst/0/

       if(ifirst.eq.0) then
         ifirst=1
         if(ipom(1).eq.2) then
          a0=1.1110
          b0=pomi(19)
          ap=pomi(20)
         endif
       endif

ccc       write(6,*) '*********** ',a0,b0,ap

       AAA=2.D0*(1.D0-a0) 

 3     continue
       xp=(PYR(0)*(pomcut(2)**AAA-pomcut(1)**AAA)
     X    +pomcut(1)**AAA)**(1.D0/AAA)
       t= LOG(PYR(0)*( exp(b0*pomcut(4))-exp(b0*pomcut(3)) ) 
     X    +exp(b0*pomcut(3)) )/b0
       
       f=exp(-2.D0*ap*LOG(xp)*t)
       if(PYR(0).gt.f) go to 3

       return
       end
c ---------------------------
c
c -------------------------------
       subroutine regflux(xp,t)
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       common/pomcut1/ pomcut(20)
       common/pompyt1/ pomi(20),ipom(20)

c fit B reggeon
       data a0 /0.5/, b0 /1.6/, ap /0.3/
       data ifirst/0/

       if(ifirst.eq.0) then
         ifirst=1
         if(ipom(1).eq.2) then
          a0=1.1110
          b0=pomi(19)
          ap=pomi(20)
         endif
       endif

ccc       write(6,*) '*********** ',a0,b0,ap

       AAA=2.D0*(1.D0-a0) 

 3     continue
       xp=(PYR(0)*(pomcut(2)**AAA-pomcut(1)**AAA)
     X    +pomcut(1)**AAA)**(1.D0/AAA)
       t= LOG(PYR(0)*( exp(b0*pomcut(4))-exp(b0*pomcut(3)) ) 
     X    +exp(b0*pomcut(3)) )/b0
       
       f=exp(-2.D0*ap*LOG(xp)*t)
       if(PYR(0).gt.f) go to 3

       return
       end
c -----------------------------------------
c
c ------------------------------------------
       subroutine POMOUT
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
cc       INTEGER PYK,PYCHGE,PYCOMP
       common/pompyt1/ pomi(20),ipom(20)
        call pystat(1)
       if(ipom(1).eq.1) then
        write(6,*) 'Integrated flux ', pomi(9)
       endif
       if(ipom(1).eq.2) then
        write(6,*) 'No reason to believe in these cross sections'
       endif

       call PDFsta
       return
       end
c ----------------------------------------
c
c -----------------------------------------
      SUBROUTINE POMOVE(NSTART,NSHIFT)

C...Shift line numbers NSTART to N in event record NSHIFT lines forward ,
C...keeping line pointer information intact.
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       IMPLICIT INTEGER(I-N)
cc       INTEGER PYK,PYCHGE,PYCOMP

       COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
       COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)



      DO 100 I=N,NSTART,-1
      DO 10 J=1,5
      K(I+NSHIFT,J)=K(I,J)
      P(I+NSHIFT,J)=P(I,J)
   10 V(I+NSHIFT,J)=V(I,J)
C...Correct line pointers in K(i,3), K(i,4), K(i,5).
      IF(K(I,3).GT.0) K(I+NSHIFT,3)=K(I,3)+NSHIFT
      DO 20 J=4,5
      IF(K(I,J).GE.NSTART.AND.K(I,J).LT.MSTU(5)) THEN
        K(I+NSHIFT,J)=K(I,J)+NSHIFT
      ELSEIF(K(I,J).GE.MSTU(5)) THEN
        ICFR=MOD(K(I,J),MSTU(5))
        IF(ICFR.GE.NSTART) ICFR=ICFR+NSHIFT
        ICTO=MOD(K(I,J)/MSTU(5),MSTU(5)**2)
        IF(ICTO.GE.NSTART) ICTO=ICTO+NSHIFT
        MCTO=MOD(K(I,J)/MSTU(5)**2,2*MSTU(5)**2)
        IF(MCTO.GE.NSTART) MCTO=MCTO+NSHIFT
        MCFR=    K(I,J)/(2*MSTU(5)**2)
        IF(MCFR.GE.NSTART) MCFR=MCFR+NSHIFT
        K(I+NSHIFT,J)=2*MSTU(5)**2*MCFR+MSTU(5)**2*MCTO
     &  +MSTU(5)*ICTO+ICFR
      ENDIF
   20 CONTINUE
      DO 30 J=1,5
      K(I,J)=0
      P(I,J)=0.
   30 V(I,J)=0.
  100 CONTINUE
      N=N+NSHIFT

      IF(MSTU(70).GT.0) THEN
C...Shift pointers for `======' lines.
        DO 200 I=71,MAX(80,MSTU(70))
  200   IF(MSTU(I).GE.NSTART) MSTU(I)=MSTU(I)+NSHIFT
      ENDIF
      IF(MSTU(90).GT.0) THEN
C...Shift pointers for heavy flavours.
        DO 210 I=91,98
  210   IF(MSTU(I).GE.NSTART) MSTU(I)=MSTU(I)+NSHIFT
      ENDIF

      RETURN
      END
c -----------------------------------------------------------------
c
c ----------------------------------------------------------------

C#######################################################################
C
C   One- and two-dimensional adaptive Gaussian integration routines.
C
C **********************************************************************

      SUBROUTINE DGADAP(A0,B0,F,EPS,SUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X)
C   METHOD            - ADAPTIVE GAUSSIAN
C   USAGE             - CALL GADAP(A0,B0,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)
C               B0    - UPPER LIMIT (INPUT,REAL)
C               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE
C                       SUPPLIED BY THE USER. (INPUT,REAL FUNCTION)
C               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS
C                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - F
C   AUTHOR            - T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
C   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
C
CAE Changed to double precision
C
      COMMON/DGADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
    1 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)

ccc      write(6,*) 'III',A0,B0,eps

      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.3
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(0.5*(A0+B0))
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
  100 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(B(I)-A(I+2)+W1)
      F1(I+2)=F(U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(B(I+2)+A(I+2)-U2)
      F1(I+3)=F(A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 130
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 120
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
  110 CONTINUE
      IF(L.EQ.1) GOTO 130
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 100
      I=I-1
      L=L-1
      GOTO 110
  120 WRITE(6,1)
  130 RETURN
      END
c -------------------------------------------------------
c
c --------------------------------------------------------
      function pomfluxb(xp)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      common/pomcut1/ pomcut(20)
      common/pompyt1/ pomi(20),ipom(20)

      if(ipom(19).eq.1) then  !!! pomeron dlux
       tmax=-3.0
       call h12006flux(xp,tmax,1,2,1,flux)
      else   !!! reggeon flux
       tmax=-3.0
       call h12006flux(xp,tmax,1,2,2,flux)
      endif

cccc      write(6,*) xp,flux

      pomfluxb=flux
      return

c fit B pomeron
c       data a0 /1.1110/, b0 /5.5/, ap /0.06/
c       data ifirst/0/
c
c       if(ifirst.eq.0) then
c         ifirst=1
c         if(ipom(1).eq.2) then
c          a0=1.1110
c          b0=pomi(19)
c          ap=pomi(20)
c         endif
c       endif
c 
c      dp=2.*a0-1.
c      B=b0+2*ap*log(1./xp)
c      fff=(1./xp)**dp *(1/B)*(EXP(B*pomcut(4))-EXP(B*pomcut(3)))
c      B003=b0+2*ap*log(1./0.003)
c      fff003=(1./0.003)**dp *(1/B003)*
c     X        (EXP(B*pomcut(4))-EXP(B*pomcut(3)))
c
c      pomfluxb=fff/(0.003*fff003)
c      return
      end


