!DIR$ FREEFORM
!CC*******************************************************************************
!C  UMAT: developed by S.C. Soare. 
!C  IT USES ISOTROPIC ELASTICITY COUPLED WITH A Neural Network model of ANISOTROPIC YIELD FUNCTION. 
!C  THE SUBROUTINE IS DESIGNED FOR PLANE (2D) STRESS STATES. 
!C  IT IS BASED ON THE FULLY IMPLICIT RETURN MAPPING ALGORITHM (with quadratic line search) 
!C  ONE STATE VARIABLE : HARDENING PARAMETER (THE EQUIVALENT PLASTIC STRAIN).
!C  IT IS ASSUMED THAT THE USER HAS DEFINED (USING THE *ORIENTATION OPTION IN ABAQUS)
!C  A LOCAL COORDINATE SYSTEM THAT ROTATES WITH THE MATERIAL (SO THAT DEROT IS ALWAYS
!C  EQUAL WITH THE IDENTITY MATRIX).
!C*****************************************************************************
!C The vector PROPS(NPROPS) contains the material properties as defined in the 
!C  *MATERIAL=USER option in ABAQUS in the following order
!C  PROPS(1) = EMOD  (Elasticity: Young modulus)
!C  PROPS(2) = MU  (Elasticity: Poisson ratio)
!C  Hardening laws defined: 
!C  Swift (power-)law: sigma^bar = a*(b + ep^bar)**c
!C  Voce (exp-)law: sigma^bar = a - b*exp(-c*ep^bar)  (default)
!C  Read further naming/renaming convention in the HARDENING section of this code 
!C  (more specific hardening laws can be implemented in the same section) 
!C  PROPS(3) = a 
!C  PROPS(4) = b
!C  PROPS(5) = c
!C  PROPS(6),...,PROPS(NPROPS): PARAMETERS OF THE YIELD FUNCTION (Neural Network model) 
!C  Yield func implemented: nnYS
!C  The parameters of the nnYS model are as follow:
!C  PROPS(6) = number of  layers (must be >=2)
!C  PROPS(7) = number of computational units/neurons in the first (hidden) layer 
!C  PROPS(8) = homogeneity degree of the activation func of the first layer 
!C  ......
!C  PROPS(7+2*PROPS(7)) =  number of computational units/neurons in the last (hidden) layer
!C  PROPS(7+2*PROPS(7)+1),..., PROPS(NPROPS) = the weights of the NN-model 
!C!************************************************************************************** 

!CCCC-----NOTE: the UMAT interface may vary with the ABAQUS version 	
!C	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!C      RPL,DDSDDT,DRPLDE,DRPLDT, &
!C      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
!C      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROTT,PNEWDT, &
!C      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

!CCC---NOTE: the INCLUDE directive enforces implicit casting in conflict with 'IMPLICIT NONE'
!CC          use 'IMPLICIT NONE' in the testing/implementation phase and then comment it out 
!    IMPLICIT NONE      	
	INCLUDE 'ABA_PARAM.INC'
!C    INTEGER, PARAMETER :: PREC =  SELECTED_REAL_KIND(15,307)
    INTEGER, PARAMETER :: PREC = 8
!C******************************************************************************
!C  VARIABLES REQUIRED BY ABAQUS (THE ARGUMENTS OF THE SUBROUTINE)
!C  FOR A DESCRIPTION OF THE LIST OF VARIABLES SEE ABAQUS MANUAL (VOL. VI)

!C    !!!CHARACTER(80)::  CMNAME
	CHARACTER*80 CMNAME
    REAL(PREC)::SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT
    INTEGER::NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC
    REAL(PREC),DIMENSION(NTENS):: STRESS,DDSDDT,DRPLDE,STRAN,DSTRAN
    REAL(PREC),DIMENSION (NTENS, NTENS) :: DDSDDE 
    REAL(PREC),DIMENSION(NSTATV) :: STATEV
    REAL(PREC),DIMENSION(NPROPS) :: PROPS
    REAL(PREC),DIMENSION(3,3) :: DFGRD0, DFGRD1, DROT
    REAL(PREC),DIMENSION(3) :: COORDS
    REAL(PREC),DIMENSION(2) :: TIME
    REAL(PREC),DIMENSION(1) :: PREDEF, DPRED
    INTEGER,DIMENSION(4)::JSTEP
		
!C    CHARACTER*80 CMNAME
!C	INTEGER::NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC
!C	REAL(PREC)::SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT
!C    REAL(PREC),DIMENSION::STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS), &
!C	DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS), &
!C	TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS), &
!C	COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!C	INTEGER,DIMENSION::JSTEP(4)	
	
!C*******************************************************************************
!C  INTERNAL VARIABLES OF THE SUBROUTINE
    REAL(PREC), PARAMETER::ZERO=0.0D0

!C  elastic constants
    REAL(PREC) :: EMOD, ENU
!C  COMPLIANCE TENSOR
    REAL(PREC),DIMENSION(NTENS,NTENS)::SCOMP

!C  HARDENING PARAMETERS
    REAL(PREC)::AA,BB,CC
!C  HARDENING VALUES
    REAL(PREC):: HF, HPF

!C  STRESS TENSOR AND ITS INCREMENTS
    REAL(PREC),DIMENSION(NTENS)::SIGMA, DSIGMA, D2SIGMA

!C  EQUIVALENT PLASTIC STRAIN AND ITS INCREMENTS
    REAL(PREC):: EPBAR, DEPBAR, D2EPBAR

!C  YIELD FUNCTION VALUE, GRADIENT AND HESSIAN
    REAL(PREC):: YF
	REAL(PREC),DIMENSION(NTENS)::GYF
	REAL(PREC),DIMENSION(NTENS,NTENS)::HYF

!C  CONVERGENCE TOLERANCES
!    REAL(PREC),PARAMETER::TOL1=1.0E-006, TOL2=1.0E-008
	REAL(PREC),PARAMETER::TOL1=1.0E-007

!C  TEMPORARY HOLDERS
    REAL(PREC)::TT, TTA, TTB, ZALPHA, F1, FZERO, TDEPBAR
	REAL(PREC),DIMENSION(NTENS)::YVECTOR, F2, TDSIGMA
	REAL(PREC),DIMENSION(NTENS,NTENS)::XIMAT
	REAL(PREC),DIMENSION(6)::BV
	REAL(PREC),DIMENSION(3)::ZZ

!C  LOOP COUNTERS
    INTEGER::K1,K2,NRK,JJ,KK

!C  NEWTON-RAPHSON MAXIMUM NUMBER OF ITERATIONS
    INTEGER,PARAMETER:: NRMAX=100  
	
!C NN-model interface variables 
!C    integer, parameter::nnZero=6
!C	integer,parameter::inputDim=NTENS
    INTEGER::nnLayers,nUnits, nAllUnits, firstWeight	
	integer,allocatable::vUnits(:)
    real(prec),allocatable::vDegree(:)
	real(prec):: nnHom, homDegReciproc  

!C*****************************************************
	EMOD = PROPS(1)
    ENU = PROPS(2)
    AA = PROPS(3)
    BB = PROPS(4)
    CC = PROPS(5)
	
	nnLayers=int(PROPS(6))
	allocate(vUnits(nnLayers),vDegree(nnLayers))
	vUnits(1)=NTENS !!! input dimension
	vDegree(1)=0.0d0 !! not used
	nAllUnits=NTENS
	nnHom=dble(1.0)
	JJ=7
	do KK=2,nnLayers
	  nUnits=int(PROPS(JJ))
	  vUnits(KK)=nUnits
	  nAllUnits=nAllUnits+nUnits
	  vDegree(KK)=PROPS(JJ+1)
	  nnHom=nnHom*vDegree(KK)
	  JJ=JJ+2
	end do  
	firstWeight=JJ
	homDegReciproc=1.0d0/nnHom
    
			
!C!********************************************
!C RECOVER THE EQUIVALENT PLASTIC STRAIN AT THE BEGINING OF THE INCREMENT
    EPBAR = STATEV(1)
      
!C!********************************************
!C INITIALIZE THE STIFFNESS TENSOR (IT WILL BE STORED IN DDSDDE)
    DDSDDE = ZERO

!C ELASTIC PROPERTIES (plane stress)
	TT = 0.5D0*EMOD/(1.0D0+ENU)
	TTA=EMOD/(1.0D0-ENU*ENU)
	TTB=ENU*TTA

    DDSDDE(1,1) = TTA
    DDSDDE(2,1) = TTB
    DDSDDE(1,2) = TTB
    DDSDDE(2,2) = TTA
	DDSDDE(3,3) = TT

!C COMPUTE THE TRIAL STRESS : SIGMA_{N+1} = SIGMA_{N} + C[DELTA_EPSILON]
    DO K1=1,NTENS,1
        TT=DDSDDE(K1,1)*DSTRAN(1)+DDSDDE(K1,2)*DSTRAN(2)+DDSDDE(K1,3)*DSTRAN(3)
	    DSIGMA(K1) = TT
	    SIGMA(K1)=STRESS(K1)+TT
    END DO

!C CHECK YIELDING CONDITION
    CALL KHARD(HF,HPF,EPBAR,AA,BB,CC)
!!(vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDegReciproc,nAllUnits,firstWeight,Fval)	
    CALL YFUNCTION(SIGMA,NTENS,PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc,&
	               nAllUnits,firstWeight,YF)

!C  ELASTIC STEP :  UPDATE STRESS
	IF (YF <= HF) THEN
	    STRESS = SIGMA
!C  DDSDDE HAS BEEN DEFINED ABOVE
!C  THE EQUIVALENT PLASTIC STRAIN, STATEV(1), REMAINS UNCHANGED		
        RETURN
	END IF

!C***********************************************
!C MAIN LOOP : RETURN MAPPING ALGORITHM

!C  DEFINE COMPLIANCE (note that it outputs ENGINEERING shears)
	TTA=1.0D0/EMOD
	TTB=-ENU/EMOD 
	SCOMP(1,1)=TTA
	SCOMP(2,1)=TTB
	SCOMP(3,1)=ZERO
	SCOMP(1,2)=TTB
	SCOMP(2,2)=TTA
	SCOMP(3,2)=ZERO
	SCOMP(1,3)=ZERO
	SCOMP(2,3)=ZERO
	SCOMP(3,3)=2.0D0*(1.0D0+ENU)/EMOD

!C    FIRST Newton-Raphson step (no Hessian required)
!C**************************************************************      
!C  DEPBAR=ZERO
!!(vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDegReciproc,nAllUnits,firstWeight,Fval,FGval)
	CALL GYFUNCTION(SIGMA,NTENS,PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc,&
	                nAllUnits,firstWeight,YF,GYF)
	F1=YF-HF

!C  ASSEMBLE XIMAT MATRIX AND Y-VECTOR
    DO K1=1,NTENS,1
	    YVECTOR(K1)=-F1*GYF(K1)
	    DO K2=K1,NTENS,1
	        TT=HPF*SCOMP(K1,K2)+GYF(K1)*GYF(K2)
	        XIMAT(K1,K2)=TT
	        XIMAT(K2,K1)=TT
	    END DO
	END DO

!C  SOLVE FOR STRESS NR-INCREMENT
    BV(1)=DSQRT(XIMAT(1,1))
	BV(2)=XIMAT(2,1)/BV(1)
	BV(3)=XIMAT(3,1)/BV(1)
	BV(4)=DSQRT(XIMAT(2,2)-BV(2)*BV(2))
	BV(5)=(XIMAT(3,2)-BV(2)*BV(3))/BV(4)
	BV(6)=DSQRT(XIMAT(3,3)-(BV(3)*BV(3)+BV(5)*BV(5)))
	ZZ(1)=YVECTOR(1)/BV(1)
	ZZ(2)=(YVECTOR(2)-BV(2)*ZZ(1))/BV(4)
	ZZ(3)=(YVECTOR(3)-(BV(3)*ZZ(1)+BV(5)*ZZ(2)))/BV(6)
	D2SIGMA(3)=ZZ(3)/BV(6)
	D2SIGMA(2)=(ZZ(2)-BV(5)*D2SIGMA(3))/BV(4)
	D2SIGMA(1)=(ZZ(1)-(BV(2)*D2SIGMA(2)+BV(3)*D2SIGMA(3)))/BV(1)

!C  CALCULATE EQUIVALENT PLASTIC STRAIN NR-INCREMENT 
    D2EPBAR=F1
	DO K1=1,NTENS,1
	    D2EPBAR=D2EPBAR+GYF(K1)*D2SIGMA(K1)
	END DO
	D2EPBAR=D2EPBAR/HPF
    
!C  DO LINE SEARCH (along the full NR-step)
    TDEPBAR=D2EPBAR
	TDSIGMA=DSIGMA+D2SIGMA	
    FZERO=F1
    CALL LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO,SCOMP, &
	                   AA,BB,CC,ZALPHA, &
	                   PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc, &
	                   nAllUnits,firstWeight) 
!C  UPDATE
    DEPBAR=ZALPHA*D2EPBAR
	DSIGMA=DSIGMA+ZALPHA*D2SIGMA
	
!C    THE REST OF N-R ITERATIONS
!C******************************************************	     
    DO NRK=1,NRMAX,1

!C      CALCULATE NEW VALUES ASSOCIATED WITH NEW STATE
        CALL KHARD(HF,HPF,EPBAR+DEPBAR,AA,BB,CC)
	    SIGMA=STRESS+DSIGMA
!! (vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDeg,homDegReciproc,nAllUnits,firstWeight,Fval,FGval,FHval)		
	    CALL HYFUNCTION(SIGMA,NTENS,PROPS,NPROPS,vUnits,vDegree,nnLayers,nnHom,homDegReciproc,&
	                nAllUnits,firstWeight,YF,GYF,HYF)
		
	    F1=YF-HF
	    FZERO=F1*F1
	    DO K1=1,NTENS,1
	        TT=DEPBAR*GYF(K1)-DSTRAN(K1)
	        DO K2=1,NTENS,1
	            TT=TT+SCOMP(K1,K2)*DSIGMA(K2)
	        END DO
	        F2(K1)=TT
	        FZERO=FZERO+TT*TT
	    END DO
        FZERO=DSQRT(FZERO)
!C      CHECK TOLERANCES
!C        IF ((DABS(F1)<TOL1).AND.(DSQRT(TTB)<TOL2)) EXIT
        IF(FZERO<TOL1) EXIT

!C      ASSEMBLE XIMAT MATRIX AND Y-VECTOR
        DO K1=1,NTENS,1
	        YVECTOR(K1)=-(F1*GYF(K1)+HPF*F2(K1))
	        DO K2=K1,NTENS,1
	            TT=HPF*(SCOMP(K1,K2)+DEPBAR*HYF(K1,K2))+GYF(K1)*GYF(K2)
	            XIMAT(K1,K2)=TT
                XIMAT(K2,K1)=TT
            END DO
	    END DO

!C      SOLVE FOR STRESS NR-INCREMENT
        BV(1)=DSQRT(XIMAT(1,1))
		BV(2)=XIMAT(2,1)/BV(1)
		BV(3)=XIMAT(3,1)/BV(1)
		BV(4)=DSQRT(XIMAT(2,2)-BV(2)*BV(2))
		BV(5)=(XIMAT(3,2)-BV(2)*BV(3))/BV(4)
		BV(6)=DSQRT(XIMAT(3,3)-(BV(3)*BV(3)+BV(5)*BV(5)))
		ZZ(1)=YVECTOR(1)/BV(1)
		ZZ(2)=(YVECTOR(2)-BV(2)*ZZ(1))/BV(4)
		ZZ(3)=(YVECTOR(3)-(BV(3)*ZZ(1)+BV(5)*ZZ(2)))/BV(6)
		D2SIGMA(3)=ZZ(3)/BV(6)
		D2SIGMA(2)=(ZZ(2)-BV(5)*D2SIGMA(3))/BV(4)
		D2SIGMA(1)=(ZZ(1)-(BV(2)*D2SIGMA(2)+BV(3)*D2SIGMA(3)))/BV(1)

!C      CALCULATE EQUIVALENT PLASTIC STRAIN NR-INCREMENT 
        D2EPBAR=F1
	    DO K1=1,NTENS,1
	        D2EPBAR=D2EPBAR+GYF(K1)*D2SIGMA(K1)
	    END DO
	    D2EPBAR=D2EPBAR/HPF

!C      DO LINE SEARCH
        TDEPBAR=DEPBAR+D2EPBAR
	    TDSIGMA=DSIGMA+D2SIGMA
        CALL LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO,SCOMP, &
		             AA,BB,CC,ZALPHA, &
	                 PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc, &
	                 nAllUnits,firstWeight) 

!C      UPDATE
        DEPBAR=DEPBAR+ZALPHA*D2EPBAR
	    DSIGMA=DSIGMA+ZALPHA*D2SIGMA

	END DO !!! END OF NEWTON-RAPHSON ITERATIONS
        
!C  UPDATE STATE VARIABLE
    STATEV(1)=EPBAR+DEPBAR

!C  UPDATE STRESS
    STRESS = STRESS+DSIGMA

!C************************************** COMPUTE TANGENT MODULUS: DDSDDE

!C  COMPUTE XIMAT MATRIX 
    DO K1=1,NTENS,1
	    DO K2=K1,NTENS,1
	        TT=SCOMP(K1,K2)+DEPBAR*HYF(K1,K2)
	        XIMAT(K1,K2)=TT
	        XIMAT(K2,K1)=TT
	    END DO
	END DO

!C  INVERT XIMAT AND STORE XIMAT^(-1) INTO SCOMP (NO LONGER NEEDED)
    BV(1)=DSQRT(XIMAT(1,1))
	BV(2)=XIMAT(2,1)/BV(1)
	BV(3)=XIMAT(3,1)/BV(1)
	BV(4)=DSQRT(XIMAT(2,2)-BV(2)*BV(2))
	BV(5)=(XIMAT(3,2)-BV(2)*BV(3))/BV(4)
	BV(6)=DSQRT(XIMAT(3,3)-(BV(3)*BV(3)+BV(5)*BV(5)))
	
	ZZ(1)=1.0D0/BV(1)
	ZZ(2)=-(BV(2)*ZZ(1))/BV(4)
	ZZ(3)=-(BV(3)*ZZ(1)+BV(5)*ZZ(2))/BV(6)
	SCOMP(3,1)=ZZ(3)/BV(6)
	SCOMP(2,1)=(ZZ(2)-BV(5)*SCOMP(3,1))/BV(4)
	SCOMP(1,1)=(ZZ(1)-(BV(2)*SCOMP(2,1)+BV(3)*SCOMP(3,1)))/BV(1)

!	ZZ(1)=ZERO
	ZZ(2)=1.0D0/BV(4)
	ZZ(3)=-(BV(5)*ZZ(2))/BV(6)
	SCOMP(3,2)=ZZ(3)/BV(6)
	SCOMP(2,2)=(ZZ(2)-BV(5)*SCOMP(3,2))/BV(4)
	SCOMP(1,2)=SCOMP(2,1)
	
!	ZZ(1)=ZERO
!	ZZ(2)=ZERO
	ZZ(3)=1.0D0/BV(6)
	SCOMP(3,3)=ZZ(3)/BV(6)
	SCOMP(2,3)=SCOMP(3,2)
	SCOMP(1,3)=SCOMP(3,1)
	
!C  CALCULATE  SCOMP[GYF] AND STORE IT INTO DSIGMA
!C  DSIGMA=(/ZERO,ZERO,ZERO/)
	DSIGMA=ZERO
    DO K1=1,NTENS,1
	    DO K2=1,NTENS,1
	        DSIGMA(K2)=DSIGMA(K2)+SCOMP(K2,K1)*GYF(K1)
	    END DO
	END DO

!C  CALCULATE 1/K
    TT=HPF
	DO K1=1,NTENS,1
	    TT=TT+GYF(K1)*DSIGMA(K1)
	END DO

!C  UPDATE DDSDDE
    DO K1=1,NTENS,1
	    DO K2=K1,NTENS,1
	        TTB=SCOMP(K1,K2)-DSIGMA(K1)*DSIGMA(K2)/TT
	        DDSDDE(K1,K2)=TTB
	        DDSDDE(K2,K1)=TTB
	    END DO
	END DO
	DO K1=1,NTENS,1
	    DDSDDE(K1,K1)=SCOMP(K1,K1)-DSIGMA(K1)*DSIGMA(K1)/TT
	END DO
	DDSDDE(2,3)=SCOMP(2,3)-DSIGMA(2)*DSIGMA(3)/TT
	DDSDDE(1,3)=SCOMP(1,3)-DSIGMA(1)*DSIGMA(3)/TT
	DDSDDE(1,2)=SCOMP(1,2)-DSIGMA(1)*DSIGMA(2)/TT
	DDSDDE(2,1)=DDSDDE(1,2)
	DDSDDE(3,1)=DDSDDE(1,3)
	DDSDDE(3,2)=DDSDDE(2,3)
	
	
	deallocate(vUnits,vDegree)
    RETURN
    END SUBROUTINE  UMAT


!C**************************HARDENING***************************
!C*****NOTE:
!C*****THIS UMAT IDENTIFIES THE HARDENING SUB BY THE NAME 'KHARD'
!C*****(DEACTIVATE THE OTHER BY RENAMING)
     	
!C****: Swift (Power hardening law)
    SUBROUTINE   KHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ
    HF  = AAZ*((BBZ+EPBAR)**CCZ)
    HPF =  (CCZ/(BBZ+EPBAR))*HF
    RETURN
    END SUBROUTINE  KHARD

!C****: Voce (Exponential hardening law)
    SUBROUTINE  vcKHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ
    HF=BBZ/EXP(CCZ*EPBAR)
	HPF = CCZ*HF
    HF  = AAZ-HF
    RETURN
    END SUBROUTINE  vcKHARD

!C**********************************************************
	SUBROUTINE LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO,SCOMP, &
	                    AAZ,BBZ,CCZ,ZALPHA, &
	                   PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc, &
	                   nAllUnits,firstWeight) 

	IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    INTEGER::NTENS, nprops
	REAL(PREC),DIMENSION(NTENS)::STRESS,TDSIGMA,DSTRAN
	REAL(PREC),DIMENSION(NTENS,NTENS)::SCOMP
	REAL(PREC)::EPBAR,TDEPBAR,FZERO,AAZ,BBZ,CCZ,ZALPHA
    real(prec)::props(nprops)
	integer::nnLayers	
	real(prec):: vDegree(nnLayers)
	integer,dimension(nnLayers)::vUnits
	integer::nAllUnits, firstWeight
	real(prec)::homDegReciproc
	
!C     INTERNAL VARIABLES
    REAL(PREC),DIMENSION(NTENS)::TSIGMA,GYF  
	REAL(PREC)::HF,HPF,TEPBAR,YF,TT,FONE
	INTEGER::KK,JJ

    TSIGMA=STRESS+TDSIGMA
	TEPBAR=EPBAR+TDEPBAR
	
	CALL KHARD(HF,HPF,TEPBAR,AAZ,BBZ,CCZ)
	CALL GYFUNCTION(TSIGMA,NTENS,PROPS,NPROPS,vUnits,vDegree,nnLayers,homDegReciproc,&
	                nAllUnits,firstWeight,YF,GYF)
    FONE=(YF-HF)
	FONE=FONE*FONE
	DO KK=1,NTENS,1
	    TT=TDEPBAR*GYF(KK)-DSTRAN(KK)
	    DO JJ=1,NTENS,1
	        TT=TT+SCOMP(KK,JJ)*TDSIGMA(JJ)
	    END DO
	FONE=FONE+TT*TT
	END DO
	FONE=DSQRT(FONE)
	ZALPHA=1.0D0
	IF(FONE<=0.5D0*FZERO) RETURN	
	
!!	ZALPHA=0.75D0*(FZERO/(2.0D0*FONE))
    ZALPHA=0.375D0*(FZERO/FONE)	
    RETURN
    END SUBROUTINE LSEARCH	
!C*********************************************************************************

!CCCCCCC********** YIELD FUNCTION CALCULATIONS ***********************************
!CCCCCCC NOTE:   YFUNCTION RETURNS JUST YIELD FUNCTION VALUE
!CCCCCCC        GYFUNCTION RETURNS YIELD FUNCTION VALUE AND GRADIENT
!CCCCCCC        HYFUNCTION RETURNS YIELD FUNCTION VALUE, GRADIENT AND HESSIAN
!C
!CCCCCC NOTE: PLANE STRESS

    SUBROUTINE YFUNCTION(vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDegReciproc,&
	                     nAllUnits,firstWeight,Fval)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    integer::nvx,nvw,nvDegree
    real(prec)::vx(nvx), vw(nvw), vDegree(nvDegree)
    integer,dimension(nvDegree)::vUnits
    integer::nAllUnits, firstWeight
    real(prec)::homDegReciproc,Fval

	integer::KKw,JJ,KKunit, KKlayer, KKpreviousLayer, KKcurrentLayer
	real(prec)::ssum, vAllUnits(nAllUnits)

	KKpreviousLayer=0
	KKcurrentLayer=vUnits(1)
	vAllUnits(1:KKcurrentLayer)=vx
	KKw=firstWeight
	do KKlayer=2,nvDegree 
	  do KKunit=1,vUnits(KKlayer)
		ssum=0.0d0
		do JJ=1,vUnits(KKlayer-1)
		  ssum=ssum+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
		  KKw=KKw+1
		end do  
		vAllUnits(KKcurrentLayer+KKunit)=ssum**vDegree(KKlayer)
	  end do
	  KKpreviousLayer=KKcurrentLayer
	  KKcurrentLayer=KKcurrentLayer+vUnits(KKlayer)   
	end do
	Fval=0.0d0
	do JJ=1,vUnits(nvDegree)  !!!! Calculate output
	  Fval=Fval+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
	  KKw=KKw+1
	end do
	Fval=Fval**homDegReciproc !! output value	
	RETURN
	END SUBROUTINE YFUNCTION

!CC************************************************************************
    SUBROUTINE GYFUNCTION(vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDegReciproc,&
	                      nAllUnits,firstWeight,Fval,FGval)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

	integer::nvx,nvw,nvDegree
	real(prec)::vx(nvx), FGval(nvx), vw(nvw), vDegree(nvDegree)
	integer,dimension(nvDegree)::vUnits
	integer::nAllUnits, firstWeight
	real(prec)::homDegReciproc,Fval

	integer::KKw,JJ,KKunit, KKlayer, KKpreviousLayer, KKcurrentLayer
	real(prec)::ssum, gsum
	real(prec),dimension(nAllUnits):: vAllUnits
	real(prec),dimension(nAllUnits,3)::gAllUnits
	real(prec),dimension(nvDegree)::gDegree


	KKpreviousLayer=0
	KKcurrentLayer=vUnits(1)
	vAllUnits(1:KKcurrentLayer)=vx
	gAllUnits(1:KKcurrentLayer,:)=dble(0.0)
	do JJ=1,3
	gAllUnits(JJ,JJ)=dble(1.0)
	FGval(JJ)=0.0d0
	end do
	gDegree(:)=vDegree(:)-dble(1.0)

	KKw=firstWeight
	do KKlayer=2,nvDegree !!!total number of layers
	  do KKunit=1,vUnits(KKlayer)
		ssum=0.0d0; FGval(:)=0.0d0
		do JJ=1,vUnits(KKlayer-1)
		  ssum=ssum+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
		  FGval(1)=FGval(1)+gAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
		  FGval(2)=FGval(2)+gAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
		  FGval(3)=FGval(3)+gAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
		  KKw=KKw+1
		end do  
		vAllUnits(KKcurrentLayer+KKunit)=ssum**vDegree(KKlayer)
		gsum=vDegree(KKlayer)*ssum**gDegree(KKlayer)
		gAllUnits(KKcurrentLayer+KKunit,:)=gsum*FGval(:)
	  end do
	  KKpreviousLayer=KKcurrentLayer
	  KKcurrentLayer=KKcurrentLayer+vUnits(KKlayer)   
	end do
	ssum=0.0d0; FGval(:)=0.0d0
	do JJ=1,vUnits(nvDegree)  !!!! Calculate output and its gradient
	  ssum=ssum+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
	  FGval(1)=FGval(1)+gAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
	  FGval(2)=FGval(2)+gAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
	  FGval(3)=FGval(3)+gAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
	  KKw=KKw+1
	end do  
	Fval=ssum**homDegReciproc !! output value
	ssum=homDegReciproc*Fval/ssum
	FGval(:)=ssum*FGval(:)
	RETURN
	END SUBROUTINE GYFUNCTION

!CC***************************************************************************
!!Recall: nvx=3 for plane stress 
    SUBROUTINE HYFUNCTION(vx,nvx,vw,nvw,vUnits,vDegree,nvDegree,homDeg,homDegReciproc,&
	                      nAllUnits,firstWeight,Fval,FGval,FHval)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

	integer::nvx,nvw,nvDegree
	real(prec)::vDegree(nvDegree),vw(nvw),vx(nvx),Fval,FGval(nvx), FHval(nvx,nvx)
	integer,dimension(nvDegree)::vUnits
	integer::nAllUnits, firstWeight
	real(prec)::homDeg, homDegReciproc

	integer::KKw,JJ,KKunit, KKlayer, KKpreviousLayer, KKcurrentLayer
	real(prec)::ssum, gsum, hsum
	real(prec):: vAllUnits(nAllUnits), gAllUnits(nAllUnits,3), hAllUnits(nAllUnits,6)
	real(prec)::gDegree(nvDegree),hDegree(nvDegree)
	real(prec)::vFHval(6)

	KKpreviousLayer=0
	KKcurrentLayer=vUnits(1)
	vAllUnits(1:KKcurrentLayer)=vx
	gAllUnits(1:KKcurrentLayer,:)=dble(0.0)
	hAllUnits(1:KKcurrentLayer,:)=0.0d0
	do JJ=1,3
	gAllUnits(JJ,JJ)=dble(1.0)
	FGval(JJ)=0.0d0
	end do
	gDegree(:)=vDegree(:)-dble(1.0)
	hDegree(:)=gDegree(:)-1.0d0
	KKw=firstWeight
	do KKlayer=2,nvDegree !!!total number of layers
	  do KKunit=1,vUnits(KKlayer)
		ssum=0.0d0; FGval(:)=0.0d0; vFHval(:)=0.0d0
		do JJ=1,vUnits(KKlayer-1)
		  ssum=ssum+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
		  FGval(1)=FGval(1)+gAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
		  FGval(2)=FGval(2)+gAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
		  FGval(3)=FGval(3)+gAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
		  vFHval(1)=vFHval(1)+hAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
		  vFHval(2)=vFHval(2)+hAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
		  vFHval(3)=vFHval(3)+hAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
		  vFHval(4)=vFHval(4)+hAllUnits(KKpreviousLayer+JJ,4)*vw(KKw)
		  vFHval(5)=vFHval(5)+hAllUnits(KKpreviousLayer+JJ,5)*vw(KKw)
		  vFHval(6)=vFHval(6)+hAllUnits(KKpreviousLayer+JJ,6)*vw(KKw)
		  KKw=KKw+1
		end do  
		vAllUnits(KKcurrentLayer+KKunit)=ssum**vDegree(KKlayer)
		gsum=vDegree(KKlayer)*ssum**hDegree(KKlayer)
		hsum=gDegree(KKlayer)*gsum
		gsum=ssum*gsum
		gAllUnits(KKcurrentLayer+KKunit,:)=gsum*FGval(:)
		hAllUnits(KKcurrentLayer+KKunit,1)=gsum*vFHval(1)+hsum*FGval(1)*FGval(1)
		hAllUnits(KKcurrentLayer+KKunit,2)=gsum*vFHval(2)+hsum*FGval(1)*FGval(2)
		hAllUnits(KKcurrentLayer+KKunit,3)=gsum*vFHval(3)+hsum*FGval(1)*FGval(3)
		hAllUnits(KKcurrentLayer+KKunit,4)=gsum*vFHval(4)+hsum*FGval(2)*FGval(2)
		hAllUnits(KKcurrentLayer+KKunit,5)=gsum*vFHval(5)+hsum*FGval(2)*FGval(3)
		hAllUnits(KKcurrentLayer+KKunit,6)=gsum*vFHval(6)+hsum*FGval(3)*FGval(3)
	  end do
	  KKpreviousLayer=KKcurrentLayer
	  KKcurrentLayer=KKcurrentLayer+vUnits(KKlayer)   
	end do
	ssum=0.0d0; FGval(:)=0.0d0;  vFHval(:)=0.0d0
	do JJ=1,vUnits(nvDegree)  !!!! Calculate output and its gradient
	  ssum=ssum+vAllUnits(KKpreviousLayer+JJ)*vw(KKw)
	  FGval(1)=FGval(1)+gAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
	  FGval(2)=FGval(2)+gAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
	  FGval(3)=FGval(3)+gAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
	  vFHval(1)=vFHval(1)+hAllUnits(KKpreviousLayer+JJ,1)*vw(KKw)
	  vFHval(2)=vFHval(2)+hAllUnits(KKpreviousLayer+JJ,2)*vw(KKw)
	  vFHval(3)=vFHval(3)+hAllUnits(KKpreviousLayer+JJ,3)*vw(KKw)
	  vFHval(4)=vFHval(4)+hAllUnits(KKpreviousLayer+JJ,4)*vw(KKw)
	  vFHval(5)=vFHval(5)+hAllUnits(KKpreviousLayer+JJ,5)*vw(KKw)
	  vFHval(6)=vFHval(6)+hAllUnits(KKpreviousLayer+JJ,6)*vw(KKw)
	  KKw=KKw+1
	end do  
	Fval=ssum**homDegReciproc !! output value
	ssum=Fval/(homDeg*ssum)
	FGval(:)=ssum*FGval(:)
	hsum=(homDeg-1.0d0)/Fval
	FHval(1,1)=ssum*vFHval(1)-hsum*FGval(1)*FGval(1)
	FHval(2,1)=ssum*vFHval(2)-hsum*FGval(1)*FGval(2)
	FHval(3,1)=ssum*vFHval(3)-hsum*FGval(1)*FGval(3)
	FHVal(1,2)=FHval(2,1)
	FHval(2,2)=ssum*vFHval(4)-hsum*FGval(2)*FGval(2)
	FHval(3,2)=ssum*vFHval(5)-hsum*FGval(2)*FGval(3)
	FHval(1,3)=FHval(3,1)
	FHval(2,3)=FHval(3,2)
	FHval(3,3)=ssum*vFHval(6)-hsum*FGval(3)*FGval(3)
	
	RETURN
	END SUBROUTINE HYFUNCTION

!CC**************************** END OF NN-model IMPLEMENTATION
