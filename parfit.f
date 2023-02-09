c======================================================================c
c|..optional...Register...optional...Register...optional...Register...|c
c|--------------------------------------------------------------------|c
c| You have choosen to download the attached source code for our      |c
c| Fortran program betaFIT. I would appreciate it if you would please |c
c| go to the www address                                              |c
c|            http://scienide2.uwaterloo.ca/~rleroy/dParFit16/   and  |c
c| fill in the registration form there if you wish to be accessible   |c
c| so that I can send you possible future updates and/or corrections  |c
c| for this code.  This address list will be held securely by me and  |c
c| used for no other purpose................. Robert J. Le Roy .......|c
c|..Register...optional....Register...optional...Register...optional..|c
c======================================================================c

c*********************************************************************** 
        PROGRAM dParFit16
c***********************************************************************
c** Program "D(iatomic)Par(ameter)Fit" (DParFit) performs least-squares
c fits of a data set made up of any combination of MW, IR or electronic
c vibrational bands, PAS data, fluorescence series and/or Bv values, 
c involving one or more singlet or doublet-sigma electronic states and
c one or more isotopomers, to parameters defining the observed levels 
c of each state.  Those levels may be described by band constants {Gv, 
c Bv, Dv, etc.}, by Dunham expansions, by Near-Dissociation Expansions
c (NDE's), or by Tellinghuisen's mixed Dunham/NDE 'MXS' functions, and
c different expressions may be used for different states, while Lambda 
c or spin-rotation doubling is described by band constants or by Dunham-
c type expansions in (v+1/2).  Dunham, NDE or MXS function fits 
c automatically take account of normal (1'st order semiclassical) multi-
c isotopomer mass scaling, and allow for inclusion of atomic-mass-
c dependent Born-Oppenheimer and 1'st order JWKB breakdown corrections 
c (collectively called BOB corrections).
c++++++++++++++++++++ Version of  04 April 2016 ++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c              COPYRIGHT 2000-2016  by  Robert J. Le Roy               +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4051 +
c++++++++ by e-mail to: leroy@UWaterloo.ca , or write me at: +++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Uses least-squares subroutine NLLSSRR written by R.J. Le Roy +++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c* In any of these types of fits, centrifugal distortion constants,
c  and/or Lambda doubling parameters and/or BOB corrections and/or
c  the Gv & Bv parameters, for one or more of the electronic state may 
c  be held fixed, while a limited parameter set is varied.
c* This program always reports "sensitivities" of fitted parameters, which
c  indicate the numbers of significant digits which must be retained in
c  order to ensure predictions are in optimal agreement with experiment.
c* It will also perform automatic "Sequential Rounding & Refitting" [see
c  J.Mol.Spectrosc. 191, 223 (1998)] in order to yield a final parameter 
c  set involving the smallest possible number of significant digits.
c** If desired, it will also use a set of read-in constants to make 
c  predictions or to calculate deviations [calc.-obs.] for any chosen 
c  input data set involving diatomic singlet-singlet transitions.
c** Illustrative applications of this code are found in papers on HF/DF
c  [JMS 194,189 (1999)], GeO [JMS 194, 197 (1999)] and the coinage metal
c  hydrides [JCP 110, 11756 (1999)].
c=======================================================================
c** Dimensioning parameters intrinsic to the program are input through
c     PARAMETER statements in the file/data block  'arrsizes.h'. 
c** Parameters characterizing the problem and governing the fits are
c  read on  Channel-5  while the experimental data are read on Channel-4.
c** The principle output goes to  Channel-6  while higher output channel
c  numbers are used for secondary or more detailed/voluminous output.
c***********************************************************************
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc  INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
cc    INCLUDE 'PARMBLK.h'
c=======================================================================
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c=======================================================================
c
      CHARACTER*40 DATAFILE,MAKEPRED
      CHARACTER*43 FN4,FN6,FN7,FN8,FN9,FN10
      CHARACTER*20 NAMEPARM(NPARMX),WRITFILE
      CHARACTER*12 CTYPE(3)
      CHARACTER*3  CCDC(0:10)
      CHARACTER*2  CATOM
      CHARACTER*8  SPIN(4)
      INTEGER*4 lnblnk
c
      INTEGER  EFSEL(NSTATEMX),JTRUNC(NSTATEMX),NOWIDTHS,NLR(NSTATEMX),
     1 NSIG,NDAT(0:NVIBMX,NISTPMX,NSTATEMX),IFXP(NPARMX),
     2 NTVALL(0:NSTATEMX), I,I1,I2,IV,J,L,M,NSETS,NPARM,NDEORD,MMIN,
     3 MMAX,MQ0,MQM,MQMAX,PRINP,VMAXX,NRBC(NSTATEMX),
     3 ISTATE,ISOT,IROUND,JROUND,LPRINT,NDECOUNT,ATOM,ATOM2,
     4 CHARGE,LAMIN,IDUM1,IDUM2, NEWGv,NEWBv, MKPRED, ROBUST,CYCMAX
c
      REAL*8  PV(NPARMX),PU(NPARMX),PS(NPARMX),CM(NPARMX,NPARMX),
     1 PUSAV(NPARMX),PSSAV(NPARMX),CN(NSTATEMX),RM(NDUNMX), XX,XXP,YY,
     2 VDMV,VDMVP,Sw,SwLR,FNDE,DSE,PW,VPH,TSTPS,TSTPU,ZME,ZATOM,
     3 UCUTOFF,FDUM1
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
c
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
      DATA ZME/5.4857990945d-04/,CYCMAX/30/
      DATA CCDC/' Gv',' Bv','-Dv',' Hv',' Lv',' Mv',' Nv',' Ov',' Pv',
     1          ' Qv','CDC'/
      DATA CTYPE/' Outer Pade ',' Inner Pade ','Exponential '/
      DATA SPIN/' Singlet',' Doublet',' Triplet',' Quartet'/
      DATA MAKEPRED/'MAKEPRED                                '/
c Doublet can only be used for Sigma states;
c Triplet and Quartet NOT IN USE YET!   
c** Only special Data Types 0, -1 and -3  used in dParFit
      SLABL(-6)= '   '             !! data type not yet defined
      SLABL(-5)='VAC'              !! Accoustic Virial Coefficient
      SLABL(-4)='VIR'              !! Pressure Virial Coefficients  
      SLABL(-3)='BVV'              !! Experimental Bv values
      SLABL(-2)='WID'              !! tunneling level widths
      SLABL(-1)='PAS'              !!  Photo-Association binding energies
      SLABL(0)='FLS'               !!  fluorescence series
      NOWIDTHS= 1
c=======================================================================
c** Start by reading parameters describing the overall nature of the 
c   case and placing chosen restrictions on the data set to be used.
c
c  AN(1) & AN(2) are the integer atomic numbers identifying the atoms 
c          forming the molecule.
c
c** CHARGE (+/- integer) is the charge on the molecule (=0 for neutral).
c   If(CHARGE.ne.0) use Watson's(JMS 1980) charge-modified reduced mass.
c
c  NISTP   is the number of isotopomers to be simultaneously considered.
c
c  NSTATES  is the number of electronic states associated with the data
c        set to be analysed:  NSTATES = 1  for fits to IR/MW and/or 
c        fluorescence data for a single electronic state, while  
c        NSTATES > 1  for multi-state fits.
c        Upper states of fluorescence series NOT included in this count.
c
c  DATAFILE  is the (character variable) name of the file containing the
c     experimental data to be used in the fit.  If it is not located in
c     the current directory, the name 'DATAFILE' must include the 
c     relative path.  The variable name may (currently) consist of up to
c     40 characters.  READ ON A SEPARATE LINE!
c
c !! To make predictions using a completely specified set of parameters,
c      the input value of parameter DATAFILE must be  'MAKEPRED'
c
c  WRITFILE  is the (character-variable) name of the file to which the
c     output will be written.  Channel-6 outut goes to  WRITFILE.6,
c     channel-7 output to WRITFILE.7, channel-8 to WRITFILE.8, ... etc.
c     If not in the current directory, the name 'WRITFILE' must include the 
c     relative path.  The valiable name may (currently) consist of up to
c     20 characters, enclosed in single quotes, with no leading spaces.
c-----------------------------------------------------------------------
      READ(5,*)  AN(1), AN(2), CHARGE, NISTP, NSTATES
	READ(5,*)  DATAFILE
      READ(5,*)  WRITFILE
c-----------------------------------------------------------------------
c** These statements construct and define the names of output files 
c   associated with WRITE's to channels 6-10 used by the program.
      WRITE(FN6,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.6'
      WRITE(FN7,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.7'
      WRITE(FN8,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.8'
      WRITE(FN9,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.9'
      WRITE(FN10,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.10'
      OPEN(UNIT= 6, FILE= FN6)
      OPEN(UNIT= 7, FILE= FN7)
      OPEN(UNIT= 8, FILE= FN8)
      OPEN(UNIT= 9, FILE= FN9)
      OPEN(UNIT=10, FILE= FN10)
      MKPRED= 0
      IF(DATAFILE.EQ.MAKEPRED) THEN
          MKPRED= 1
          ENDIF
c=======================================================================
c  UCUTOFF   Neglect any input data with uncertainties > UCUTOFF (cm-1)
c
c  IROUND .ne. 0  causes "Sequential Rounding & Refitting" to be
c             performed, with each parameter being rounded at the
c            |IROUND|'th sig. digit of its local uncertainty.
c         = 0  simply stops after full convergence (without rounding).
c
c  ROBUST > 0  (integer) causes "Robust" least-squares weighting (as per
c              Watson [J.Mol.Spectrosc. 219, 326 (2003)] to be used
c         = 0  uses normal data weights  1/[uncertainty(i)]**2
c
c  LPRINT  specifies the level of printing inside NLLSSRR
c        if: =  0, no print except for failed convergence.
c             < 0  only converged, unrounded parameters, PU & PS's
c            >= 1  print converged parameters, PU & PS's
c            >= 2  also print parameter change each rounding step
c            >= 3  also indicate nature of convergence
c            >= 4  also print convergence tests on each cycle
c            >= 5  also parameters changes & uncertainties, each cycle
c
c  PRINP > 0  causes a summary of the input data to be printed before
c             the fitting starts.  Normally set =0.
c-----------------------------------------------------------------------
      READ(5,*)  UCUTOFF, IROUND, ROBUST, LPRINT, PRINP
c-----------------------------------------------------------------------
c
      I= 999
      IF(NISTP.LE.NISTPMX) THEN
          IF(CHARGE.NE.0) WRITE(6,600) NISTP,CHARGE
          IF(CHARGE.EQ.0) WRITE(6,600) NISTP
        ELSE
          WRITE(6,601) NISTP,NISTPMX
          STOP
        ENDIF
      WRITE(6,602)
      DO  ISOT= 1,NISTP
c** Read the mass numbers of the atoms in each of the isotopomers
c   MN(i,ISOT)  is the mass number for atom with atomic number AN(i) 
c       [NOTE: be sure order of MN values consistent with that of AN's].
c       Choosing it .ne.(value for some known isotope) of that species
c       causes the average atomic mass to be used.
c-----------------------------------------------------------------------
          READ(5,*) MN(1,ISOT), MN(2,ISOT)
c-----------------------------------------------------------------------
          I= MIN(I,MN(1,ISOT),MN(2,ISOT))
          CALL MASSES(AN(1),MN(1,ISOT),CATOM,IDUM1,IDUM2,ZMASS(1,ISOT),
     1                                                          FDUM1)
          IF(ISOT.EQ.1) NAME(1)= CATOM
          CALL MASSES(AN(2),MN(2,ISOT),CATOM,IDUM1,IDUM2,ZMASS(2,ISOT),
     1                                                          FDUM1)
          IF(ISOT.EQ.1) NAME(2)= CATOM
          ZMASS(3,ISOT)= (ZMASS(1,ISOT)*ZMASS(2,ISOT))/
     1                      (ZMASS(1,ISOT)+ ZMASS(2,ISOT)- CHARGE*ZME)
          WRITE(6,603) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                                           (ZMASS(J,ISOT),J=1,3)
          XX= ZMASS(3,1)/ZMASS(3,ISOT)
          DO  M=0, 9
              RMUP(M,ISOT)= XX**M
              ENDDO
          RSQMU(ISOT)= DSQRT(XX)
          XX= 1.d0
          DO  L= 0,NDUNMX
              RSQMUP(L,ISOT)= XX
              XX= XX*RSQMU(ISOT)
              ENDDO
          ENDDO
      IF(CHARGE.NE.0) WRITE(6,604) CHARGE
      IF(I.EQ.0) WRITE(6,605)
      WRITE(6,599) DATAFILE
  599 FORMAT(/' Use experimental data input file:  ',a30)
      IF(IROUND.NE.0) WRITE(6,606) IABS(IROUND)
      IF(IROUND.GT.0) WRITE(6,608)
      IF(IROUND.LT.0) WRITE(6,610) 
      IF(ROBUST.GT.0) THEN
          ROBUST= 1
          WRITE(6,596)
        ELSE
          WRITE(6,598)
        ENDIF
      NDECOUNT= 0
      DO  I= 1,NPARMX
          PV(I)= 0.d0
          ENDDO
c
c=======================================================================
c** Now ... loop over the NSTATES electronic states, reading parameters
c  characterizing those states and how to represent the data for each.
c=======================================================================
      DO 60 ISTATE= 1, NSTATES
c** For each of the electronic states s=ISTATE involved in the data set,
c  read parameters characterizing that state and the data to be used, 
c  and identifying the types of parameters used to characterize its term
c  values and whether they are to be fixed or fitted.
c================================
c  SLABL(s)  is a 2-character alphameric label enclosed in single quotes
c           to identify the electronic state; e.g., 'X0', 'A1', ... etc.
c
c  IOMEG(s) .GE.0  is electronic angular momentum of singlet state with
c                    projection quantum number  Lambda= IOMEG
c          < 0  for Sigma state with spin multiplicity  |IOMEG|= -IOMEG
c               [currently only coded for  IOMEG= -2 {doublet Sigma}]
c
c  IOMEG(s) .LT.0  if it indicates a doublet SIGMA electronic state
c          [may later introduce an additional read-in parameter MULTPLT(s)
c             to label other electronic state spin-multiplicities]
c
c  JTRUNC(s):  Omit from fit electronic state-s data with  J(s) > JTRUNC
c
c  EFSEL(s)  allows a user to consider data for:
c          * ONLY the e-parity levels of this state, if EFSEL > 0
c          * ONLY the f-parity levels of this state, if EFSEL < 0
c          * BOTH e- and f-parity levels of this state, if EFSEL = 0
c
c  VMIN(s)/VMAX(s):  Neglect data for electronic state-s vibrational  
c                    levels outside the range  VMIN(s)  to  VMAX(s). 
c---------------------------------------------------------------------
          READ(5,*) SLABL(ISTATE), IOMEG(ISTATE), JTRUNC(ISTATE),
     1                       EFSEL(ISTATE), VMIN(ISTATE), VMAX(ISTATE)
c---------------------------------------------------------------------
c** MULTPLT(s)  is the spin multiplicity of electronic state-s
c  Currently only allow Doublet states for IOMEG(s) < 0  and Singlet
c  states for  IOMEG(s) .ge. 0
          MULTPLT(ISTATE)= 1
          IF(IOMEG(ISTATE).LT.0) THEN
              MULTPLT(ISTATE)=IABS(IOMEG(ISTATE))
              efREF(ISTATE)= 0
              IF(MULTPLT(ISTATE).NE.1) THEN
                  WRITE(6,594) MULTPLT(ISTATE)
  594 FORMAT(/' *** INPUT ERROR *** program not coded to handle',i3,
     1  '-Sigma  spin multiplets')
                  STOP
                  ENDIF
              ENDIF
          IF((SLABL(ISTATE).EQ.'FLS').OR.(SLABL(ISTATE).EQ.'BVV').OR.
     1                                  (SLABL(ISTATE).EQ.'PAS')) THEN
              WRITE(6,607) ISTATE, SLABL(ISTATE)
              STOP
              ENDIF
c
c  NCDC(s)  denotes the number of centrifugal distortion constants
c           (CDC's) considered for electronic state-s.
c
c  NDEGv(s) & NDEBv(s), respectively, specify whether Gv & Bv for levels
c          of electronic state-s will be represented by:
c           (a)  band constants (Gv & Bv)  when    NDEXv(s) = -1
c           (b)  pure Dunham expansions    when    NDEXv(s) = 0 
c           (c)  pure NDE expressions      when    NDEXv(s) = 1
c           (d)  MXS mixed NDE/Dunham expressions when    NDEXv(s) = 2
c      NOTE:  require  NDEGv .ge. NDEBv  
c
c  NDECDC(s)  specifies whether CDC's for levels of electronic state-s
c          will be represented by:
c            (a)  band constants (-Dv, Hv, ...)  when    NDECDC(s) = -1
c            (b)  Dunham expansions              when    NDECDC(s) = 0
c            (c)  NDE expressions                when    NDECDC(s) = 1
c
c  IFXGv(s), IFXBv(s) & IFXCDC(s), respectively, specify whether 
c        constants determining Gv & Bv for electronic state-s will be:
c           (a)  held fixed at read-in values  (IFXYv > 0 )  , or
c           (b)  determined from the fits      (IFXYv.le.0).
c                                                                              
c  BOBORD(s)  indicates whether atomic-mass-dependent Born-Oppenheimer 
c           breakdown and higher-order JWKB correction (BOB) terms are
c           to be used for this electronic state:  
c    .ge.0  implies YES, and BOBORD is the highest order in [J(J+1)] used
c      < 0  implies such terms are NOT considered for this state.
c*         These constants are fitted or held fixed as defined by the 
c          Gv/Bv parameter IFXGv(s) & IFXBv(s). 
c-----------------------------------------------------------------------
          READ(5,*) NCDC(ISTATE), NDEGv(ISTATE), NDEBv(ISTATE),
     1                   NDECDC(ISTATE), IFXGv(ISTATE), IFXBv(ISTATE),
     2                   IFXCDC(ISTATE), BOBORD(ISTATE)
c-----------------------------------------------------------------------
          IF(IOMEG(ISTATE).GE.0) WRITE(6,613) SLABL(ISTATE), 
     1           SPIN(MULTPLT(ISTATE)), IOMEG(ISTATE), IOMEG(ISTATE)**2
          IF(IOMEG(ISTATE).LT.0) WRITE(6,613) SLABL(ISTATE), 
     1           SPIN(MULTPLT(ISTATE)), IOMEG(ISTATE)
c=======================================================================
c** Ensure internal consistency among read-in parameters .....
c=======================================================================
          IF(NCDC(ISTATE).GT.8) THEN
              WRITE(6,611) ISTATE,NCDC(ISTATE)
              NCDC(ISTATE)= 8
              ENDIF
          IF((NISTP.LE.1).AND.(BOBORD(ISTATE).GE.0)) THEN
              WRITE(6,612) ISTATE
              BOBORD(ISTATE)= -1
              ENDIF
          IF(NDEBv(ISTATE).GT.NDEGv(ISTATE)) THEN
c** Bv's for a given state cannot be represented by a "more
c      sophisticated" form than its  Gv's
              WRITE(6,614) SLABL(ISTATE)
              STOP
              ENDIF
          IF((NDECDC(ISTATE).GT.NDEBv(ISTATE)).AND.
     1                                      (NDECDC(ISTATE).LE.0))THEN
c** CDC's for a given state cannot be represented by a "more
c      sophisticated" form than its  Bv's UNLESS is if fixed NDE fx.
              WRITE(6,616) SLABL(ISTATE)
              STOP
              ENDIF
          IF((BOBORD(ISTATE).GE.0).AND.(NDEGv(ISTATE).EQ.-1)) THEN
c** If use Band Constants for vib or rot constants, cannot consider BOB
              WRITE(6,620) SLABL(ISTATE)
              BOBORD(ISTATE)= -1
              ENDIF
          IF((BOBORD(ISTATE).GE.1).AND.(NDEBv(ISTATE).EQ.-1)) THEN
              WRITE(6,619) SLABL(ISTATE)
              BOBORD(ISTATE)= 0
              ENDIF
          IF((BOBORD(ISTATE).GE.2).AND.(NDECDC(ISTATE).EQ.-1)) THEN
              WRITE(6,618) SLABL(ISTATE)
              BOBORD(ISTATE)= 1
              ENDIF
c
          IF((NDECDC(ISTATE).GT.0).AND.(IFXCDC(ISTATE).LE.0)) THEN
c** If CDC's for state s=ISTATE to be represented by NDE functions, 
c     program (currently) requires them to be held fixed.
              IFXCDC(ISTATE)= 1
              WRITE(6,622) SLABL(ISTATE)
              ENDIF 
c
c** Zero various expansion coefficients
          DO  L= 0, NDUNMX
              IF(L.GT.0) THEN
                  PM0(L,ISTATE)= 0.d0
                  QM0(L,ISTATE)= 0.d0
                  PM1(L,ISTATE)= 0.d0
                  RM(L)= 0.d0
                  ENDIF
              DO  M= 0, 9
                  YLM(L,M,ISTATE)= 0.d0
                  IF(M.GE.1) QLM(L,M,ISTATE)= 0.d0
                  DO  ATOM= 1,2
                      DELTA(ATOM,L,M,ISTATE)= 0.d0
                      ENDDO
                  ENDDO
              ENDDO
c** Zero Rotational Constant counters
          DO  ISOT= 1, NISTP
              DO  IV= VMIN(ISTATE),VMAX(ISTATE)
                  FITGV(IV,ISTATE,ISOT)= 0
                  NRC(IV,ISTATE,ISOT)= 0
                  NQC(IV,ISTATE,ISOT)= 0
                  NPAR(IV,ISTATE,ISOT)= 0
                  NQPAR(IV,ISTATE,ISOT)= 0
                  ENDDO
              ENDDO
c** Zero initial values of band constants and some counters
          DO  ISOT= 1,NISTP
              DO  IV= -1,VMAX(ISTATE)
                  IF(IV.GE.0) NDAT(IV,ISOT,ISTATE)= 0
                  DO   M= 0,9
                      ZK(M,IV,ISTATE,ISOT)= 0.d0
                      ENDDO
                  DO   M= 1,9
                      ZQ(M,IV,ISTATE,ISOT)= 0.d0
                      ENDDO
                  ENDDO
              ENDDO
          Te(ISTATE)= 0.d0
c
  596 FORMAT(/" Fit uses Watson's",' "Robust" data weighting [J.Mol/Spec
     1trosc. 219, 326 (2003)] '/20x,'1/[{unc(i)}^2 + {calc.-obs.}^2/3]')
  598 FORMAT(/' Fit uses standard  1/[uncertainty(i)]**2  data weighting
     1')
  600 FORMAT(' Input data for',I3,'  isotopomer(s)': ' of a species with
     1 net charge',SP,i3)
  601 FORMAT(' *** Array Dimensioning Problem:  NISTP=',i2,' > NISTPMX='
     1   , I3/10x,'Need to increase NISTPMX & Recompile')
  602 FORMAT(1x,16('**')/'    Isotopomer      Mass of atom-1   Mass of a
     1tom-2     Reduced mass'/1x,'--------------- ',
     2 3('   --------------'))
  603 FORMAT(1x,A2,'(',i3,')-',A2,'(',I3,') ',3(2x,F15.10))
  604 FORMAT(1x,67('-')/' Since this is an ion with charge',SP,i3,
     1  ", use Watson's charge-modified reduced mass.")
  605 FORMAT(2x,77('-')/2x,'Note that   (Mass Number) = 0   causes the a
     1verage atomic mass to be used.')
  606 FORMAT(/' Apply "Sequential Rounding & Refitting" at digit-',
     1  i1,' of the (local) parameter')
  607 FORMAT(/' *** ERROR *** State-',I2,'   Label ',A3,"  uses one of t
     1he reserved names"/26x,"'FLS', 'BVV' or 'PAS', so change its name!
     2 ")
  608 FORMAT(4x,'uncertainty, selecting remaining parameter with largest
     1 relative uncertainty')
  610 FORMAT(4x,'uncertainty, proceeding sequentially from the LAST para
     1meter to the FIRST.')
  611 FORMAT(/' *** CAUTION *** program array dimensions restrict   NCDC
     1 < 9'/17x,'so change read-in  NCDC(',i1,')  from',i3,'  to  8'/)
  612 FORMAT(/ ' *** INPUT ERROR *** CANNOT have BOB DELTA corrections f
     1or only ONE isotopomer!!'/10x,'... so reset BOBORD(',I1,') = -1')
  613 FORMAT(/' State ',A3,' is a',A8,' with  Omega=',i2/1x,7('***'):
     1 '  so rotational energies depend on powers of  [J(J+1)-',i2,']')
  614 FORMAT(" *** CAUTION - FIXUP State ",A3,"  INPUT ***  Only allow",
     1 '  NDEBv .le. NDEGv')
  616 FORMAT(" *** CAUTION - FIXUP State ",A3,"  INPUT ***  Only allow",
     1 '  NDECDC .le. NDEBv')
  618 FORMAT(/" *** INCONSISTENT INPUT:  When using band constants for C
     1DC's of State ",A3/5x,'CANNOT consider BOB corrections for them.')
  619 FORMAT(/" *** INCONSISTENT INPUT:  When using band constants for B
     1v's & CDC's of State ",A3/5x, 'CANNOT consider BOB corrections for
     2them.')
  620 FORMAT(/' *** INCONSISTENT INPUT:  Since use band constants for St
     1ate ',A3/5x,'CANNOT consider any BOB corrections for it.')
  622 FORMAT(" *** CAUTION - FIXUP State ",A3,"  INPUT ***  If CDC's rep
     1resented by NDE's, REQUIRE them to be held fixed")
c=======================================================================
c** Begin to explicitly identify nature of data representations 
c  to be used, and to read in any constants to be held fixed.
c=======================================================================
c** NOTE *** the program assumes the vibrational energy zero for each
c  state with Dunham or MXS vibrational energies [NGEGv=0 or 2] is the
c  hypothetical (v=-1/2,J=0) level of its 1-st isotopomer, while the 
c  reference energy for each NDE-described state [for which NGEGB>0] is
c  its asymptote.  Absolute energy of all upper states defined relative
c  to that for the first.  If vib. energies of state-1 treated with band
c  constants, absolute energy zero is lowest vib. level.
c** In the current code, if the Gv's and Bv's for a given state are
c  to be held fixed, then its CDC's, and (if BOBORD > 0) its BOB 
c  correction coefficients MUST also be fixed.
c
          IF(NDEGv(ISTATE).LE.-2) GOTO 60
          MMIN= -1
          MMAX= 0
          IF((IFXGv(ISTATE).GT.0).OR.(NDEGv(ISTATE).EQ.-1))
     1                                                  IFXD(ISTATE)=1
c
c=======================================================================
c** If Gv for this state represented by band constants (then Bv's and 
c  CDC's also represented the same way) .....
c=======================================================================
          IF(NDEGv(ISTATE).EQ.-1) THEN
              IF(IFXGv(ISTATE).GT.0) THEN
c** If the Gv's to be held fixed at read-in numerical values, then so
c  are Bv's and CDC's, and we read them in here.
                  IFXBv(ISTATE)= 1
                  IFXCDC(ISTATE)= 1
                  WRITE(6,630) SLABL(ISTATE)
                  MMAX= 1
                  IF(IFXCDC(ISTATE).GT.0) MMAX= NCDC(ISTATE)+ 1
c*  Read-in vibrational energies are assumed to be absolute energies,
c   (including the value of Te for that state).
c ... Looping over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... Looping over vibrational levels ...
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
c** For all vibrational levels  v=IV  from 0 to VMAX, read the 
c  vibrational energy  Gv = ZK(0,IV,ISTATE,ISOT),  and
c  Bv = ZK(1,IV,ISTATE,ISOT),  and  -Dv = ZK(2,IV,ISTATE,ISOT),  and
c  Hv = ZK(3,IV,ISTATE,ISOT),  and   Lv = ZK(4,IV,ISTATE,ISOT), ... etc.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** MAKE SURE THAT (!!) read-in values of leading CDC are -Dv, to make
c  sign convention for ZK(2,IV,ISTATE,ISOT) consistent with that for
c  other CDC's (s.th. +ve value makes +ve contribution to the energy!!)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
                          READ(5,*) I, (ZK(M,IV,ISTATE,ISOT), M=0, MMAX)
c-----------------------------------------------------------------------
                          IF(I.NE.IV) THEN
                              WRITE(6,632)  SLABL(ISTATE),I,IV
                              STOP
                              ENDIF
                          ENDDO
                      ENDDO
                ELSE
c======================================================
c** If fitting to band constants for all parameters ...
c======================================================
                  MMAX= 0
c ... Loop over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... and then loop over vibrational levels ...
                      DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c ... read variables specifying how many free parameters for each level:
c  FITGV(v,es,isot) = 0  for the lowest vibrational level (v) of each
c               separate "connected set";  otherwise set it  = 1
c  NRC(v,es,isot) = number of rotational constants (Bv, -Dv, Hv, Lv, ...)
c               to be fitted for that vibrational level ( =0 if no data)
c-----------------------------------------------------------------------
                          READ(5,*,END=200) I,FITGV(IV,ISTATE,ISOT),
     1                                             NRC(IV,ISTATE,ISOT)
c-----------------------------------------------------------------------
                          IF(I.NE.IV) THEN
                              WRITE(6,632)  SLABL(ISTATE),I,IV
                              STOP
                              ENDIF
                          IF(FITGV(IV,ISTATE,ISOT).LT.0) 
     1                                        FITGV(IV,ISTATE,ISOT)= 0
                          IF(FITGV(IV,ISTATE,ISOT).GT.0) 
     1                                        FITGV(IV,ISTATE,ISOT)= 1
c** Check consistency with "fixed constant" constraints read earlier
                          IF((ISTATE.EQ.1).AND.(IV.EQ.VMIN(ISTATE))
     1                      .AND.(ISOT.EQ.1)) FITGV(IV,ISTATE,ISOT)= 0
                          IF((IFXCDC(ISTATE).GT.0).AND.
     1                                     (NRC(IV,ISTATE,ISOT).GT.1))
     2                                          NRC(IV,ISTATE,ISOT)= 1
                          MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
                          ENDDO
                      ENDDO
                ENDIF
              IF(MMAX.EQ.1) GOTO 20
              IF(IOMEG(ISTATE).NE.0) GO TO 30
              GOTO 60
              ENDIF
c
          IF(NDEGv(ISTATE).EQ.2) THEN
c=======================================================================
c** If using Tellinghuisen-style Mixed Representations ... read 
c  VS(s) ... the isotopomer-1 v=v_S value where Dunham switches to NDE &
c  DVS(s) ... the width parameter for the switching function
c      (NOTE: be sure to put d0 on these REAL*8 numbers!)
c  IFXVS(s) .le.0  if Gv's fitted AND want to fit to VS; else IFXVS > 0.
c  IFXDVS(s) .le.0  if Gv's fitted AND want to fit to DVS; else  > 0.
c=======================================================================
              READ(5,*) VS(ISTATE), DVS(ISTATE), IFXVS(ISTATE),
     1                                                  IFXDVS(ISTATE)
c-----------------------------------------------------------------------
              IF(IFXGv(ISTATE).GT.0) THEN
                  IFXVS(ISTATE)= 1
                  IFXDVS(ISTATE)= 1
                  ENDIF
              XX= VS(ISTATE)+ 0.5d0
              DO  ISOT= 1, NISTP
                    VSISO(ISTATE,ISOT)= XX/RSQMU(ISOT) - 0.5d0
                    DVSISO(ISTATE,ISOT)= DVS(ISTATE)/RSQMU(ISOT)
                    ENDDO
              MMIN= 0
              MMAX= 0
              ENDIF
c
c=======================================================================
c** If Gv's defined by Dunham Y(l,m)'s ...
c=======================================================================
          IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).EQ.2)) THEN
c** L=LMAX(0,s) labels highest-order (highest-power) non-zero Y(L,0) in
c               the Dunham Gv expansion for electronic state 's'
c-----------------------------------------------------------------------
              READ(5,*) LMAX(0,ISTATE)
c-----------------------------------------------------------------------
cc            IF(LMAX(0,ISTATE).LE.0) LMAX(0,ISTATE)= -1
              IF(NDEGv(ISTATE).EQ.2) WRITE(6,634) SLABL(ISTATE),CCDC(0),
     1                LMAX(0,ISTATE),VS(ISTATE),VS(ISTATE),DVS(ISTATE)
              IF(IFXGv(ISTATE).GT.0) THEN
c
c** If Gv's held fixed at values defined by Dunham constants, read  
c  Te=T(v=-1/2) & Dunham Y(l,m) expansion coefficients for isotopomer-1
c  and define other coefficients using the normal isotope relations.
c  [NOTE ... EXclude Y(0,0)  ...  and for ISTATE=1,  Te = 0 !]
c** Program internally sets Y(0,0,s) = 'Te(s)' = T(v=-1/2)
c-----------------------------------------------------------------------
                  READ(5,*) Te(ISTATE)
                  IF(LMAX(0,ISTATE).GE.1) READ(5,*) (YLM(L,0,ISTATE), 
     1                                            L=1, LMAX(0,ISTATE))
c-----------------------------------------------------------------------
                  WRITE(6,635) SLABL(ISTATE), Te(ISTATE)
c ... first consider expansions for Gv and Bv ...
c** Set vib. energy zero for each state at  Te=T(-1/2)  for that state.
c  For lowest (s=1) state this defines absolute energy zero. 
                  YLM(0,0,ISTATE)= Te(ISTATE)
                  WRITE(6,636) SLABL(ISTATE),CCDC(0),(YLM(L,0,ISTATE),
     1                                            L= 1,LMAX(0,ISTATE))
                  MMIN= 0
                  MMAX= 0
                  ENDIF
              ENDIF
c
          IF(NDEGv(ISTATE).GE.1) THEN
c=======================================================================
c** If NDE functions used to represent the Gv's, or Gv & Bv, and/or 
c  the CDC's of state  s=ISTATE, read the limiting long-range (inverse)
c  power and (initial trial) values of D(limit), vD and Cn:
c=======================================================================
c  NLR(ISTATE) asymptotically-dominant (inverse) power associated with
c       the long-range potential defining the NDE function.
c  DLIMIT(ISTATE)  is absolute energy at the dissociation limit relative
c    to the reference energy for the first state considered.  For a
c    Dunham-described first state [NDEGv(1)=0], that reference energy 
c    is its  Gv(v=-1/2);  if the first state is described by an NDE or 
c    MXS Gv function [NDEGv(1) > 0], the reference energy is its 
c    asymptote [for which DLIMIT(1)=0].  If NDEGv(1)=-1, the reference
c    energy is the lowest vib. level of the first "connected set".  In
c    a multi-NDE/MXS case, the differences between the DLIMIT values are
c    simply the atomic level spacings.
c  VD(ISTATE)  non-integer effective vibrational index at dissociation
c              for isotopomer-1
c  CN(ISTATE)  is the coefficient of the asymptotically-dominant long-
c       range inverse-power potential term in units [(cm-1)*Angst**NLR]
c  NSIG  is the number of significant digits in CN, to be retained in
c       calculating rounded-off limiting theory coefficient X(Cn,n,zmu)
c-----------------------------------------------------------------------
              READ(5,*) NLR(ISTATE), DLIMIT(ISTATE), VD(ISTATE),
     1                                                CN(ISTATE), NSIG
c-----------------------------------------------------------------------
              NDECOUNT= NDECOUNT+ 1
              NUMNDE(ISTATE)= NDECOUNT
              NDEORD= 0
              IF(NDEBv(ISTATE).GT.0) NDEORD= 1
              IF(NDECDC(ISTATE).GT.0) NDEORD= NCDC(ISTATE)+1
              CALL NDEXM(ISTATE,SLABL(ISTATE),NLR(ISTATE),
     1           DLIMIT(ISTATE),VD(ISTATE),CN(ISTATE),NSIG,ZMASS(3,1),
     2                                     NDEORD,NISTP,RSQMU,PNDE,XM)
c** Define form of NDE expansion for Gv, and read initial trial/fixed
c    expansion parameters.
c* ITYPE  identifies type of NDE function for state 's':
c    (i) ITYPE=1  for an "outer" Pade expansion;  (ii)  ITYPE=2  for an
c    "inner" Pade, and  (iii) ITYPE=3  for an exponent polynomial NDE.
c* The leading non-zero contribution to the NP0-term numerator polynomial
c    is  (vD-v)**IP0 ,  while (for ITYPE=1 or 2) the corresponding
c    leading term in the NQ0-term denominator polynomial is  (vD-v)**IQ0.
c** NP0  is # terms in vib. exponent polynomial for case (iii)
c** PM0 & QM0  are vibrational numerator & denominator polynomial coeffts.
c-----------------------------------------------------------------------
              READ(5,*) ITYPE(ISTATE), NP0(ISTATE), NQ0(ISTATE), 
     1                                        IP0(ISTATE), IQ0(ISTATE)
              IF(NP0(ISTATE).GT.0) READ(5,*) (PM0(I,ISTATE), I=1,
     1                                                    NP0(ISTATE))
              IF(NQ0(ISTATE).GT.0) READ(5,*) (QM0(I,ISTATE), I=1,
     1                                                    NQ0(ISTATE))
c-----------------------------------------------------------------------
              IF(ITYPE(ISTATE).EQ.3) THEN
                  IP0(ISTATE)= 1
                  IQ0(ISTATE)= 1
                  ENDIF
c** If fitting NDE or mixed MXS function to Gv (IFXGv.le.0), specify 
c             whether or not  DLIMIT  and/or  vD  are to be held fixed:
c  If  IFXD or IFXVD > 0  hold  DLIMIT or VD, respectively, fixed at 
c          read-in value;  else it is to be varied in the fit.
c-----------------------------------------------------------------------
              IF(IFXGv(ISTATE).LE.0) READ(5,*)IFXD(ISTATE),IFXVD(ISTATE)
c-----------------------------------------------------------------------
c
c** Printout, for input NDE Gv parameters
              WRITE(6,638) SLABL(ISTATE),CCDC(0),NP0(ISTATE),
     1                   NQ0(ISTATE),CTYPE(ITYPE(ISTATE)),IP0(ISTATE),
     2                           IQ0(ISTATE),VD(ISTATE),DLIMIT(ISTATE)
              IF(NP0(ISTATE).GT.0) WRITE(6,640) (PM0(I,ISTATE),I=1,
     1                                                    NP0(ISTATE))
              IF(NQ0(ISTATE).GT.0) WRITE(6,642) (QM0(I,ISTATE),I=1,
     1                                                    NQ0(ISTATE))
              IP0(ISTATE)= IP0(ISTATE)- 1
              IQ0(ISTATE)= IQ0(ISTATE)- 1
              NEWGv= 1
ccc           NEWBv= 0
ccc           I= 0
ccc           CALL NDEDGB(ISTATE,1,NEWGv,NEWBv,RSQMU,I)
              ENDIF
c
c=======================================================================
c** If Bv's defined by Band Constants while Gv's are Not -- so are CDC's 
c=======================================================================
          IF((NDEBv(ISTATE).EQ.-1).AND.(NDEGv(ISTATE).GE.0)) THEN
              IF(IFXBv(ISTATE).GT.0) THEN
c** If the Bv's to be held fixed at read-in numerical values, then so
c  are CDC's, and we read them in here.
                  IFXCDC(ISTATE)= 1
                  WRITE(6,631) SLABL(ISTATE)
                  MMAX= 1
                  IF(IFXCDC(ISTATE).GT.0) MMAX= NCDC(ISTATE)+ 1
c ... Looping over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... Looping over vibrational levels ...
                      DO  I= VMIN(ISTATE), VMAX(ISTATE)
c** For all vibrational levels  v=IV  from 0 to VMAX, read in v and
c  Bv = ZK(1,IV,ISTATE,ISOT),  and  -Dv = ZK(2,IV,ISTATE,ISOT),  and
c  Hv = ZK(3,IV,ISTATE,ISOT),  and   Lv = ZK(4,IV,ISTATE,ISOT), ... etc.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** MAKE SURE THAT (!!) read-in values of leading CDC are -Dv, to make
c  sign convention for ZK(2,IV,ISTATE,ISOT) consistent with that for
c  other CDC's (s.th. +ve value makes +ve contribution to the energy!!)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
                          READ(5,*) IV, (ZK(M,I,ISTATE,ISOT), M=1, MMAX)
c-----------------------------------------------------------------------
                          IF(IV.NE.I) THEN
                              WRITE(6,632)  SLABL(ISTATE),I,IV
                              STOP
                              ENDIF
                          ENDDO
                      ENDDO
                ELSE
c** If FITTING to band constant Bv's, SAME for CDC's 
                  IFXCDC(ISTATE)= 0
                  IF(IOMEG(ISTATE).NE.0) IFXLD(ISTATE)= 0
                  MMAX= 1
c ... Loop over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... and then loop over vibrational levels ...
                      DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c ... read variables specifying how many free parameters for each level:
c  NRC(v,es,isot) = number of rotational constants (Bv, -Dv, Hv, Lv, ...)
c               to be fitted for vibrational level IV ( =0 if no data)
c-----------------------------------------------------------------------
                          READ(5,*,END=200) I ,NRC(IV,ISTATE,ISOT) 
c-----------------------------------------------------------------------
                          IF(IV.NE.I) THEN
                              WRITE(6,632)  SLABL(ISTATE),I,IV
                              STOP
                              ENDIF
                          FITGV(IV,ISTATE,ISOT)= 0
                          MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
                          ENDDO
                      ENDDO
                  IF(BOBORD(ISTATE).GT.0) BOBORD(ISTATE)= 0
                ENDIF
              GO TO 30
              ENDIF
c
c=======================================================================
c** If Bv's defined by Dunham Y(l,m)'s or by MXS Dunham+NDE functions:
c=======================================================================
          IF((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).GE.2)) THEN
c** L=LMAX(1,s) labels highest-order non-zero Y(L,1) in the Bv expansion
c-----------------------------------------------------------------------
              READ(5,*) LMAX(1,ISTATE)
c-----------------------------------------------------------------------
              IF(LMAX(1,ISTATE).LT.-1) LMAX(0,ISTATE)= -1
              IF(NDEBv(ISTATE).GE.2) WRITE(6,634) SLABL(ISTATE),CCDC(1),
     1                                       LMAX(1,ISTATE),VS(ISTATE)
              IF(IFXBv(ISTATE).GT.0) THEN
c** If Bv's held fixed at values defined by Dunham constants, read 
c  Dunham Y(l,m) expansion coefficients for isotopomer-1 and define
c  other coeffts. using the normal isotope relations.
c-----------------------------------------------------------------------
                  READ(5,*) (YLM(L,1,ISTATE), L=0, LMAX(1,ISTATE))
c-----------------------------------------------------------------------
c ... first consider expansions for Gv and Bv ...
c** Set vib. energy zero for each state at its isotopomer-1 (v=0,J=0)
c  level.  For lowest (s=1) state this defines absolute energy zero.
                  WRITE(6,636) SLABL(ISTATE),CCDC(1),(YLM(L,1,ISTATE),
     1                                            L= 0,LMAX(1,ISTATE))
                  IF(MMIN.LT.0) MMIN= 1
                  MMAX= 1
                  ENDIF
              ENDIF
c
          IF(NDEBv(ISTATE).GE.1) THEN
c=======================================================================
c** If Bv's represented by NDE or MXS functions ... read parameters
c  defining form of the NDE expansion and initial expansion parameters.
c* ITYPB  identifies type of NDE  Bv function for state 's':
c    (i) ITYPB=1  for an "outer" Pade expansion;  (ii)  ITYPB=2  for an
c    "inner" Pade, and  (iii) ITYPB=3  for an exponent polynomial NDE.
c* Leading non-zero contribution to the NP1-term numerator polynomial
c    is  (vD-v)**IP1 ,  while (for ITYPE=1 or 2) the corresponding
c    leading term in the NQ1-term denominator polynomial is (vD-v)**IQ1.
c* NP1 & NQ1  are # terms in rot. numerator & denom rational polynomials
c* PM1 & QM1  are vibrational numerator & denominator polynomial coeffts.
c-----------------------------------------------------------------------
              READ(5,*) ITYPB(ISTATE), NP1(ISTATE), NQ1(ISTATE),
     1                                        IP1(ISTATE), IQ1(ISTATE)
              IF(NP1(ISTATE).GT.0) READ(5,*) (PM1(I,ISTATE), I=1,
     1                                                    NP1(ISTATE))
              IF(NQ1(ISTATE).GT.0) READ(5,*) (QM1(I,ISTATE), I=1,
     1                                                    NQ1(ISTATE))
c-----------------------------------------------------------------------
              IF(ITYPB(ISTATE).EQ.3) THEN
                  IP1(ISTATE)= 1
                  IQ1(ISTATE)= 1
                  ENDIF
              IF(NP1(ISTATE).LT.0) NP1(ISTATE)= 0
              WRITE(6,638) SLABL(ISTATE),CCDC(1),NP1(ISTATE),
     1        NQ1(ISTATE),CTYPE(ITYPB(ISTATE)),IP1(ISTATE),IQ1(ISTATE)
              IF(NP1(ISTATE).GT.0) WRITE(6,640) (PM1(I,ISTATE),I=1,
     1                                                    NP1(ISTATE))
              IF(NQ1(ISTATE).GT.0) WRITE(6,642) (QM1(I,ISTATE),I=1,
     1                                                    NQ1(ISTATE))
              IP1(ISTATE)= IP1(ISTATE)- 1
              IQ1(ISTATE)= IQ1(ISTATE)- 1
              NEWBv= 1
              ENDIF
c** Call subroutine to generate predicted vibrational energies & Bv's
c [returned & held in COMMON/PARMBLK as ZK(0,v,s,ISTP) & ZK(1,v,s,ISTP)]
          IF((NDEGv(ISTATE).GT.0).OR.(NDEBv(ISTATE).GT.0))
     1        CALL NDEDGB(ISTATE,NISTP,NEWGv,NEWBv,RSQMU,VMAX(ISTATE))
c
c=======================================================================
c** Now ... begin consideration of CDC's
c=======================================================================
   20     IF((NDECDC(ISTATE).EQ.-1).AND.(IFXCDC(ISTATE).GT.0)
     1       .AND.((NDEBv(ISTATE).GE.0).OR.(IFXBv(ISTATE).LE.0))) THEN
c=======================================================================
c** If CDC's to be fixed at read-in band constant values (and Gv's & 
c  Bv's not fixed at read-in band constant values!), input CDC's here:
c  [NOTE: if Gv or Bv also fixed this way, CDC's already read in above!]
c=======================================================================
              WRITE(6,652) SLABL(ISTATE)
              MMAX= NCDC(ISTATE)+ 1
c ... Loop over isotopomers ...
              DO  ISOT= 1,NISTP
c ... Loop over vibrational levels ...
                  DO  I= VMIN(ISTATE), VMAX(ISTATE)
c** For all vibrational levels  v=IV  from 0 to VMAX, read
c  -Dv = ZK(2,IV,ISTATE,ISOT)  &  Hv = ZK(3,IV,ISTATE,ISOT) & ... etc.
c  [NOTE - input format assumed to include Gv & Bv, which are ignored!]
c%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** MAKE SURE THAT (!!) read-in values of leading CDC are -Dv, to make
c  sign convention for ZK(2,IV,ISTATE,ISOT) consistent with that for
c  other CDC's (s.th. +ve value makes +ve contribution to the energy!!)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
                      READ(5,*) IV, (ZK(M,I,ISTATE,ISOT), M= 2,MMAX)
c-----------------------------------------------------------------------
                      IF(IV.NE.I) THEN
                          WRITE(6,632)  SLABL(ISTATE),I,IV
                          STOP
                          ENDIF
                      ENDDO
                  ENDDO
              ENDIF              
c
          IF((NDECDC(ISTATE).EQ.-1).AND.(IFXCDC(ISTATE).LE.0)
     1       .AND.((NDEBv(ISTATE).GE.0))) THEN
c=======================================================================
c** If CDC's to be FITTED as band constants [but Gv's & Bv's are not]
c=======================================================================
              MMAX= 0
c ... Loop over isotopomers ...
              DO  ISOT= 1,NISTP
c ... and then loop over vibrational levels ...
                  DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c ... read variables specifying # free parameters for each level  v = I
c  NRC(v,es,isot) = total number of rotational constants for vib level v
c      of isotopomer isot in electronic state # es.  The number of CDC's 
c      {-Dv, Hv, Lv, ...} to be fitted is, of course  (NRC - 1).
c      [Set it =0 if insufficient data.]
c-----------------------------------------------------------------------
                      READ(5,*,END=200) I, NRC(IV,ISTATE,ISOT)
c-----------------------------------------------------------------------
                      IF(I.NE.IV) THEN
                          WRITE(6,632)  SLABL(ISTATE),I,IV
                          STOP
                          ENDIF
                      FITGV(IV,ISTATE,ISOT)= 0
                      MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
                      ENDDO
                  ENDDO
c** Ensure that BOB corrections not invoked for these Band Constants
              IF(BOBORD(ISTATE).GE.2) BOBORD(ISTATE)= 1
              GO TO 30
              ENDIF
c
          IF(NDECDC(ISTATE).EQ.0) THEN
c=======================================================================
c** If CDC's represented by Dunham expansions ...
c=======================================================================
c  L=LMAX(m,s) labels highest-order Y(L,m) in the Dunham expansion for 
c       CDC number 'm-1', for State-s [m=2 for -Dv, m=3 for Hv, ...]    
c-----------------------------------------------------------------------
              READ(5,*) (LMAX(M,ISTATE), M= 2, NCDC(ISTATE)+1)
c-----------------------------------------------------------------------
c** If CDC's are to be fixed at values defined by Dunham coefficients,
c  read in Ylm's for isotopomer-1 (defining others by mass scaling).
c** L=LMAX(m,s) labels highest-order Y(L,m) in expansion for ZK(m,v,s,1)
              IF(IFXCDC(ISTATE).GT.0) THEN
                  IF(MMIN.LT.0) MMIN= 2
                  MMAX= NCDC(ISTATE)+ 1
                  DO  M= 2,MMAX
                      IF(LMAX(M,ISTATE).GE.0) THEN
c-----------------------------------------------------------------------
                         READ(5,*) (YLM(L,M,ISTATE),L= 0,LMAX(M,ISTATE))
c-----------------------------------------------------------------------
                         WRITE(6,654) SLABL(ISTATE),M-1,CCDC(M),
     1                           (YLM(L,M,ISTATE),L= 0,LMAX(M,ISTATE))
                         ENDIF
                      ENDDO
                  ENDIF
              ENDIF
          IF(NDECDC(ISTATE).GT.0) THEN
c=======================================================================
c** If CDC's represented by NDE functions [ONLY allowed if CDC's fixed]:
c=======================================================================
              MMAX= NCDC(ISTATE)
              DO  M= 1,NCDC(ISTATE)
c ... drop NCDC to omit constants with no expansion terms
                  IF(LMAX(M,ISTATE).LT.0) MMAX= MMAX-1
                  ENDDO
              NCDC(ISTATE)= MMAX
              MMAX= NCDC(ISTATE)+ 1
c**  LMAX(m,s) No. of terms in exponential expansion for each CDC(m)
c**  RM(i)  are Isotope=1 exponent expansion coefficients for CDC(M-1)
c  (Note ... these coefficients not saved past this program segment)
c-----------------------------------------------------------------------
              READ(5,*) (LMAX(M,ISTATE), M= 2,MMAX)
              DO  M= 2,MMAX
                  IF(LMAX(M,ISTATE).GT.0) THEN
                      READ(5,*)  (RM(I),I= 1,LMAX(M,ISTATE))
c-----------------------------------------------------------------------
                      WRITE(6,656) SLABL(ISTATE),M-1,CCDC(M),(RM(I),
     1                                            I= 1,LMAX(M,ISTATE))
c** Now ... use read-in coefficients to generate actual fixed CDC values
                      DO ISOT= 1,NISTP
                          DO  IV= 0, VMAX(ISTATE) 
                              VDMV= (VD(ISTATE)- IV)*RSQMU(ISOT)
                              VDMVP= 1.d0
                              FNDE= 0.d0
                              DO  I= 1,LMAX(M,ISTATE)
                                  VDMVP= VDMVP*VDMV
                                  FNDE= FNDE+ VDMVP*RM(I)
                                  ENDDO
                              ZK(M,IV,ISTATE,ISOT)= XM(M,ISTATE,ISOT)*
     1                                 DEXP(FNDE)*VDMV**PNDE(M,ISTATE)
                              ENDDO
                          ENDDO
                      ENDIF
                  ENDDO
              ENDIF
c=======================================================================
c** If  IOMEG(s) > 0 ,  allow for a Lambda doubling parameter expansion.
c   If  IOMEG(s) < 0 ,  allow for a Gamma doubling parameter expansion.
c=======================================================================
   30     IF(IOMEG(ISTATE).NE.0) THEN
c=======================================================================
c  NLDMX(s) is the number of Lambda/Doublet-Sigma doubling parameters to
c     be considered for each vibrational levels.  For IOMEG > 0, leading
c     term proportional to [J(J+1)]^{IOMEG} and higher terms add powers
c     of [J(J+1)- IOMEG**2].  For IOMEG < 0 leading term proportional to
c     J (for e parity) or  (J+1) (for f parity), and higher terms add 
c     powers of [J(J+1)].
c  efREF(s)  identifies the `reference' (zero-shift) level in Lambda 
c      doubling as being the e sublevels [efREF= +1],  the f sublevels 
c      [for efREF= -1], or the mid-point [efREF= 0];  e.g., for efREF=-1
c      the f levels treated as unperturbed and the e levels shifted by
c      +q*[J(J+1)].  For  efREF= 0  attribute half of splitting to each.
c*    For 2\Sigma states set  efREF(s)= 0
c  NDELD(s) specifies form used for the Lambda doubling parameters:
c           a)  NDELD(s) < 0  ...  use band constant form 
c           b)  NDELD(s) = 0  ...  use conventional Dunham form
c   !!!!!   c)  NDELD(s) > 0  ...  use  (Bv)**2 times (v+1/2) polynomial
c   !!!!!    NDELD > 0  option NOT yet implemented!
c  IFXLD(s) specifies whether these for electronic state-s will be:
c           (a)  held fixed at read-in values  (IFXLD > 0 )  , or
c           (b)  determined from the fits      (IFXLD.le.0).
c-----------------------------------------------------------------------
              READ(5,*) NLDMX(ISTATE), efREF(ISTATE), NDELD(ISTATE),
     1                                                   IFXLD(ISTATE)
c-----------------------------------------------------------------------
              IF(NDELD(ISTATE).GT.0) NDELD(ISTATE)= 0
              MQ0= MAX0(0,IOMEG(ISTATE)-1)
              IF(IOMEG(ISTATE).LT.0) efREF(ISTATE)= 0
c=======================================================================
c** If using "band-constant" treatment of Lambda/Gamma-doubling constants
c for this state *** Loop over isotopomers and over vibrational levels
c to read parameters governing treatment of each level.
c=======================================================================
              IF((NLDMX(ISTATE).GT.0).AND.(NDELD(ISTATE).EQ.-1)) THEN
                  MQMAX= 0
                  DO  ISOT= 1, NISTP
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
c=======================================================================
c** ZQ(M, ...)  are FIXED values of the Lambda/Gamma doubling parameters 
c        for vibrational level I=IV of state ISTATE for isotopomer ISOT.
c  NQC(IV,es,isot) = number of Lambda/Gamma doubling constants to be
c    fitted for vibrational level v=IV [set =0 if no data or IOMEG(s)=0]
c-----------------------------------------------------------------------
                          IF(IFXLD(ISTATE).GT.0) READ(5,*) I, 
     1                    (ZQ(MQ0+M,I,ISTATE,ISOT),M= 1,NLDMX(ISTATE))
                          IF(IFXLD(ISTATE).LE.0) THEN
                              READ(5,*) I,NQC(IV,ISTATE,ISOT)
c-----------------------------------------------------------------------
                              MQMAX= MAX(MQMAX,NQC(IV,ISTATE,ISOT))
                              ENDIF
                          IF(I.NE.IV) THEN
                              WRITE(6,633) SLABL(ISTATE),IV,I
                              STOP
                              ENDIF 
                          ENDDO
                      ENDDO
                  IF(IFXLD(ISTATE).LE.0) NLDMX(ISTATE)= MQMAX
                  ENDIF
              IF((NLDMX(ISTATE).GT.0).AND.(NDELD(ISTATE).GE.0)) THEN
c=======================================================================
c  LDMAX(M,s) labels the highest-order (highest-power in v+1/2) non-zero
c       q(L,M) coefficient in the Dunham polynomial representation for   
c       the M-th order (power M in [J(J+1)]) Lambda/Gamma doubling coeffts.
c     If NO coefficients of this type, input:   LDMAX(M,s) < 0
c  If NDELD > 0  polynomial premultiplies  (Bv)**2 [NOT YET IMPLEMENTED]
c-----------------------------------------------------------------------
                  READ(5,*) (LDMAX(M+MQ0,ISTATE), M= 1,NLDMX(ISTATE))
c-----------------------------------------------------------------------
                  IF((IFXLD(ISTATE).GT.0).AND.(NLDMX(ISTATE).GT.0)) THEN
c** If Lambda/Gamma doubling coefficients are to be held fixed, read in 
c  the Dunham-type expansion coefficients to be used to define them.
                      DO  M= 1,NLDMX(ISTATE)
                          MQM= M+ MQ0
                          IF(LDMAX(MQM,ISTATE).GE.0) THEN
c-----------------------------------------------------------------------
                              READ(5,*) (QLM(L,MQM,ISTATE),L= 0,
     1                                              LDMAX(MQM,ISTATE))
c-----------------------------------------------------------------------
                              IF(IOMEG(ISTATE).GE.0) WRITE(6,658)
     1               SLABL(ISTATE),CCDC(MQM),(L,MQM,QLM(L,MQM,ISTATE),
     2                                         L= 0,LDMAX(MQM,ISTATE))
                              IF(IOMEG(ISTATE).LT.0) WRITE(6,659)
     1 SLABL(ISTATE),CCDC(M),(L,M,QLM(L,M,ISTATE),L= 0,LDMAX(M,ISTATE))
c** Now ... generate values of the fixed Lambda/Gamma doubling coefft.
c   for each vibrational level of each isotopomer.
                              DO  IV= VMIN(ISTATE),VMAX(ISTATE)
                                  DO  ISOT= 1,NISTP
                                      XX= (IV+0.5d0)*RSQMU(ISOT)
                                      YY= QLM(0,MQM,ISTATE)
                                      IF(LDMAX(MQM,ISTATE).GE.1) THEN
                                          XXP= 1.d0
                                          DO  L= 1,LDMAX(MQM,ISTATE)
                                             XXP= XXP*XX
                                             YY=YY+QLM(L,MQM,ISTATE)*XXP
                                             ENDDO
                                          ENDIF
                                      IF(IOMEG(ISTATE).GT.0)
     1          ZQ(MQM,IV,ISTATE,ISOT)= YY*(RMUP(1,ISOT))**(MQ0+MQM+1)
                                      IF(IOMEG(ISTATE).LT.0)
     1                           ZQ(M,IV,ISTATE,ISOT)= YY*RMUP(M,ISOT)
                                      ENDDO
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF
c
          IF(BOBORD(ISTATE).GE.0) THEN
c=======================================================================
c** If allowing for mass-dependent BOB/JWKB breakdown corrections ...
c=======================================================================
c* If fitting to BOB parameters, for ISTATE=1 must specify whether
c    (BOB00 > 0) or not (BOB00.le.0) the {0,0} coefficient is to be 
c    included when fitting to vibrational BOB corrections for this state
c    [Always included for other (ISTATE > 1) states!]
c* LAMAX(atom,i,s)  labels the highest power of the (v+1/2) expansion 
c      used to represent the atomic-mass-dependent BOB correction term
c      for  [J(J+1)]**i  (i=0 for pure vibration, etc.);  atom-1 for the
c      first atom of the pair & atom=2 for the second.
c      Power series in (v+1/2) starts at linear term (L=1) for i=0,
c      and at constant [L=0] for  i > 0  (rotational) terms
c*   Set LAMAX.lt.0  if no correction terms to be considered for this i.
c* For Homonuclear molecule, only read powers for the onee type of atom.
c-----------------------------------------------------------------------
              IF((IFXGv(ISTATE).LE.0).AND.(ISTATE.EQ.1)) READ(5,*) BOB00
              READ(5,*) (LAMAX(1,M,ISTATE), M= 0, BOBORD(ISTATE))
              IF(AN(1).NE.AN(2))
     1                READ(5,*) (LAMAX(2,M,ISTATE), M= 0,BOBORD(ISTATE)) 
c-----------------------------------------------------------------------
              IF(IFXGv(ISTATE).GT.0) THEN
c=======================================================================
c** If BOB correction terms are to be used and held fixed ... 
c ... read-em for each atom in turn
c=======================================================================
c** NOTE that leading (constant) vibrational correction coefficient 
c  normally fixed at 0.0 for lowest (ISTATE=1) state.
                  CATOM= NAME(1)
                  ATOM2= 2
                  IF(AN(1).EQ.AN(2)) ATOM2= 1
                  DO  ATOM= 1,ATOM2
                      DO  M= 0, BOBORD(ISTATE)
                          IF(LAMAX(ATOM,M,ISTATE).GE.0) THEN
c ... if included, read and write the fixed LeRoy-type BOB correction
c   expansion coefficients associated with [J(J+1)]**M, for each atom
c-----------------------------------------------------------------------
                              READ(5,*) (DELTA(ATOM,L,M,ISTATE),
     1                                  L= 0, LAMAX(ATOM,M,ISTATE))
c-----------------------------------------------------------------------
                              WRITE(6,660) CATOM,CCDC(M),SLABL(ISTATE),
     1                              (CATOM,L,M,DELTA(ATOM,L,M,ISTATE),
     2                                   L= 0,LAMAX(ATOM,M,ISTATE))
                              ENDIF
                          ENDDO
                      CATOM= NAME(2)
                      ENDDO
                  ENDIF
              ENDIF
c
c=======================================================================
c** As appropriate, now use any read-in Dunham coefficients to generate
c  desired fixed  ZK(M,v,ISTATE,ISOT) vib. & rot. energy band constants.
c=======================================================================
          IF(MMIN.GE.0) THEN 
              MMAX= 1
              IF((NDECDC(ISTATE).EQ.0).AND.(IFXCDC(ISTATE).GT.0)) 
     1                                           MMAX= NCDC(ISTATE)+ 1
              Sw= 1.d0
              SwLR= 0.d0
              DO  IV=0, VMAX(ISTATE)
                  DO  ISOT= 1,NISTP
                      XX= (IV+0.5d0)*RSQMU(ISOT)
                      IF((NDEGv(ISTATE).EQ.2).OR.(NDEBv(ISTATE).EQ.2)) 
     1                                                            THEN
                          SwLR= dexp((IV- VSISO(ISTATE,ISOT))/
     1                                            DVSISO(ISTATE,ISOT))
                          Sw= 1.d0/(1.d0+ SwLR)
                          SwLR= SwLR*Sw
                          ENDIF
                      DO  44 M= MMIN,MMAX
                          IF((M.EQ.1).AND.((IFXBv(ISTATE).LE.0)
     1                             .OR.(NDEBv(ISTATE).EQ.1))) GO TO 44
                          YY= 0.d0
                          IF(LMAX(M,ISTATE).GE.0) THEN
                              DO  L= LMAX(M,ISTATE),0,-1
                                  YY= YY*XX + YLM(L,M,ISTATE)
                                  ENDDO
                              ENDIF
                          YY= YY*RMUP(M,ISOT)
                          IF(((M.EQ.0).AND.(NDEGv(ISTATE).EQ.2)).OR.
     1                            ((M.EQ.1).AND.(NDEBv(ISTATE).GE.2)))
     2                           YY=  Sw*YY +SwLR*ZK(M,IV,ISTATE,ISOT)
                          ZK(M,IV,ISTATE,ISOT)= YY
   44                     CONTINUE
                      ENDDO
                  ENDDO
              ENDIF
c=======================================================================
c** If appropriate, add BOB correction contributions to fixed band constants
c=======================================================================
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).GT.0)) THEN
              ATOM2= 2
              IF(AN(1).EQ.AN(2)) ATOM2= 1
              LAMIN= 0
              DO  ISOT= 1,NISTP
                  DO  IV= 0, VMAX(ISTATE)
                      XX= (IV+0.5d0)*RSQMU(ISOT)
                      DO 50 M= 0,BOBORD(ISTATE)
                          IF(((M.EQ.0).AND.(IFXGv(ISTATE).LE.0)).OR.
     1                       ((M.EQ.1).AND.(IFXBv(ISTATE).LE.0)).OR.
     2                      ((M.GE.2).AND.(IFXCDC(ISTATE).LE.0)))GOTO 50
                          YY= 0.d0
                          ZATOM= 1.d0 - ZMASS(1,1)/ZMASS(1,ISOT)
                          IF(AN(1).EQ.AN(2)) ZATOM= ZATOM+ 1.d0 -
     1                                        ZMASS(2,1)/ZMASS(2,ISOT)
                          DO  ATOM= 1,ATOM2
                              IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                                  XXP= 1.d0
                                  IF(LAMIN.EQ.0) XXP= 1.d0/XX
                                  DO  L= LAMIN,LAMAX(ATOM,M,ISTATE)
                                      XXP= XXP*XX
                                      YY= YY+XXP*
     1                                    DELTA(ATOM,L,M,ISTATE)*ZATOM
                                      ENDDO
                                  ENDIF
                              ZATOM= 1.d0 - ZMASS(2,1)/ZMASS(2,ISOT)
                              ENDDO
                          ZK(M,IV,ISTATE,ISOT)= ZK(M,IV,ISTATE,ISOT) +
     1                                                 YY*RMUP(M,ISOT)
   50                     CONTINUE
                      ENDDO   
                  ENDDO
              ENDIF
          IF((IFXGv(ISTATE).GT.0).OR.(IFXBv(ISTATE).GT.0).OR.
     1                                     (IFXCDC(ISTATE).GT.0)) THEN
c=======================================================================
c** Print values of vib/rotation band constants held fixed in the fits.
c=======================================================================
              MMIN= 2
              IF(IFXBv(ISTATE).GT.0) MMIN= 1
              IF(IFXGv(ISTATE).GT.0) THEN
                  MMIN= 0
                  IF(IFXBv(ISTATE).LE.0) MMIN= -1
                  ENDIF
              MMAX= 0
              IF(IFXBv(ISTATE).GT.0) MMAX= 1
              IF(IFXCDC(ISTATE).GT.0) MMAX= NCDC(ISTATE)+ 1
              DO  ISOT= 1,NISTP
                  IF(MMIN.EQ.0) THEN
                      IF(MMAX.LE.1) WRITE(6,662) 
     1                       SLABL(ISTATE),NAME(1),MN(1,ISOT),NAME(2),
     2                       MN(2,ISOT),(CCDC(M),M= 0,MMAX)
                      IF(MMAX.GT.1) WRITE(6,662) 
     1                       SLABL(ISTATE),NAME(1),MN(1,ISOT),NAME(2),
     2                       MN(2,ISOT),(CCDC(M),M=0,1),(M-1,CCDC(M),
     3                                                       M=2,MMAX)
                      ENDIF
                  IF(MMIN.LT.0) WRITE(6,667) 
     1                       SLABL(ISTATE),NAME(1),MN(1,ISOT),NAME(2),
     2                       MN(2,ISOT),CCDC(0),(M-1,CCDC(M),M=2,MMAX)
                  IF(MMIN.EQ.1) WRITE(6,667) 
     1                       SLABL(ISTATE),NAME(1),MN(1,ISOT),NAME(2),
     2                       MN(2,ISOT),CCDC(1),(M-1,CCDC(M),M=2,MMAX)
                  IF(MMIN.EQ.2) WRITE(6,668) 
     1                       SLABL(ISTATE),NAME(1),MN(1,ISOT),NAME(2),
     2                       MN(2,ISOT), (M-1,CCDC(M),M=2,MMAX)
                  IF(MMIN.LT.0) WRITE(6,664) ('-',M=1,MMAX)
                  IF(MMIN.GE.0) WRITE(6,664) ('-',M=MMIN,MMAX)
                  DO  IV= VMIN(ISTATE),VMAX(ISTATE)
                      IF(MMIN.EQ.0) WRITE(6,666) IV,
     1                                (ZK(M,IV,ISTATE,ISOT),M= 0,MMAX)
                      IF(MMIN.LT.0) THEN
                          IF(MMAX.LT.2) WRITE(6,665) IV,
     1                                            ZK(0,IV,ISTATE,ISOT)
                          IF(MMAX.GE.2) WRITE(6,665) IV,
     1           ZK(0,IV,ISTATE,ISOT),(ZK(M,IV,ISTATE,ISOT),M= 2,MMAX)
                          ENDIF
                      IF(MMIN.EQ.1) WRITE(6,671) IV,
     1                                (ZK(M,IV,ISTATE,ISOT),M= 1,MMAX)
                      IF(MMIN.EQ.2) WRITE(6,672) IV,
     1                                (ZK(M,IV,ISTATE,ISOT),M= 2,MMAX)
                      ENDDO
                  ENDDO
              ENDIF
c
          IF((IOMEG(ISTATE).NE.0).AND.(IFXLD(ISTATE).GT.0).AND.
     1                                      (NLDMX(ISTATE).GT.0)) THEN
c=======================================================================
c** Print values of Lambda/Gamma Doubling constants held fixed in the fits.
c=======================================================================
              DO  ISOT= 1,NISTP
                  IF(IOMEG(ISTATE).GE.0) WRITE(6,663) SLABL(ISTATE),
     1 NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),(M+MQ0, M=1,NLDMX(ISTATE))
                  IF(IOMEG(ISTATE).LT.0) WRITE(6,6663) SLABL(ISTATE),
     1    NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),(M, M=1,NLDMX(ISTATE))
                  WRITE(6,664) ('---',M=1,NLDMX(ISTATE))
                  DO  IV= VMIN(ISTATE),VMAX(ISTATE)
                      WRITE(6,672) IV, (ZQ(M+MQ0,IV,ISTATE,ISOT), 
     1                                              M=1,NLDMX(ISTATE))
                      ENDDO
                  ENDDO
              IF(IOMEG(ISTATE).GT.0) THEN
                  IF(efREF(ISTATE).NE.0) WRITE(6,716) SLABL(ISTATE),
     1                                                   efREF(ISTATE)
                  IF(efREF(ISTATE).EQ.0) WRITE(6,718) SLABL(ISTATE)
                  ENDIF
              ENDIF
   60     CONTINUE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Call subroutine to input experimental data in specified band-by-band, 
c  format & do bookkeeping to document amounts of data or each type.
c
c For each "band", read in:  (i) upper/lower vibrational quantum numbers
c   VP & VPP,  (ii) a two-character electronic-state alphameric label
c   {enclosed in single quotes; e.g., 'X0' or 'A1'} for the upper
c   (LABLP) and lower (LABLP) state, and  (iii) integers NM1 & NM2 are
c   the mass numbers [corresponding to input atomic numbers AN(1) &
c   AN(2)] identifying the particular isotopomer.  Note that LABLP also
c   identifies the type of data in the 'band' or data-group (see below).
c
c** LABLP = LABLPP  and  VP = VPP  for a microwave band
c   LABLP = LABLPP  and  VP.ne.VPP  for an infrared band
c   LABLP = 'FLS'  identifies this data group/band as a fluorescence
c           series from a single emitting level into vibrational levels
c           of electronic state LABLPP.  In this case: VP is the quantum
c           number v' for the emitting level, while VPP is actually the
c           rotational quantum number J' for the emitting level and JP
c           [see below] the lower state vibrational quantum number v".
c   LABLP = 'BVV'  identifies this data group/band as a set of Bv values
c           for electronic state LABLPP.  In this case, parameters  VP
c           & VPP, and EFPP are dummy variables, as is JP and EFP [see
c           below],  JPP is actually the vibrational quantum number v",
c           EFPP the parity p", FREQ the Bv value & UFREQ its uncertainty
c** STOP reading when run out of bands OR when read-in VPP is negative
c-----------------------------------------------------------------------
cc    READ(4,*,END=20) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1,MN2
c-----------------------------------------------------------------------
cc      IF(VPP(IBAND).LT.0) GO TO 20
c** For each of the lines in a given band/series, read upper level
cc rotational quantum number (JP) and e/f parity (EFP= +1 for e, = -1 
cc for f, and =0 when e/f splitting not resolved), and lower level
cc rotational quantum number (JPP) and parity (EFPP), the transition
c  frequency  FREQ, and its uncertainty UFREQ.
c-----------------------------------------------------------------------
cc    5 READ(4,*) JP(COUNT), EFP(COUNT), JPP(COUNT), EFPP(COUNT), 
cc                                           FREQ(COUNT), UFREQ(COUNT)
c-----------------------------------------------------------------------
c** At end of a band, exit from implicit loop
cc      IF((JPP(COUNT).LT.0).OR.(JP(COUNT).LT.0)) GOTO 9
cc  -----------------------------------
cc  Sample IR band data of HF for the '.4' file:
cc  --------------------------------------------
cc  1 0  'X0' 'X0'  1 19             % VP VPP LABLP LABLPP MN1 MN2
cc  1 0  1 1   1 19                  % VP VPP IEP IEPP MN1 MN2
cc  8 1   9 1  266.0131002  0.005    % JP EFP JPP EFPP FREQ UFREQ
cc  9 1  10 1  265.8885896  0.003
cc 10 1  11 1  265.7716591  0.002
cc  .    .      .            . 
cc  .    .      .            . 
cc  [end of a band indicated by -ve JP and/or JPP value(s)]
cc -1 1  -1 1  -1.1         -1.1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      IF(MKPRED.LE.0) OPEN(UNIT= 4, STATUS= 'OLD', FILE= DATAFILE)
      IF(MKPRED.GT.0) THEN
          WRITE(FN4,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.4'
          OPEN(UNIT= 4, FILE= FN4)
          IF(UCUTOFF.LT.1.d0) UCUTOFF= 1.d0
          CALL MKPREDICT(NSTATES,NDAT)
          REWIND(4)
          ENDIF
      CALL READATA(NSTATES,NDEGV,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,
     1                                            NDAT,NOWIDTHS,PRINP)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Prepare and store powers of (v+1/2) for use in partial derivatives
      VMAXX= 0
      DO  ISTATE= 1,NSTATES
          VMAXX= MAX(VMAXX,VMAX(ISTATE))
          ENDDO
      IF(VMAXX.GT.NVIBMX) THEN
          WRITE(6,674) NVIBMX
          STOP
          ENDIF
      DO  IV= 0,VMAXX
          VPH= IV+0.5d0
          PW= 1.d0/VPH      
          DO L= 0,NDUNMX
              PW= PW*VPH     
              VPHPW(IV,L)= PW
              ENDDO
          ENDDO
c** Zero band origins of any fluorescence series.  To fix them, a USER
c   would need to change the code to add parameters & conditions!!]
      DO  IBAND= 1,NBANDMX
          ORIGIN(IBAND)= 0.d0
          ENDDO
  630 FORMAT(/" Fix State ",A3,"  Gv's,  Bv's  and  CDC's at read-in val
     1ues")
  631 FORMAT(/" Fix State ",A3,"  Bv's  and  CDC's at read-in values")
  632 FORMAT(/" **** ERROR  when reading State ",A3,"  Band Constants,",
     1 "   I=",i3,"   when  IV=v=",i3)
  633 FORMAT(/" **** ERROR  when reading State ",A3,"  Q(Lambda)'s.",
     1  "   I=",i3,"   when  IV=v=",i3)
  634 FORMAT(/" Represent State ",A3,1x,A3, "'s  by Tellinghuisen-type M
     1XS mixed representation:"/1x,9('=='),I3,' order Dunham for  v .le.
     2 VS=',F9.4,'  &  NDE for  v > VS':/9x,   'with switching function:
     3   Sw = 1/[1 + exp{(v -',F11.6,')/',F9.6,'}]')
  635 FORMAT(/' Fixed State ',A3,  "  Gv for Isot.-1 defined relative to
     1  T(v=-1/2)=",F12.5,' cm-1')
  636 FORMAT(' Fixed State ',A3,1x,A3,'  values defined by isotopomer-1'
     1 ,' Dunham coefficients:'/(1P4D20.12:))
  638 FORMAT(/' State ',A3,' NDE for',A3,' represented by  (NP=',I2,
     1 '/NQ=',I2,') ',A12,'NDE  in  (vD-v)'/1x,19('-'):'  with leading n
     2umerator and denominator powers',I3,' &',I3:/
     3 8x,'where for Isotopomer-1   vD=',F10.5,'   and   DLIMIT=',F12.4)
  640 FORMAT(5x,'Input numerator coefficients:  ',2(1PD21.13:)/
     1  (15X,3D21.13:))
  642 FORMAT(5x,'Input denominator coefficients:',2(1PD21.13:)/
     1  (15X,3D21.13:))
  652 FORMAT(/" Fix State ",A3,"  CDC's  at read-in values")
  654 FORMAT(/' Fixed State ',A3, '  CDC(',i1,')=',A3, '  values defined
     1 by isotopomer-1 Dunham coefficients:'/(1P4D20.12:))
  656 FORMAT(/" Fixed State ",A3,"  CDC(",I1,")=",A3,"'s defined by isot
     1opomer-1 NDE Exponent coefficients:"/(1P4D20.12:))
  658 FORMAT(/" Fixed State ",A3, '  "',A3,'"-Lambda Doubling defined by
     1 isotopomer-1 Dunham coeffts:'/(2x,3('   Q(',i2,',',i1,')=',
     2  1PD15.8:)))
  659 FORMAT(/" Fixed State ",A3,'  "',A3,'"-Gamma Doubling defined by i
     1sotopomer-1 Dunham coeffts.:'/(2x,3('   Q(',i2,',',i1,')=',
     2  1PD15.8:)))
  660 FORMAT(/' Fixed ',A2,' "',A3,'" BOB corrections for State ',A3,
     1  '  defined by expansion coefficients:'/
     2  (2x,2('   delta(',A2,',',i2,',',i2,')=',1PD19.11:)))
  662 FORMAT(/' State ',A3,'  Constants for  ',A2,'(',i3,')-',A2,'(',
     1  i3,')  held fixed at values:'/'   v',6x,A3:9x,A3:2x,
     2 6('   CDC(',i1,')=',A3:))
  663 FORMAT(/' State ',A3,'  Lambda Doubling Constants for  ',A2,'(',
     1  i3,')-',A2,'(',i3,')  fixed at values:'/
     2  '   v   ',6('     q(',i1,')    ':))
 6663 FORMAT(/' State ',A3,'  Gamma  Doubling Constants for  ',A2,'(',
     1  i3,')-',A2,'(',i3,')  fixed at values:'/
     2  '   v   ',6('     q(',i1,')    ':))
  664 FORMAT(2x,'--------------',A1:'----------',A1:
     1  6('------------',A1:))
  665 FORMAT(I4,F12.5,1P6D13.5)
  666 FORMAT(I4,F12.5,F12.8,1P6D13.5)
  667 FORMAT(/' State ',A3,'  Constants for  ',A2,'(',i3,')-',A2,'(',
     1  i3,')  held fixed at values:'/'   v',5x,A3:4x,6('   CDC(',
     2  i1,')=',A3:))
  668 FORMAT(/' State ',A3,'  Constants for  ',A2,'(',i3,')-',A2,'(',
     1 i3,')  held fixed at values:'/'   v',4('   CDC(',i1,')=',A3),
     2 2('  CDC(',i1,')=',A3:))
  671 FORMAT(I4,F12.8,1P4D13.5,2D12.4)
  672 FORMAT(I4,1P4D13.5,2D12.4)
  674 FORMAT(/' *** DIMENSIONING PROBLEM *** maximum vibrational range o
     1f input data EXCEEDS array dimension  NVIBMX=',i3/10x,'which is se
     2t in "included"  arrsizes.h  program file')
c=======================================================================
c=======================================================================
c** Loop over the NSTATES electronic states, preparing counters and 
c  labels for the free parameters to be determined by the fit.
c=======================================================================
c** For a global fit to data for one or more isotopomers with energies 
c  represented by (i) band constants, (ii) Dunham or (iii) NDE functions
c++ For each electronic state in turn, the fitted parameters are ordered
c++ in the following way:  +++++++++
c
c   1. For vibrational band-constant fits:  for each isotopomer, and for
c       each vibrational level of that isotopomer (for which there is
c       data), in turn, Gv, Bv, the free CDC's:  then SKIP to #5
c   2. The vibrational energy expansion parameters:
c       a) For pure Dunham (or MXS) expansions the parameters are (in
c          order) T(-1/2), Y_{1,0}, Y_{2,0}, Y_{3,0}, ... etc., where 
c          we always fix T(-1/2)=0  if  ISTATE=1  to define energy zero.
c       b) For mixed MXS functions, VS & DVS next (if they are fitted)
c       c) For NDE or mixed MXS functions, parameter order continues as:  
c          D=DLIMIT, vD, the PM's and then the QM's (where approptiate
c          for that ITYPE).
c   3) The Bv constant expansion parameters:
c       a) For rotational band-constant fits [when Gv not fitted by band
c          constants]:  for each isotopomer, and for each vibrational 
c          level of that isotopomer (for which there is rotational data),
c          in turn: Bv & the free CDC's:  then SKIP to #5
c       b) For Dunham expansions,  Y_{0,1}, Y_{1,1}, Y_{2,1}, ... etc.
c       c) For NDE functions, the exponent polynomial coefficients
c           p^1_1, p^1_2, p^1_3, p^1_4, ... etc.
c   4) The CDC Dunham expansion parameters (if free):
c       a) For fitted CDC band constants [when Bv's not fitted as band 
c         constants], the CDC's for each vib level of each isotopomer 
c       b) For Dunham expansions, Y_{0,m}, Y_{1,m}, ... etc., for  m=2
c         to NCDC(s)+1  [NOTE: fitting CDCs to NDE or MXS not implemented]
c   5) If  IOMEG > 0, and LAMBDA doubling considered, or ...
c      If  IOMEG < 0, and GAMMA  doubling considered, then
c       a) band-constant or Dunham-like doubling parameters for ISOT=1
c       b) band-constant or Dunham-like doubling parameters for ISOT=2
c       c) ... etc
c   6) If  BOBORD.ge.0  s.th. BOB correction  delta  coefficients used:
c      For each of atom-A and atom-B in turn, consider for  M=0  to
c      M=BOBORD  the  delta  expansion:  delta{atom,0,M},
c      delta{atom,1,M}, delta{atom,2,M}, delta{atom,3,M}, ... etc.
c+ After all electronic state constants taken into account, consider
c   7) the ORIGIN of each fluorescence series, in turn (all isotopomers)
c-----------------------------------------------------------------------
      ISOT= 1
      NPARM= 0
      NSETS= 0
c** IPSTATE(s) is a parameter counter s.th. [IPSTATE(s)+1] is the first 
c  free parameter for state (s)
      IPSTATE(1)= NPARM
      NTVALL(0)= 0
      DO 90 ISTATE= 1,NSTATES
c** Count parameters and prepare final printout ..
          IF(NDEGv(ISTATE).EQ.-2) THEN
c** If fitting to term values Tv(v,J,p) for this state, then ...
c=================================================================
              CALL TVSORT(ISTATE,NPARM,VMAX,NTVALL)
              NTVALL(0)= NTVALL(0)+ NTVALL(ISTATE)
              NDEBv(ISTATE)= -2
              NDECDC(ISTATE)= -2
              BOBORD(ISTATE)= -1
              GOTO 90
c!! Go to end of the ISTATES loop if using Term Values for this state
              ENDIF
          IF((NDEGv(ISTATE).EQ.-1).AND.(IFXGv(ISTATE).LE.0)) THEN
c** If fitting vib-rot levels of this state using band constants:
c=================================================================
              WRITE(6,676) SLABL(ISTATE),CCDC(0),CCDC(1),CCDC(10)
              WRITE(6,677)
              MMAX= 0
              DO  ISOT= 1,NISTP
                  DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c** Try to ensure that the # free parameters doesn't exceed # data for 
c   that v.  These checks are NOT rigorous, and problems yielding
c   underflows and  nan's  can still arise if too few independent data.
                      IF(NDAT(IV,ISOT,ISTATE).LT.(FITGV(IV,ISTATE,ISOT)
     1                                     + NRC(IV,ISTATE,ISOT))) THEN
                          IF(NDAT(IV,ISOT,ISTATE).GT.
     1                                     FITGV(IV,ISTATE,ISOT)) THEN
                              NRC(IV,ISTATE,ISOT)= NDAT(IV,ISOT,ISTATE)
     1                                          - FITGV(IV,ISTATE,ISOT)
                            ELSE
                              NRC(IV,ISTATE,ISOT)= 0
                              IF(NDAT(IV,ISOT,ISTATE).LT.
     1                 FITGV(IV,ISTATE,ISOT)) FITGV(IV,ISTATE,ISOT)= 0
                            ENDIF
                          WRITE(6,678) NDAT(IV,ISOT,ISTATE),
     1                    SLABL(ISTATE),ISOT,IV,FITGV(IV,ISTATE,ISOT),
     2                                             NRC(IV,ISTATE,ISOT)
                      ENDIF
                      MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
c** Use parameter counter NPAR to store the # free parameters preceeding
c  the first free vib/rot band constants for level IV of state ISTATE of
c  isotopomer ISOT, &  Accumulate total No. free parameters:  NPARM
                      NPAR(IV,ISTATE,ISOT)= NPARM
                      NPARM= NPARM+ FITGV(IV,ISTATE,ISOT)+ 
     1                                             NRC(IV,ISTATE,ISOT)
                      IF(FITGV(IV,ISTATE,ISOT).GT.0) 
     1                                    WRITE(7,761) CCDC(0),IV,ISOT
                      IF(NRC(IV,ISTATE,ISOT).GT.0) WRITE(7,761) 
     1                      (CCDC(M),IV,ISOT,M= 1,NRC(IV,ISTATE,ISOT))
                      IF((FITGV(IV,ISTATE,ISOT).LE.0).AND.
     1                                (NRC(IV,ISTATE,ISOT).GT.0)) THEN
                          NSETS= NSETS+ 1
                          WRITE(6,680) NSETS,SLABL(ISTATE),IV,ISOT
                          ENDIF
                      IF(FITGV(IV,ISTATE,ISOT).GT.0) WRITE(6,682) 
     1                                     IV,ISOT,NRC(IV,ISTATE,ISOT)
                      IF(FITGV(IV,ISTATE,ISOT).LE.0) WRITE(6,683) 
     1                                     IV,ISOT,NRC(IV,ISTATE,ISOT)
                      ENDDO
                  ENDDO
              NRBC(ISTATE)= NPARM- IPSTATE(ISTATE)
              NCDC(ISTATE)= MMAX-1
              GO TO 70
c====end of preparation for vib/rot-band-constant fit to this state=====
              ENDIF
c
          IF(IFXGv(ISTATE).LE.0) THEN
c=======================================================================
c** If Gv's are to be fitted to Dunham, NDE or mixed MXS functions ...
c=======================================================================
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).EQ.2)) THEN
c** If representing Gv for this state by fitted Dunham or MXS function:
c ... for state-1, absolute zero of energy fixed by setting  Te = 0 
                  IF(ISTATE.EQ.1) THEN
                      Te(ISTATE)= 0.d0
                      WRITE(6,686) SLABL(ISTATE)
                    ELSE
c ... for ISTATE > 1 ,  Te  is always fitted for Dunham or MXS  Gv's
                      NPARM= NPARM+ 1
                      PV(NPARM)= Te(ISTATE)
                      WRITE(7,765) SLABL(ISTATE)
                    ENDIF
                  ENDIF
                  IF(NDEGv(ISTATE).EQ.0) WRITE(6,688) SLABL(ISTATE),
     1                                          CCDC(0),LMAX(0,ISTATE)
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).EQ.2)) THEN
                  IF(LMAX(0,ISTATE).GT.0) THEN
                      NPARM= NPARM+ LMAX(0,ISTATE)
                      WRITE(7,766) (L,0,L= 1,LMAX(0,ISTATE))
                      ENDIF
                  ENDIF
c
              IF(NDEGv(ISTATE).GT.0) THEN
c** If representing Gv for this state by fitted NDE or MXS function:
c=======================================================================
                  IF(NDEGv(ISTATE).EQ.2) THEN
                      WRITE(6,689) SLABL(ISTATE),CCDC(0),LMAX(0,ISTATE),
     1                               VS(ISTATE),VS(ISTATE),DVS(ISTATE)
                      IF(IFXVS(ISTATE).LE.0) THEN
                          NPARM= NPARM+ 1
                          PV(NPARM)= VS(ISTATE)
                          WRITE(6,690)
                          WRITE(7,768) SLABL(ISTATE)
                          ENDIF
                      IF(IFXDVS(ISTATE).LE.0) THEN
                          NPARM= NPARM+ 1
                          PV(NPARM)= DVS(ISTATE)
                          WRITE(6,691)
                          WRITE(7,769) SLABL(ISTATE)
                          ENDIF
                      ENDIF
                  IF(ISTATE.EQ.1) THEN
                      IF(NDEGv(ISTATE).EQ.1) THEN
c** If State-1 fitted to (or fixed by) a pure NDE function, absolute 
c  zero of energy defined by fixing its asymptote at read-in value
                          IFXD(ISTATE)= 1
                          ENDIF
                    ELSEIF(NUMNDE(ISTATE).GT.1) THEN
c** For second or higher NDE-defined state, dissociation limit MUST be
c  fixed by the known atomic limit spacings (so override read-in IFXD).
                      IFXD(ISTATE)= 1
                    ENDIF
                  IF(IFXD(ISTATE).GT.0) THEN
                      WRITE(6,692) SLABL(ISTATE),DLIMIT(ISTATE)
                    ELSE
                      WRITE(6,694) SLABL(ISTATE)
                      NPARM= NPARM+1
                      PV(NPARM)= DLIMIT(ISTATE)
                      WRITE(7,770) SLABL(ISTATE)
                    ENDIF
                  IF(IFXVD(ISTATE).LE.0) THEN
                      NPARM= NPARM+ 1
                      PV(NPARM)= VD(ISTATE)
                      WRITE(7,771) SLABL(ISTATE)
                      ENDIF
c** Accumulate parameter count
                  WRITE(6,698) SLABL(ISTATE),CCDC(0),NP0(ISTATE),
     1                 NQ0(ISTATE),CTYPE(ITYPE(ISTATE)),IP0(ISTATE)+1,
     2                                        IQ0(ISTATE)+1,VD(ISTATE)
                  IF(NP0(ISTATE).GT.0) THEN
                      DO  I= 1,NP0(ISTATE)
                          NPARM= NPARM+ 1
                          PV(NPARM)= PM0(I,ISTATE)
                          ENDDO
                      WRITE(6,640) (PM0(I,ISTATE),I=1,NP0(ISTATE))
                      WRITE(7,773) (0,I+IP0(ISTATE),I= 1,NP0(ISTATE))
                      ENDIF
                  IF(NQ0(ISTATE).GT.0) THEN
                      DO  I= 1,NQ0(ISTATE)
                          NPARM= NPARM+ 1
                          PV(NPARM)= QM0(I,ISTATE)
                          ENDDO
                      WRITE(6,642) (QM0(I,ISTATE),I=1,NQ0(ISTATE))
                      WRITE(7,775) (0,I+IQ0(ISTATE),I= 1,NQ0(ISTATE))
                      ENDIF
                  ENDIF
              ENDIF
c
          IF(IFXBv(ISTATE).LE.0) THEN
c=======================================================================
c** If Bv's are to be fitted ..................
c=======================================================================
              IF((NDEBv(ISTATE).LT.0).AND.(NDEGv(ISTATE).GE.0)) THEN
c-----------------------------------------------------------------------
c** If fit Bv's (& hence CDC's) as band constants while Gv's treated as
c   Dunham, NDE or MXS expansions, then ...
c-----------------------------------------------------------------------
                  WRITE(6,676) SLABL(ISTATE),CCDC(1),CCDC(10)
                  WRITE(6,677)
                  MMAX= 0
c** Set parameter counter before beginning with band constants ...
                  NEBC(ISTATE)= NPARM
                  DO  ISOT= 1,NISTP
                      DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c** Try to ensure that the # free parameters doesn't exceed # data for
c   that v.  These checks are NOT rigorous, and problems yielding
c   underflows and  nan's  can still arise if too few independent data.
                          IF(NDAT(IV,ISOT,ISTATE).LT.
     1                                       NRC(IV,ISTATE,ISOT)) THEN
                              NRC(IV,ISTATE,ISOT)= 
     1                                            NDAT(IV,ISOT,ISTATE)
                              WRITE(6,678) NDAT(IV,ISOT,ISTATE),
     1                    SLABL(ISTATE),ISOT,IV,FITGV(IV,ISTATE,ISOT),
     2                                             NRC(IV,ISTATE,ISOT)
                              ENDIF
                          MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
c** Use parameter counter NPAR to store the # free parameters preceeding
c  the free rotational band constants for level IV of state ISTATE of 
c  isotopomer ISOT, &  Accumulate total No. free parameters:  NPARM
                          NPAR(IV,ISTATE,ISOT)= NPARM
                          IF(NRC(IV,ISTATE,ISOT).GT.0) THEN 
                              NPARM= NPARM+ NRC(IV,ISTATE,ISOT)
                              WRITE(7,761) (CCDC(M),IV,ISOT,M= 
     1                                          1,NRC(IV,ISTATE,ISOT))
                              ENDIF
                          IF(NRC(IV,ISTATE,ISOT).GT.1) WRITE(6,683)
     1                                   IV,ISOT,NRC(IV,ISTATE,ISOT)
                          ENDDO
                      ENDDO
c** Set no. fitted Rotational band constants for this state
                  NRBC(ISTATE)= NPARM- NEBC(ISTATE)
                  NEBC(ISTATE)= NPARM
                  NCDC(ISTATE)= MMAX-1
                  GO TO 70
c= end of preparation for all-rotational band-constant fit to this state
                  ENDIF
              IF((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).GE.2)) THEN
c** If representing Bv for this state by a fitted Dunham or MXS function
c=======================================================================
                  IF(NDEBv(ISTATE).EQ.0) 
     1               WRITE(6,688) SLABL(ISTATE),CCDC(1),LMAX(1,ISTATE)
                  NPARM= NPARM+ LMAX(1,ISTATE)+ 1
                  IF(LMAX(1,ISTATE).GE.0) 
     1                          WRITE(7,766) (L,1,L= 0,LMAX(1,ISTATE))
                  ENDIF
c
c** If representing Bv for this state by fitted NDE or MXS functions:
c=======================================================================
              IF(NDEBv(ISTATE).GE.1) THEN
                  IF(NDEBv(ISTATE).GE.2) WRITE(6,689) SLABL(ISTATE),
     1                               CCDC(1),LMAX(1,ISTATE),VS(ISTATE)
                  WRITE(6,698) SLABL(ISTATE),CCDC(1),NP1(ISTATE),
     1    NQ1(ISTATE),CTYPE(ITYPB(ISTATE)),IP1(ISTATE)+1,IQ1(ISTATE)+1
                  IF(NP1(ISTATE).GT.0) THEN
                      DO  I= 1,NP1(ISTATE)
                          NPARM= NPARM+ 1
                          PV(NPARM)= PM1(I,ISTATE)
                          ENDDO
                      WRITE(6,640) (PM1(I,ISTATE),I=1,NP1(ISTATE))
                      WRITE(7,773) (1,I,I= 1,NP1(ISTATE)) 
                      ENDIF
                  IF(NQ1(ISTATE).GT.0) THEN
                      DO  I= 1,NQ1(ISTATE)
                          NPARM= NPARM+ 1
                          PV(NPARM)= QM1(I,ISTATE)
                          ENDDO
                      WRITE(6,642) (QM1(I,ISTATE),I=1,NQ1(ISTATE))
                      WRITE(7,775) (1,I+IQ1(ISTATE),I= 1,NQ1(ISTATE))
                      ENDIF
                  ENDIF
              ENDIF
c
c=================================
c** If fitting to CDC's .....
c=================================
          IF((IFXCDC(ISTATE).LE.0).AND.(NDECDC(ISTATE).EQ.-1)
     1                                 .AND.(NDEBv(ISTATE).GE.0)) THEN
c** If fit CDC's as band constants while & Bv's (& Gv's) represented by
c   Dunham, NDE or MXS functions, then ...
              WRITE(6,676) SLABL(ISTATE),CCDC(10)
              WRITE(6,677)
              MMAX= 0
c** Set parameter counter before beginning with band constants ...
              NEBC(ISTATE)= NPARM
              DO  ISOT= 1,NISTP
                  DO  IV= VMIN(ISTATE),VMAX(ISTATE)
c** Try to ensure that the # free parameters doesn't exceed # data for
c   that v.  These checks are NOT rigorous, and problems yielding
c   underflows and  nan's  can still arise if too few independent data.
                      IF(NDAT(IV,ISOT,ISTATE).LT.NRC(IV,ISTATE,ISOT))
     1                                                              THEN
                          NRC(IV,ISTATE,ISOT)= NDAT(IV,ISOT,ISTATE)
                          WRITE(6,678) NDAT(IV,ISOT,ISTATE),
     1                    SLABL(ISTATE),ISOT,IV,FITGV(IV,ISTATE,ISOT),
     2                                             NRC(IV,ISTATE,ISOT)
                          ENDIF
                      MMAX= MAX(MMAX,NRC(IV,ISTATE,ISOT))
c** Use parameter counter NPAR to store the # free parameters preceeding
c  the free CDC band constants for level IV of state ISTATE of 
c  isotopomer ISOT, &  Accumulate total No. free parameters:  NPARM
                      NPAR(IV,ISTATE,ISOT)= NPARM
                      IF(NRC(IV,ISTATE,ISOT).GT.1) THEN 
                          NPARM= NPARM+ NRC(IV,ISTATE,ISOT) - 1
                          WRITE(7,761) (CCDC(M),IV,ISOT,M= 
     1                                          2,NRC(IV,ISTATE,ISOT))
                          ENDIF
                      IF(NRC(IV,ISTATE,ISOT).GT.1) WRITE(6,681)
     1                                   IV,ISOT,NRC(IV,ISTATE,ISOT)-1
                      ENDDO
                  ENDDO
c** Set number of fitted CDC band constants for this state
              NRBC(ISTATE)= NPARM- NEBC(ISTATE)
              NEBC(ISTATE)= NPARM
              NCDC(ISTATE)= MMAX-1
              GO TO 70
c=======end of preparation for CDC band-constant fit to this state======
              ENDIF
          IF((NDECDC(ISTATE).EQ.0).AND.(IFXCDC(ISTATE).LE.0)) THEN
c** If representing CDC's for this state by fitted Dunham functions
c=======================================================================
              WRITE(6,700) SLABL(ISTATE),(M-1,CCDC(M),
     1                             LMAX(M,ISTATE),M= 2,NCDC(ISTATE)+1)
              DO  M= 2,NCDC(ISTATE)+1
                  IF(LMAX(M,ISTATE).GE.0) THEN
                      NPARM= NPARM+ LMAX(M,ISTATE)+ 1
                      IF(LMAX(M,ISTATE).GE.0) THEN
                          WRITE(7,766) (L,M,L= 0,LMAX(M,ISTATE))
                          ENDIF
                      ENDIF
                  ENDDO
              ENDIF
c
   70     IF((IOMEG(ISTATE).NE.0).AND.(IFXLD(ISTATE).LE.0)
     1                                 .AND.(NLDMX(ISTATE).GT.0)) THEN
c=======================================================================
c** If Omega.ne.0  and Lambda or Gamma doubling constants being fitted 
c=======================================================================
              MQ0= MAX0(0,IOMEG(ISTATE)- 1)
              IF(NDELD(ISTATE).LT.0) THEN
c** If doubling constants are to be fitted using band constants ... 
c=======================================================================
                  DO  ISOT= 1, NISTP
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
c** Try to ensure that the # free parameters doesn't exceed # data for
c   that v.  These checks are NOT rigorous, and problems yielding
c   underflows and  nan's  can still arise if too few independent data.
                          IF(NDAT(IV,ISOT,ISTATE) .LT.
     1                 (NRC(IV,ISTATE,ISOT)+NQC(IV,ISTATE,ISOT))) THEN
                              NQC(IV,ISTATE,ISOT)= NDAT(IV,ISOT,ISTATE)
     1                                           - NRC(IV,ISTATE,ISOT)
                              WRITE(6,702) NDAT(IV,ISOT,ISTATE),
     1   SLABL(ISTATE),ISOT,IV,NRC(IV,ISTATE,ISOT),NQC(IV,ISTATE,ISOT)
                              ENDIF
                          IF(NQC(IV,ISTATE,ISOT).GT.0) THEN
c** Use parameter counter NQPAR to store the # free parameters preceeding
c  the first free doubling band constants for level IV of state ISTATE of 
c  isotopomer ISOT, &  Accumulate total No. free parameters:  NPARM
                              NQPAR(IV,ISTATE,ISOT)= NPARM
                              NPARM= NPARM+ NQC(IV,ISTATE,ISOT)
                              IF(IOMEG(ISTATE).GT.0) THEN
                                  WRITE(7,763) (CCDC(M+MQ0),IV,ISOT,
     1                                       M= 1,NQC(IV,ISTATE,ISOT))
                                  WRITE(6,684) IV,ISOT,
     1                                             NQC(IV,ISTATE,ISOT)
                                  ENDIF
                              IF(IOMEG(ISTATE).LT.0) THEN
                                  WRITE(7,764) (CCDC(M),IV,ISOT,
     1                                       M= 1,NQC(IV,ISTATE,ISOT))
                                  WRITE(6,685) IV,ISOT,
     1                                             NQC(IV,ISTATE,ISOT)
                                  ENDIF
                              ENDIF
                          ENDDO
                      ENDDO
                  NQPAR(VMAX(ISTATE)+1,ISTATE,NISTP)= NPARM
                  ENDIF
c
              IF(NDELD(ISTATE).GE.0) THEN
c** If fit to doubling parameters using use Dunham-type expansions ...
c=======================================================================
c** Here NQPAR stores # free parameters preceeding 1'st free Dunham-type
c  doubling expansion parameter for this state
                  NQPAR(0,ISTATE,1)= NPARM
                  DO  M= 1,NLDMX(ISTATE)
                      MQM= MQ0+ M
                      IF(LDMAX(MQM,ISTATE).GE.0) THEN
                          NPARM= NPARM+ LDMAX(MQM,ISTATE)+1
                          IF(IOMEG(ISTATE).GT.0) THEN 
                              WRITE(6,706) SLABL(ISTATE),CCDC(MQM),
     1                                               LDMAX(MQM,ISTATE)
                              WRITE(7,779)(L,MQM,L= 0,LDMAX(MQM,ISTATE))
                            ELSE
                              WRITE(6,707) SLABL(ISTATE),CCDC(M),
     1                                                 LDMAX(M,ISTATE)
                              WRITE(7,780) (L,M,L= 0,LDMAX(M,ISTATE))
                            ENDIF
                          ENDIF
                      ENDDO
                  ENDIF
              IF(IOMEG(ISTATE).GT.0) THEN
                  IF(efREF(ISTATE).NE.0) WRITE(6,716) SLABL(ISTATE),
     1                                                   efREF(ISTATE)
                  IF(efREF(ISTATE).EQ.0) WRITE(6,718) SLABL(ISTATE)
                  ENDIF
              ENDIF
c=======================================================================
c** If fitting to BOB correction expansion coefficients for this state
c=======================================================================
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).LE.0)) THEN
              CATOM= NAME(1)
              ATOM2= 2
              IF(AN(1).EQ.AN(2)) ATOM2= 1
              DO ATOM= 1,ATOM2
                  LAMIN= 0
                  IF((ISTATE.EQ.1).AND.(BOB00.LE.0)) LAMIN= 1
                  DO   M= 0, BOBORD(ISTATE)
                      IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                          NPARM= NPARM+ LAMAX(ATOM,M,ISTATE)+ 1- LAMIN
                          WRITE(6,708) SLABL(ISTATE),CCDC(M),CATOM,
     1                                            LAMAX(ATOM,M,ISTATE)
                          IF((M.EQ.0).AND.(ISTATE.EQ.1).AND.
     1                                      (BOB00.LE.0)) WRITE(6,709)
                          WRITE(7,781) 
     1                       (CATOM,L,M,L= LAMIN,LAMAX(ATOM,M,ISTATE))
                          ENDIF
                      LAMIN= 0
                      ENDDO
                  CATOM= NAME(2)
                  ENDDO
              ENDIF
   90     IF(ISTATE.LT.NSTATEMX) IPSTATE(ISTATE+1)= NPARM
c
c** Add fluorescence series origin levels to parameter count.
c============================================================
      IF(NFSTOT.GT.0) THEN
          WRITE(7,783) (VP(FSBAND(I)),VPP(FSBAND(I)),
     1             EFP(IFIRST(FSBAND(I))),ISTP(FSBAND(I)),I= 1,NFSTOT)
          M= min(20,NFSTOT)
          WRITE(6,710) NFSTOT,(VP(FSBAND(I)),VPP(FSBAND(I)),
     1                  EFP(IFIRST(FSBAND(I))),ISTP(FSBAND(I)),I= 1,M)
          IF(NFSTOT.GT.M) WRITE(6,714) NFSTOT-20
          IF((NPARM+NFSTOT).GT.NPARMX) THEN
              WRITE(6,711) NPARM,NFSTOT,NPARMX
              STOP
              ENDIF
          DO  I= 1,NFSTOT
              PV(NPARM+I)= ORIGIN(I)
              ENDDO
c         write(6,715) (origin(i),I= 1,nfstot)
c 715 format('   with origins:',4F12.3)
          NPARM= NPARM+ NFSTOT
          ENDIF
c** Rewind channel-7 and read parameter names for final printout.
      REWIND(7)
      DO  I= 1,NPARM
          READ(7,785) NAMEPARM(I)
          IFXP(I)= 0
          ENDDO
      REWIND(7)
      IF(NPARM.GT.NPARMX) THEN
c** If need more parameters than dimensioning allows ... STOP
          WRITE(6,712) NPARM,NPARMX
          STOP
          ENDIF
c** Call NLLSSRR to do actual fit 
      JROUND= IROUND
      IF((IROUND.NE.0).AND.(NFSTOT+NTVALL(0).GT.0)) JROUND= 0
c***********************************************************************
c***********************************************************************
      CALL NLLSSRR(COUNTOT,NPARM,NPARMX,CYCMAX,JROUND,ROBUST,LPRINT,
     1              IFXP,FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
      IF(JROUND.NE.IROUND) THEN
c** Perform group rounding of fitted term values and/or fluorescence 
c   series origins. 
          DO  I= 1, NPARM
              PUSAV(I)= PU(I)
              PSSAV(I)= PS(I)
              ENDDO
          JROUND= IABS(IROUND)+ 1
c** Round all term values for each state in a single step 
          IF(NTVALL(0).GT.0) THEN
              DO  ISTATE= 1, NSTATES
                  IF(NTVALL(ISTATE).GT.0) CALL GPROUND(JROUND,NPARM,
     1        NPARMX,IPSTATE(ISTATE)+1,IPSTATE(ISTATE)+NTVALL(ISTATE),
     2                                              LPRINT,IFXP,PV,PU)
                  ENDDO
              ENDIF
c** Round all fluorescence series origins in a single step 
          IF(NFSTOT.GT.0) THEN
              I= NPARM- NFSTOT+ 1
              CALL GPROUND(JROUND,NPARM,NPARMX,I,NPARM,LPRINT,IFXP,
     1                                                          PV,PU)
              ENDIF
c ... and then call NLLSSRR again to sequentially round remaining parm.
          CALL NLLSSRR(COUNTOT,NPARM,NPARMX,CYCMAX,IROUND,ROBUST,LPRINT,
     1              IFXP,FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c ... and finally, reset all parameter uncertainties at original values
          DO  I= 1, NPARM
              PU(I)= PUSAV(I)
              PS(I)= PSSAV(I)
              ENDDO
          DSE= DSE*DSQRT(DFLOAT(COUNTOT- NPARM+ NFSTOT+ NTVALL(0))/
     1                                         DFLOAT(COUNTOT- NPARM))
          ENDIF
c***********************************************************************
c***********************************************************************
c** Now ... print results of the fit & final parameter values ...
      IF(IROUND.NE.0) WRITE(6,720) NPARM,COUNTOT,DSE
      IF(IROUND.EQ.0) WRITE(6,722) NPARM,COUNTOT,TSTPS,DSE,TSTPU
c
      DO 110 ISTATE= 1,NSTATES
          I2= IPSTATE(ISTATE)
          I1= I2+ 1
          IF(IFXGv(ISTATE).LE.0) THEN
c** Write fitted parameters for Gv Expansion
              IF(NDEGv(ISTATE).EQ.-2) THEN
                  WRITE(6,675) SLABL(ISTATE),NTVALL(ISTATE)
                  I2= IPSTATE(ISTATE)+ NTVALL(ISTATE)
                  ENDIF
              IF(NDEGv(ISTATE).EQ.-1) THEN
                  WRITE(6,724) SLABL(ISTATE),VMIN(ISTATE),VMAX(ISTATE)
                  I2= I2+ NRBC(ISTATE) 
                  ENDIF
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).EQ.2)) THEN
                  IF(LMAX(0,ISTATE).GT.0) I2= I2+ LMAX(0,ISTATE)
                  IF(ISTATE.GT.1) I2= I2+1
                  IF(NDEGv(ISTATE).EQ.0) WRITE(6,726) SLABL(ISTATE),
     1                                                         CCDC(0)
                  ENDIF
              IF(NDEGv(ISTATE).GT.0) THEN
                  IF(NDEGv(ISTATE).EQ.2) THEN
                      WRITE(6,730) SLABL(ISTATE),CCDC(0),LMAX(0,ISTATE),
     1      VS(ISTATE),VS(ISTATE),DVS(ISTATE),NP0(ISTATE),NQ0(ISTATE),
     2                CTYPE(ITYPE(ISTATE)),IP0(ISTATE)+1,IQ0(ISTATE)+1
                      IF((IFXVS(ISTATE).LE.0).OR.(IFXDVS(ISTATE).LE.0))
     1                                                            THEN
                          IF(IFXVS(ISTATE).LE.0) I2= I2+1
                          IF(IFXDVS(ISTATE).LE.0) I2= I2+1
                          WRITE(6,731)
                          ENDIF
                      ENDIF
                  IF(IFXD(ISTATE).LE.0) I2= I2+1
                  IF(IFXVD(ISTATE).LE.0) I2= I2+ 1
                  I2= I2+ NP0(ISTATE)+ NQ0(ISTATE)
                  IF(NDEGv(ISTATE).EQ.1) WRITE(6,728) SLABL(ISTATE),
     1           CCDC(0),NP0(ISTATE),NQ0(ISTATE),CTYPE(ITYPE(ISTATE)),
     2                                     IP0(ISTATE)+1,IQ0(ISTATE)+1
c???              I= VD(ISTATE)
c???              VMAX(ISTATE)= MIN(I,NVIBMX)
c???          CALL NDEDGB(ISTATE,NISTP,NEWGv,NEWBv,RSQMU,VMAX(ISTATE))
                  ENDIF
              IF(I2.GE.I1) THEN
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
c** For NON term value fits, write vibrational energies to channel-7
                  IF(NDEGv(ISTATE).GE.-1) THEN
                      DO  ISOT= 1,NISTP
                          WRITE(7,791)NAME(1),MN(1,ISOT),NAME(2),
     1MN(2,ISOT),SLABL(ISTATE),(I,ZK(0,I,ISTATE,ISOT),I= 0,VMAX(ISTATE))
                          ENDDO
                      ENDIF
                  ENDIF
              IF(NDEGv(ISTATE).EQ.-2) GO TO 108
              IF(NDEGv(ISTATE).EQ.-1) GO TO 100
              ENDIF
          IF(IFXBv(ISTATE).LE.0) THEN
c-----------------------------------------------------------------------
c** Write fitted parameters for Bv Expansion
c-----------------------------------------------------------------------
              IF((NDEBv(ISTATE).EQ.-1).AND.(NDEGv(ISTATE).GE.0).AND.
     1                                       (NRBC(ISTATE).GT.0)) THEN
                  I1= I2+1
                  I2= I2+ NRBC(ISTATE)
                  WRITE(6,734) SLABL(ISTATE)
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  ENDIF
              IF(NDEBv(ISTATE).EQ.0) WRITE(6,726) SLABL(ISTATE),CCDC(1)
              IF(NDEBv(ISTATE).EQ.1) WRITE(6,728) SLABL(ISTATE),
     1            CCDC(1),NP1(ISTATE),NQ1(ISTATE),CTYPE(ITYPB(ISTATE))
              IF(NDEBv(ISTATE).GE.2) WRITE(6,730) SLABL(ISTATE),
     1       CCDC(1),LMAX(1,ISTATE),VS(ISTATE),VS(ISTATE),DVS(ISTATE),
     2                    NP1(ISTATE),NQ1(ISTATE),CTYPE(ITYPB(ISTATE))
              IF(((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).GE.2))
     1                                .AND.(LMAX(1,ISTATE).GE.0)) THEN
                  I1= I2+1
                  I2= I2+ LMAX(1,ISTATE)+1
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  ENDIF
              IF((NDEBv(ISTATE).GE.1).AND.(NP1(ISTATE).GT.0)) THEN
                  I1= I2+1
                  I2= I2+ NP1(ISTATE)
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  ENDIF
              IF((NDEBv(ISTATE).GE.1).AND.(NQ1(ISTATE).GT.0)) THEN
                  I1= I2+1
                  I2= I2+ NQ1(ISTATE)
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  ENDIF
c** Write inertial rotational constants to channel-7
              DO  ISOT= 1,NISTP
                  WRITE(7,793)NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1         SLABL(ISTATE),(I,ZK(1,I,ISTATE,ISOT),I= 0,VMAX(ISTATE))
                  ENDDO
              ENDIF
          IF(IFXCDC(ISTATE).LE.0) THEN
c-----------------------------------------------------------------------
c** Write fitted parameters for CDC expansion(s)
c-----------------------------------------------------------------------
              IF((NDECDC(ISTATE).EQ.-1).AND.(NDEBv(ISTATE).GE.0).AND.
     1                                       (NRBC(ISTATE).GT.0)) THEN
                  I1= I2+1
                  I2= I2+ NRBC(ISTATE)
                  WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  ENDIF
              IF(NDECDC(ISTATE).EQ.0) THEN
                  DO  M= 2,NCDC(ISTATE)+1
                      I1= I2+ 1
                      IF(LMAX(M,ISTATE).GE.0) I2= I2+ LMAX(M,ISTATE)+1
                      IF(I2.GE.I1) THEN
                          WRITE(6,736) SLABL(ISTATE),M-1
                          WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),
     1                                                      I=  I1,I2)
                          ENDIF
                      ENDDO
                  ENDIF
              DO  M= 1,NCDC(ISTATE)
                  DO  ISOT= 1,NISTP
                      WRITE(7,795) CCDC(M+1),NAME(1),MN(1,ISOT),NAME(2),
     1                     MN(2,ISOT),SLABL(ISTATE),CCDC(M+1),
     2                     (I,ZK(M+1,I,ISTATE,ISOT),I= 0,VMAX(ISTATE))
                      ENDDO
                  ENDDO
              DO  ISOT= 1,NISTP
                  DO  I= 0,VMAX(ISTATE)
                      WRITE(7,797)  I,(ZK(M,I,ISTATE,ISOT),
     1                                            M= 0,NCDC(ISTATE)+1)
                      ENDDO
                  ENDDO
              ENDIF 
c-----------------------------------------------------------------------
c** Write fitted parameters for Lambda/Gamma-doubling expansions.
c??????????????????????????????????????????????????????????????????????
c If IOMEG < 0 has been reset to 0 (only sigma case considered so far!!)
c , thus need to test the value of MULTPLT when considering gamma-doubling
c-----------------------------------------------------------------------
c          IF(((IOMEG(ISTATE).GT.0).OR.(MULTPLT(ITSTATE).EQ.2)).AND.
c     1                                     (IFXLD(ISTATE).LE.0)) THEN
c          IF((MULTPLT(ISTATE).EQ.2).AND.(IFXLD(ISTATE).LE.0)) THEN
c          IF(((IOMEG(ISTATE).GT.0).OR.(MULTPLT(ISTATE).EQ.2)).AND.
c          IF(((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).LT.0)).AND.
c     1            (IFXLD(ISTATE).LE.0).AND.(NDELD(ISTATE).GE.0)) THEN
c??????????????????????????????????????????????????????????????????????
  100     IF((IOMEG(ISTATE).NE.0).AND.(IFXLD(ISTATE).LE.0)
     1                                 .AND.(NLDMX(ISTATE).GT.0)) THEN
              MQ0= MAX0(0,IOMEG(ISTATE)- 1)
              IF(NDELD(ISTATE).EQ.-1) THEN
c ... if doubling parameters represented by isotopic band constants
                  WRITE(6,704) SLABL(ISTATE)
                  I1= I2+1
                  DO  ISOT= 1, NISTP
                      DO IV= VMIN(ISTATE), VMAX(ISTATE)
                          I2= I2+ NQC(IV,ISTATE,ISOT)
                          ENDDO
                      ENDDO
                  IF(I2.GE.I1) WRITE(6,732) (NAMEPARM(I),PV(I),
     1                                           PU(I),PS(I),I= I1,I2)
                  ENDIF
c
              IF(NDELD(ISTATE).GE.0) THEN
c ... if doubling constants represented by Dunham-type expansions ...
                  I1= I2+ 1
                  IF(IOMEG(ISTATE).GT.0) THEN
                      WRITE(6,742) SLABL(ISTATE)
                      IF(efREF(ISTATE).NE.0) WRITE(6,716) SLABL(ISTATE),
     1                                                   efREF(ISTATE)
                      IF(efREF(ISTATE).EQ.0) WRITE(6,718) SLABL(ISTATE)
                      ENDIF
                  IF(IOMEG(ISTATE).LT.0) WRITE(6,743) SLABL(ISTATE)
                  DO  M= 1,NLDMX(ISTATE)
                      MQM= MQ0+ M
                      IF(LDMAX(MQM,ISTATE).GE.0) I2= I2 + 
     1                                             LDMAX(MQM,ISTATE)+1
                      ENDDO
                  IF(I2.GE.I1)
     1            WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
                  DO  M= 1,NLDMX(ISTATE)
                      MQM= MQ0+ M
                      DO  ISOT= 1, NISTP
                          WRITE(7,796) CCDC(MQM),ISOT,NAME(1),
     1          MN(1,ISOT),NAME(2),MN(2,ISOT),SLABL(ISTATE),CCDC(MQM),
     2                     (I,ZQ(MQM,I,ISTATE,ISOT),I= 0,VMAX(ISTATE))
                          ENDDO
                      ENDDO
                  ENDIF
              ENDIF
c** Write fitted parameters for B-O-Breakdown expansion(s)
c=======================================================================
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).LE.0)) THEN
              I1= I2+ 1
              WRITE(6,744) SLABL(ISTATE)
              DO  ATOM= 1, ATOM2
                  DO  M= 0, BOBORD(ISTATE)
                      IF(LAMAX(ATOM,M,ISTATE).GE.0) THEN
                         I2= I2+ LAMAX(ATOM,M,ISTATE)+ 1
c** NOTE: for ISTATE=1  must specify whether include vib. BOB {0.0} term
                         IF((M.EQ.0).AND.(ISTATE.EQ.1).AND.(BOB00.LE.0))
     1                                                        I2= I2-1
                         ENDIF
                      ENDDO
                  ENDDO
              IF(I2.GE.I1)
     1           WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
              ENDIF
  108     WRITE(6,*)
  110     CONTINUE
c** Print the correlation matrix to Channel-10 ... ignoring correlation
c  to any fluorescence series origins
      IF(MKPRED.LE.0) THEN
          WRITE(10,748) CM(1,1),I1,CM(2,1),CM(2,2)
          IF(I1.GE.3) THEN 
              DO  J= 3, I1
                  WRITE(10,750) (CM(j,i), i= 1,J)
                  ENDDO
              ENDIF
          ENDIF
c** Now print fluorescence series origins (if any)
      I1= I2+ 1
      I2= NPARM
      IF(I2.GE.I1) THEN
          WRITE(6,752) (I2-I1+1)
          WRITE(6,732) (NAMEPARM(I),PV(I),PU(I),PS(I),I= I1,I2)
          ENDIF
c
c** Calculate Y00(semiclass) & zero point energy and their uncertainties
c  and for multi-isotopomer Dunham case, generate and print rounded YLM 
c  parameters (with their proper uncertainties) for minority isotopomers
      CALL PPISOT(NISTP,AN,MN,PV,PU,PS,CM,ZMASS,RSQMUP,RMUP,NAME,
     1                                                 SLABL,NAMEPARM)
      DO  ISTATE= 1, NSTATES
c** Calculate and output values of derived parameters (for NDE case)
          IF((IFXGv(ISTATE).LE.0).AND.(NDEGv(ISTATE).EQ.1)) THEN
              I=0
              CALL NDEDGB(ISTATE,NISTP,NEWGv,NEWBv,RSQMU,I)
              CALL NDEDUN(ISTATE,NISTP,ZMASS,RSQMU,PU,CM)
              ENDIF
          ENDDO
c=======================================================================
  675 FORMAT(/' State ',A3,i6,'  fitted term values   Tv{state:v,J,p;iso
     1t}'/1x,6('======'))
  676 FORMAT(/' Use a Band-Constant fit for State ',A3,2x,a3,"'s":
     1  2x,A3,"'s": '  and  ',A3,"'s")
  677 FORMAT(1x,6('******'))
  678 FORMAT('* WARNING: find',i2,' data for State ',A3,'  ISOT-',i1,
     1  '  v=',I3,'.  Set (FITGV,NRC)=(',I1,',',i1,')')
  680 FORMAT(' Base connected-level set #',I2, '  at  State ',A3,
     1  '   v=',I3,'   of isotopomer-',i1)
  681 FORMAT(5x,'For  v=',i3,'  of  ISOT=',i2,'  fit to',i3,' CDC band c
     1onstants')
  682 FORMAT(5x,'For  v=',i3,'  of  ISOT=',i2,'  fit to the energy and',
     1  i3,' rotational band constants')
  683 FORMAT(5x,'For  v=',i3,'  of  ISOT=',i2,'  fit to',i3,' rotational
     1 band constants')
  684 FORMAT(5x,'For  v=',i3,'  of  ISOT=',i2,'  fit to',i2,' Lambda dou
     1bling band constants')
  685 FORMAT(5x,'For  v=',i3,'  of  ISOT=',i2,'  fit to',i2,' Gamma doub
     1ling band constants')
  686 FORMAT(/' Absolute zero of energy is fixed at  G(v"=-1/2)  of  Sta
     1te ',A2/1x,12('**'))
  688 FORMAT(/' Fit for State ',A3,1x,A3,"'s  uses Dunham expansion of o
     1rder",i3/1x,7('***'))
  689 FORMAT(/" Fit for State ",A3,1x,A3,"'s  uses Tellinghuisen-type MS
     1X mixed representation:"/1x,8('**'),' order',i3,' Dunham for  v .l
     2e. VS=',F11.6,'  &  NDE for  v > VS':/9x, 'with switching function
     3:   Sw = 1/[1 + exp{(v -',F11.6,')/',F9.6,'}]')
  690 FORMAT(9x,'and treat  VS  as a free parapeter in the fit')
  691 FORMAT(9x,'and treat  DVS  as a free parapeter in the fit')
  692 FORMAT(/" Fit 3tate ",A3,"  Gv's to NDE or MXS function while fixi
     1ng   DLIMIT=",F12.5/1x,43('*')) 
  694 FORMAT(/" Fit to State ",A3, "  DLIMIT  and to NDE or MXS function
     1 for  Gv"/)
  698 FORMAT(" State ",A3,1x,A3,"'s  initially defined by  (NP=",I2,
     1 '/NQ=',I2,') ',A12,'NDE  in  (vD-v)'/1x,6('--'),'  with leading n
     2umerator and denominator powers',I3,' &',I3:/15X,
     3  'where for Isotopomer-1   vD=',F13.8)
  700 FORMAT(/' Dunham Fit for State ',A3,'  uses  CDC(',i1,')=',A3,'  e
     1xpansion of order',i3:/' =======================',8x,'CDC(',i1,
     2  ')=',A3,'  expansion of order',I3:/(32x,'CDC(',i1,')=',A3,
     3  '  expansion of order',I3:))
  702 FORMAT('* WARNING: find',i2,' data for State ',A3,'  ISOT-',i1,
     1  '  v=',I3,'  with  NRC=',i2,', so set  NQC=',I2)
  704 FORMAT(/' Band-Constant Doubling parameters for Electronic State '
     1 ,A3)
  706 FORMAT(' Fit state ',A3,' "',A3, '"-type Lambda doubling constants
     1 to order',I3,'  Dunham expansion')
  707 FORMAT(' Fit state ',A3,' "',A3, '"-type Gamma  doubling constants
     1 to order',I3,'  Dunham expansion')
  708 FORMAT(' Fit State ',A3,'  "',A3,'" BOB corrections for ',A2,
     1  '  to Dunham expansion of order:',I3)
  709 FORMAT(21x,'while IGNORING the constant  delta(0,0)  term')
  710 FORMAT(/' Fit to the origins of',i5,' fluorescence series with ini
     1tial-state labels'/21x,"(where  p= 'parity'  &  IS= 'isotope'):"/
     2  5(2x,13('-'),1x)/5("   v'  J'  p IS ")/5(2x,13('-'),1x)/
     3  5(i5,i4,SP,I3,SS,I3,1x:))
  711 FORMAT(/' *** ERROR *** Dimension allocated for number of paramete
     1rs exceeded:'/15x,'(NPARM=',i4,') + (NFSTOT=',i4,') > (NPARMX=',
     2  i4,')')
  712 FORMAT(/' **** Case=2  G-B-C-FIT Option FAILS because # of free pa
     1rameters   NPARM=',I5,'  exceeds array dimension limit   NPARMX=',
     2  I5)
  714 FORMAT('  ........ and',i6,' others ...........')
  716 FORMAT(' => State ',A3,' parity=',SP,i3,'  sublevels treated as un
     1perturbed reference energy.')
  718 FORMAT(' => State ',A3,' parity sublevel  midpoint  treated as unp
     1erturbed reference energy')
  720 FORMAT(/' After Sequential Rounding & Refitting,  fit of',i6,
     1 ' parameters to',i6,' data'/1x,37('*'),'   yields    DSE=',
     2  G11.4/)
  722 FORMAT(/' Fit',i5,' parameters to',i6,' data:   Test(PS)=',1PD8.1,
     1  '   DSE=',0PG11.4/1x,11('***'),'    Test(PU)=',1PD8.1/)
  724 FORMAT(' State ',A3,'  Pure Band-Constant representation for level
     1s   v=',i2,' -',i3)
  726 FORMAT(' State ',A3,'  Dunham expansion',A3,' parameters:')
  728 FORMAT(' State ',A3,'  NDE ',A3,'  function is a  (NP=',I2,'/NQ='
     1 ,I2,') ',A12,'NDE  in  (vD-v)':/17x,'with leading numerator and d
     2enominator powers',i3,' &',i3)
  730 FORMAT(' State ',A3,'  Tellinghuisen-type MSX mixed representation
     1',A3,'  parameters based on:'/1x,4('=='),I5, "'th order Dunham for
     2  v .le. VS=",F11.6,'   &   NDE for  v > VS':/9x,'with switching f
     3unction:   Sw = 1/[1 + exp{(v -',F11.6,')/',F9.6,'}]'/9x,'where ND
     4E function is a  (NP=',I2,'/NQ=',I2,') ',A12,'NDE  in  (vD-v)':/
     5 17x,'with leading numerator and denominator powers',i3,' &',i3)
  731 FORMAT(6x,'and fit optimizes switching function parameters  VS and
     1/or DVS')
  732 FORMAT(1x,a20,'=',1PD20.12,' (+/-',D8.1,')    Sensitivity=',D8.1)
  734 FORMAT(' State ',A3,'  Band-Constant rotational constants:')
  736 FORMAT(' State ',A3,'  Dunham expansion CDC(',i1,') parameters:')
  742 FORMAT(' State ',A3,'  Lambda-doubling Dunham-type expansion coeff
     1icients:')
  743 FORMAT(' State ',A3,'  Gamma-doubling Dunham-type expansion coeffi
     1cients:')
  744 FORMAT(' State ',A3,'  Born-Oppenheimer breakdown parameters:')
  748 FORMAT(/f10.6,13x,'Correlation Matrix linking the first ',i4,
     1  ' parameters'/2F10.6,3x,26('=='))
  750 FORMAT(8F10.6:/(10x,7F10.6:))
  752 FORMAT(" Energy origins  FS(v', j', p'; isotopomer)  of the",
     1  i5,' fluorescence series')
  761 FORMAT(7x,A3,'(v=',i3,';',i2,')')
  763 FORMAT(4x,'q[',A3,'(v=',i3,';',i2,')]')
  764 FORMAT(4x,'g[',A3,'(v=',i3,';',i2,')]')
  765 FORMAT(10x,'T(v= -1/2',A2,')')
  766 FORMAT(11x,'YLM(',i2,',',i1,')')
  768 FORMAT(13x,'VS(',A2,')')
  769 FORMAT(12x,'dVS(',A2,')')
  770 FORMAT(10x,'DLimit(',A2,')')
  771 FORMAT(14x,'vD(',A2,')')
  773 FORMAT(14x,'P',i1,'(',i2,')')
  775 FORMAT(14x,'Q',i1,'(',i2,')')
  779 FORMAT(10x,'qLM(',i2,',',i1,')')
  780 FORMAT(10x,'gLM(',i2,',',i1,')')
  781 FORMAT(6x,'delta(',A2,';',i2,',',i1,')')
  783 FORMAT('  FS(',SS,i3,',',i3,',',SP,i3,';',SS,i2,')') 
  785 FORMAT(A20)
  791 FORMAT(' Vibrational energies for ',A2,'(',i3,')-',A2,'(',I3,
     1  ') in State ',A3,':  {v, Gv}'/(4(i5,f15.6)))
  793 FORMAT(' Inertial rotational constants Bv for ',A2,'(',i3,')-',
     1  A2,'(',I3,') in State ',A3,':  {v, Bv}'/(4(i5,f15.8)))              
  795 FORMAT(' Centrifugal distortion constants ',A3,' for ',A2,'(',
     1  i3,')-',A2,'(',I3,') in State ',A3,':  {v,',A3,'}'/
     2  (4(i5,1PD15.7)))
  796 FORMAT(' Lambda Doubling constants q[',A3,'(v;',i2,')] for ',
     1  A2,'(',i3,')-',A2,'(',I3,') in State ',A3,': {v,',A3,'}'/
     2  (4(i5,1PD15.7)))
  797 FORMAT(I4,f12.4,f14.10,6(1PD15.7:))
c=======================================================================
c=======================================================================
c** Now ... calculate band-by-band DSE values for output summary
      CALL DIFFSTATS(NSTATES,ROBUST,MKPRED)
  200 STOP
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,DGNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy DGNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2012 mass table [Wang, Audi, Wapstra, Kondev, MacCormick, Xu
c  & Pfeiffer, Chin.Phys.C 36, 1603-2014 (2012)] ,the proton, deuteron,
c  and triton masses are taken from the 2010 fundamental constants table 
c  [Mohr, Taylor, & Newell, Rev. Mod. Phys. 84, 1587-1591 (2012)] and other
c  quantities from Tables 6.2 and 6.3 of "Quantities, Units and Symbols in
c  Physical Chemistry", by Mills et al.(Blackwell,2'nd Edition, Oxford,1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set DGNS=-1 and ABUND=-1.
c** For Atomic number IAN=0 and isotope mass numbers IMN=1-3,  return the
c    masses of the proton, deuteron, and triton, p,d & t, respectively
c Masses and properties of selected Halo nuclei an unstable nuclei included
c                 COPYRIGHT 2005-2015  :  last  updated 10 January 2016
c** By R.J. Le Roy, with assistance from 
c                 G.T. Kraemer, J.Y. Seto and K.V. Slaughter.
c***********************************************************************
      REAL*8 zm(0:123,0:15),mass,ab(0:123,15),abund
      INTEGER i,ian,imn,gel(0:123),nmn(0:123),mn(0:123,15),
     1                                        gns(0:123,15),DGNS,gelgs
      CHARACTER*2 NAME,AT(0:123)
cc
      DATA  at(0),gel(0),nmn(0),(mn(0,i),i=1,3)/' p',1,3,1,2,3/
      DATA  (zm(0,i),i=0,3)/1.008d0,1.007276466812d0,2.013553212712d0,
     2                 3.0155007134d0/
      DATA  (gns(0,i),i=1,3)/2,3,2/
      DATA  (ab(0,i),i=1,3)/0.d0, 0.d0, 0.d0/
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503223d0, 2.01410177812d0,
     1                 3.0160492779d0/
      DATA  (gns(1,i),i=1,3)/2,3,2/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,4)/'He',1,4,3,4,6,8/
      DATA  (zm(2,i),i=0,4)/4.002602d0, 3.0160293201d0, 4.00260325413d0,
     1                                        6.0188891d0, 8.033922d0/
      DATA  (gns(2,i),i=1,4)/2,1,1,1/
      DATA  (ab(2,i),i=1,4)/0.000137d0,99.999863d0, 2*0.d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,6)/'Li',2,6,6,7,8,9,11,12/
      DATA  (zm(3,i),i=0,6)/6.941d0, 6.0151228874d0, 7.016003437d0,
     1     8.02248736d0,9.0267895d0,11.043798d0,12.05378d0/
      DATA  (gns(3,i),i=1,6)/3,4,5,4,4,1/
      DATA  (ab(3,i),i=1,6)/7.5d0, 92.5d0, 4*0.d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,8)/'Be',1,8,7,9,10,11,12,
     1                                                       14,15,16/
      DATA  (zm(4,i),i=0,8)/9.012182d0, 7.01692983d0, 9.01218307d0,
     1 10.0135338d0, 11.021658d0, 12.026921d0, 14.04289d0, 15.05346d0,
     2 16.06192d0/
      DATA  (gns(4,i),i=1,8)/4,4,3,2,1,1,2,1/
      DATA  (ab(4,i),i=1,8)/0.d0, 100.d0, 6*0.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,10)/' B',2,10,8,10,11,12,
     1                                              13,14,15,17,18,19/
      DATA (zm(5,i),i=0,10)/10.811d0, 8.0246072d0, 10.0129369d0, 
     1          11.0093054d0, 12.0143521d0, 13.0177802d0, 14.025404d0,
     2          15.031103d0, 17.04699d0, 18.05617d0,19.06373d0/
      DATA  (gns(5,i),i=1,10)/5,7,4,3,4,5,4,4,1,4/
      DATA  (ab(5,i),i=1,10)/0.d0, 19.9d0,80.1d0, 7*0.d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,14)/' C',1,14,9,10,11,12,13,
     1               14,15,16,17,18,19,20,21,22/
      DATA (zm(6,i),i=0,14)/12.011d0, 9.0310367d0, 10.0168532d0,
     1          11.0114336d0, 12.d0, 13.00335483507d0, 14.003241989d0, 
     1  15.0105993d0, 16.014701d0, 17.022586d0, 18.02676d0, 19.03481d0,
     2  20.04032d0, 21.04934d0, 22.05720d0/
      DATA  (gns(6,i),i=1,14)/4,1,4,1,2,1,2,1,4,1,2,1,2,1/
      DATA  (ab(6,i),i=1,14)/3*0.d0, 98.90d0,1.10d0, 9*0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.00307400443d0,15.0001088989d0/
      DATA (gns(7,i),i=1,2)/3,2/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461957d0, 16.9991317565d0,
     1                      17.9991596129d0/
      DATA (gns(8,i),i=1,3)/1,6,1/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.9984031627d0/
      DATA (gns(9,i),i=1,1)/2/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,4)/'Ne',1,4,17,20,21,22/
      DATA (zm(10,i),i=0,4)/20.1797d0, 17.017672d0, 19.9924401762d0, 
     1                                   20.99384669d0,21.991385115d0/
      DATA (gns(10,i),i=1,4)/2,1,4,1/
      DATA (ab(10,i),i=1,4)/0.d0, 90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692820d0/
      DATA (gns(11,i),i=1,1)/4/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041698d0, 24.98583698d0,
     1                       25.98259297d0/
      DATA (gns(12,i),i=1,3)/1,6,1/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153853d0/
      DATA (gns(13,i),i=1,1)/6/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265346d0, 28.9764946649d0,
     1                       29.973770136d0/
      DATA (gns(14,i),i=1,3)/1,2,1/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,2)/' P',4,2,26,31/
      DATA (zm(15,i),i=0,2)/30.973762d0, 26.01178d0, 30.9737619984d0/
      DATA (gns(15,i),i=1,2)/15,2/
      DATA (ab(15,i),i=1,2)/0.d0, 100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,5)/' S',5,5,27,32,33,
     1                                                          34,36/
      DATA (zm(16,i),i=0,5)/32.066d0, 27.01883d0, 31.9720711744d0,
     1                   32.9714589098d0,33.96786700d0, 35.96708071d0/
      DATA (gns(16,i),i=1,5)/6,1,4,1,1/
      DATA (ab(16,i),i=1,5)/0.d0, 95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590260d0/
      DATA (gns(17,i),i=1,2)/4,4/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545105d0, 37.96273211d0,
     1                       39.9623831237d0/
      DATA (gns(18,i),i=1,3)/1,1,1/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.963706486d0, 39.96399817d0,
     1                       40.961825258d0/
      DATA (gns(19,i),i=1,3)/4,9,4/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.962590864d0, 41.95861783d0,
     1         42.95876644d0, 43.9554816d0, 45.9536890d0, 47.95252277d0/
      DATA (gns(20,i),i=1,6)/1,1,8,1,1,1/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559083d0/
      DATA (gns(21,i),i=1,1)/8/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526277d0, 46.9517588d0,
     1         47.9479420d0, 48.9478657d0, 49.9447869d0/
      DATA (gns(22,i),i=1,5)/1,6,1,8,1/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471560d0, 50.9439570d0/
      DATA (gns(23,i),i=1,2)/13,8/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460418d0, 51.9405062d0,
     1                       52.9406481d0, 53.9388792d0/
      DATA (gns(24,i),i=1,4)/1,1,4,1/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.938049d0/
      DATA (gns(25,i),i=1,1)/6/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396090d0, 55.9349363d0,
     1                       56.9353928d0, 57.9332744d0/
      DATA (gns(26,i),i=1,4)/1,1,2,1/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331943d0/
      DATA (gns(27,i),i=1,1)/8/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353424d0, 59.9307859d0,
     1         60.9310556d0, 61.9283454d0, 63.9279668d0/
      DATA (gns(28,i),i=1,5)/1,1,4,1,1/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295977d0,64.9277897d0/
      DATA (gns(29,i),i=1,2)/4,4/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291420d0, 65.9260338d0,
     1         66.9271277d0, 67.9248446d0, 69.9253192d0/
      DATA (gns(30,i),i=1,5)/1,1,6,1,1/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255735d0, 70.9247026d0/
      DATA (gns(31,i),i=1,2)/4,4/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242488d0, 71.92207583d0,
     1         72.92345896d0, 73.921177762d0, 75.921402726d0/
      DATA (gns(32,i),i=1,5)/1,1,10,1,1/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215946d0/
      DATA (gns(33,i),i=1,1)/4/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.922475935d0, 75.919213704d0,
     1         76.91991415d0, 77.91730928d0, 79.9165218d0, 81.9166995d0/
      DATA (gns(34,i),i=1,6)/1,1,2,1,1,1/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183376d0, 80.9162897d0/
      DATA (gns(35,i),i=1,2)/4,4/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203649d0, 79.9163781d0,
     1     81.9134827d0, 82.9141272d0, 83.911497728d0, 85.910610627d0/
      DATA (gns(36,i),i=1,6)/1,1,1,10,1,1/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180532d0/
      DATA (gns(37,i),i=1,2)/6,4/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.9134191d0, 85.9092606d0,
     1                      86.9088775d0, 87.9056125d0/
      DATA (gns(38,i),i=1,4)/1,1,10,1/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058403d0/
      DATA (gns(39,i),i=1,1)/2/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9046977d0, 90.9056396d0,
     1                      91.9050347d0, 93.9063108d0, 95.9082714d0/
      DATA (gns(40,i),i=1,5)/1,6,1,1,1/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063730d0/
      DATA (gns(41,i),i=1,1)/10/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.9068080d0, 93.9050849d0,
     1        94.9058388d0, 95.9046761d0, 96.9060181d0, 97.9054048d0,
     2        99.9074718d0/
      DATA (gns(42,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907212d0/
      DATA (gns(43,i),i=1,1)/13/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.9075903d0, 97.905287d0,
     1     98.9059341d0, 99.9042143d0, 100.9055769d0, 101.9043441d0,
     2     103.9054275d0/
      DATA (gns(44,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.9054980d0/
      DATA (gns(45,i),i=1,1)/2/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.9056022d0, 103.9040305d0,
     1       104.9050796d0, 105.9034804d0, 107.9038916d0, 109.9051722d0/
      DATA (gns(46,i),i=1,6)/1,1,6,1,1,1/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.9050916d0, 108.9047553d0/
      DATA (gns(47,i),i=1,2)/2,2/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.9064599d0, 107.9041834d0, 
     1       109.9030066d0, 110.9041829d0, 111.9027629d0, 112.9044081d0,
     2       113.9033651d0, 115.90476315d0/
      DATA (gns(48,i),i=1,8)/1,1,1,2,1,2,1,1/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.9040618d0, 114.903878776d0/
      DATA  (gns(49,i),i=1,2)/10,10/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.9048239d0, 113.9027827d0,
     1    114.903344699d0, 115.90174280d0, 116.9029540d0, 117.9016066d0,
     2    118.9033112d0, 119.9022016d0, 121.9034438d0, 123.9052766d0/
      DATA (gns(50,i),i=1,10)/1,1,2,1,2,1,2,1,1,1/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.903812d0, 122.9042132d0/
      DATA (gns(51,i),i=1,2)/6,8/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904059d0, 121.9030435d0,
     1    122.9042698d0, 123.9028171d0, 124.9044299d0, 125.9033109d0,
     2    127.9044613d0, 129.906222749d0/
      DATA (gns(52,i),i=1,8)/1,1,2,1,2,1,1,1/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904472d0, 128.904984d0/
      DATA (gns(53,i),i=1,2)/6,8/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058920d0, 125.904298d0,
     1    127.9035310d0, 128.904780861d0,129.903509350d0,130.90508406d0,
     2    131.904155086d0, 133.9053947d0, 135.907214484d0/
      DATA (gns(54,i),i=1,9)/1,1,1,2,1,4,1,1,1/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451961d0/
      DATA (gns(55,i),i=1,1)/8/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063207d0, 131.9050611d0,
     1    133.90450818d0, 134.90568838d0, 135.90457573d0, 136.9058271d0,
     2    137.9052470d0/
      DATA (gns(56,i),i=1,7)/1,1,1,4,1,4,1/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907115d0, 138.9063563d0/
      DATA (gns(57,i),i=1,2)/11,8/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.9071292d0, 137.905991d0,
     1    139.9054431d0, 141.9092504d0/
      DATA (gns(58,i),i=1,4)/1,1,1,1/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076576d0/
      DATA (gns(59,i),i=1,1)/6/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077290d0, 142.9098200d0,
     1    143.9100930d0, 144.9125793d0, 145.9131226d0, 147.9168993d0,
     2    149.9209022d0/
      DATA (gns(60,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912756d0/
      DATA (gns(61,i),i=1,1)/6/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.9120065d0, 146.9149044d0,
     1    147.9148292d0, 148.9171921d0, 149.9172829d0, 151.9197397d0,
     2    153.9222169d0/
      DATA (gns(62,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198578d0, 152.9212380d0/
      DATA (gns(63,i),i=1,2)/6,6/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197995d0, 153.9208741d0,
     1    154.9226305d0, 155.9221312d0, 156.9239686d0, 157.9241123d0,
     2    159.9270624d0/
      DATA (gns(64,i),i=1,7)/1,1,4,1,4,1,1/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253547d0/
      DATA (gns(65,i),i=1,1)/4/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.9242847d0, 157.924416d0,
     1    159.9252046d0, 160.9269405d0, 161.9268056d0, 162.9287383d0,
     2    163.9291819d0/
      DATA (gns(66,i),i=1,7)/1,1,1,6,1,6,1/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303288d0/
      DATA (gns(67,i),i=1,1)/8/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.9287884d0, 163.9292088d0,
     1    165.9302995d0, 166.9320546d0, 167.9323767d0, 169.9354702d0/
      DATA (gns(68,i),i=1,6)/1,1,1,8,1,1/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342179d0/
      DATA (gns(69,i),i=1,1)/2/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.9338896d0, 169.9347664d0,
     1    170.9363302d0, 171.9363859d0, 172.9382151d0, 173.9388664d0,
     2    175.9425764d0/
      DATA (gns(70,i),i=1,7)/1,1,2,1,6,1,1/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407752d0, 175.9426897d0/
      DATA (gns(71,i),i=1,2)/6,15/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.9400461d0, 175.9414076d0,
     1    176.9432277d0, 177.9437058d0, 178.9458232d0, 179.9465570d0/
      DATA (gns(72,i),i=1,6)/1,1,8,1,10,1/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (gns(73,i),i=1,2)/17,8/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.9467108d0, 181.9482039d0,
     1    182.9502227d0, 183.9509309d0, 185.9543628d0/
      DATA (gns(74,i),i=1,5)/1,1,2,1,1/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529545d0, 186.9557501d0/
      DATA (gns(75,i),i=1,2)/6,6/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524885d0, 185.9538350d0,
     1    186.9557474d0, 187.9558352d0, 188.9581442d0, 189.9584437d0,
     2    191.9614770d0/
      DATA (gns(76,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605893d0, 192.9629216d0/
      DATA (gns(77,i),i=1,2)/4,4/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959930d0, 191.961039d0,
     1    193.9626809d0, 194.9647917d0, 195.9649521d0, 197.9678949d0/
      DATA (gns(78,i),i=1,6)/1,1,1,2,1,1/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665688d0/
      DATA (gns(79,i),i=1,1)/4/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667686d0,
     1    198.9682806d0, 199.9683266d0, 200.9703028d0, 201.9706434d0,
     2    203.9734940d0/
      DATA (gns(80,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723446d0, 204.9744278d0/
      DATA (gns(81,i),i=1,2)/2,2/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730440d0, 205.9744657d0,
     1    206.9758973d0, 207.9766525d0/
      DATA (gns(82,i),i=1,4)/1,1,2,1/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803991d0/
      DATA (gns(83,i),i=1,1)/10/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824308d0/
      DATA (gns(84,i),i=1,1)/2/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (gns(85,i),i=1,1)/11/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175782d0/
      DATA (gns(86,i),i=1,1)/1/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197360d0/
      DATA (gns(87,i),i=1,1)/4/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254103d0/
      DATA (gns(88,i),i=1,1)/1/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277523d0/
      DATA (gns(89,i),i=1,1)/4/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380558d0/
      DATA (gns(90,i),i=1,1)/1/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358842d0/
      DATA (gns(91,i),i=1,1)/4/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396355d0, 234.0409523d0,
     1    235.0439301d0, 238.0507884d0/
      DATA (gns(92,i),i=1,4)/6,1,8,1/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481736d0/
      DATA (gns(93,i),i=1,1)/6/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064205d0/
      DATA (gns(94,i),i=1,1)/1/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613815d0/
      DATA (gns(95,i),i=1,1)/6/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (gns(96,i),i=1,1)/10/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (gns(97,i),i=1,1)/4/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079589d0/
      DATA (gns(98,i),i=1,1)/2/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (gns(99,i),i=1,1)/11/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095106d0/
      DATA (gns(100,i),i=1,1)/10/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (gns(101,i),i=1,1)/17/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (gns(102,i),i=1,1)/10/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105510d0/
      DATA (gns(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (gns(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114070d0/
      DATA (gns(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118290d0/
      DATA (gns(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122970d0/
      DATA (gns(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.129793d0/
      DATA (gns(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137370d0/
      DATA (gns(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LT.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.GT.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          ENDIF
      IF((IAN.EQ.0).AND.(IMN.GT.1)) THEN
          IF(IMN.EQ.2) NAME=' d'
          IF(IMN.EQ.3) NAME=' t'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      DGNS= -1
      ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.15)  write(6,606) ian,imn,nmn(ian)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              DGNS= gns(IAN,I)
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
  606  format(/' *** ERROR *** called MASSES for atom with  AN=',I4,
     1  '  MN=',I4,'n(MN)=',I4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE READATA(NSTATES,PASok,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,
     1                                            NDAT,NOWIDTHS,PRINP)
c***********************************************************************
c** Subroutine to read, do book-keeping for, and print summary of
c  experimental data used in fits to spectroscopic data for one or more
c  electronic states and one or more isotopomers. 
c             ********* Version of 4 April 2016 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 1997-2016 by  Robert J. Le Roy & Dominique R.T. Appadoo +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The present program version can treat seven types of experimental
c   experimental data, for up to NISTPMX isotopomers of a given species.
c   The data are read in grouped as "bands", as (fluorescence) series, 
c   as binding energies (from photoassociation spectroscopy), as a set
c   of Bv values for a given electronic state, and [in a potential-fit
c   aanalysis] as tunneling predissociation level widths.  The types are
c   identified by the values of the 'electronic state label' parameters
c   IEP & IEPP.  They are:
c (i)  microwave transitions within a given electronic state;
c (ii)  infrared bands among the vibrational levels a given state;
c (iii) fluorescence series from some initial excited state level into 
c    vibration-rotation levels of a given electronic state
c (iv)  visible (electronic) absorption or emission bands between vib.
c    levels of two electronic state.
c (v)  binding energies - as from photoassociation spectroscopy
c (vi) "experimental" B_v values for vibrational levels of one of the
c    electronic states.
c (vii) Widths of tunneling predissociation quasibound levels (this 
c    option only meaningful for program DSPotFit).  
c-----------------------------------------------------------------------
c** On Entry:
c  NSTATES is the number of electronic states involved in the data set
c    considered (don't count states giving rise to fluorescence series).
c  PASok indicates how photoassociation data to be treated in analysis:
c    If(PASok(ISTATE).GE.1) treat it as proper PA binding energy data.
c    If(PASok(ISTATE).LE.0) treat PAS data as fluorescence series.
c    Set PASok= 0 if potential model has no explicit Dissoc. Energy
c  Data cutoffs:  for levels of electronic state  s , neglect data with:
c     J(s) > JTRUNC(s),  or vibrational levels lying outside the range
c     VMIN(s)  to  VMAX(s),  AND  NEGLECT any data for which the read-
c     in uncertainty is  > UCUTOFF (cm-1).  EFSEL(s) > 0 causes f-parity
c     levels to be neglected, EFSEL(s) < 0 omits e-parity levels
c     while  EFSEL(s) = 0  allows both types of parity to be included.
c  NOWIDTHS > 0  causes the program to ignore any tunneling widths in
c            the data set.
c  PRINP > 0  turns on the printing of a summary description of the data.
c** On Return:
c  UCUTOFF (cm-1)  is the smallest uncertainty in the (accepted) data
c  NDAT(v,i,s)  is the number of transitions associated with 
c    vibrational level-v of isotopomer-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in the experimental data on channel-4
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      INTEGER I,IBB,NTRANS,COUNT,IBAND,JMAX(NBANDMX),JMIN(NBANDMX),
     1  VMX(NSTATEMX),ISOT,NBND,ESP,ESPP,ISTATE,ISTATEE,MN1,MN2,PRINP,
     2  FSOMIT,VMAXesp,VMINesp,VMAXespp,VMINespp,JTRUNCesp,JTRUNCespp
      INTEGER NSTATES,NOWIDTHS,JTRUNC(NSTATEMX),EFSEL(NSTATEMX),
     1  VMIN(NSTATEMX),VMAX(NSTATEMX),NDAT(0:NVIBMX,NISTPMX,NSTATEMX),
     2  PASok(NSTATES)
      REAL*8 UCUTOFF,UMIN,TOTUFREQ
      CHARACTER*3 NEF(-1:1)
      CHARACTER*3 LABLP,LABLPP
c
c** Type statements & common block for data
cc
cc    REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
cc   1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
cc   2  RMUP(0:9,NISTPMX)
cc    INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
cc   1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
cc   2  EFP(NDATAMX),EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),
cc   3  FSBAND(NBANDMX),NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),
cc   4  ISTP(NBANDMX),IFIRST(NBANDMX),ILAST(NBANDMX),
cc   5  NTV(NSTATEMX,NISTPMX)
cc    CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX)
cc    COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,COUNTOT,
cc   1 NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,TVUP,TVLW,VP,VPP,
cc   2 FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV, NAME,SLABL
c
c** Type statements & common blocks for characterizing transitions
c
      REAL*8  AVEUFREQ(NBANDMX),MAXUFREQ(NBANDMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NBVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NBANDS(NISTPMX),
     7  YPR(NISTPMX,NSTATEMX,7,6,NBANDMX)
c
      COMMON /TYPEBLK/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NBVPP,NWIDTH,
     2  NEBPAS,NBANDS,YPR
c
      DATA NEF/'  f','   ','  e'/
c-----------------------------------------------------------------------
      WRITE(6,603) UCUTOFF 
      DO  ISTATE= 1,NSTATES
          IF(JTRUNC(ISTATE).GE.0) THEN
              WRITE(6,607) SLABL(ISTATE),JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ELSE
              WRITE(6,605) SLABL(ISTATE),-JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ENDIF 
          IF(EFSEL(ISTATE).GT.0) WRITE(6,601) NEF(-1)
          IF(EFSEL(ISTATE).LT.0) WRITE(6,601) NEF(1)
          ENDDO
      UMIN= UCUTOFF
c** Initialize counters for book-keeping on input data
      COUNT= 0
      DO  ISOT= 1,NISTP
          DO  ISTATE= 1,NSTATES
              NTRANSFS(ISOT,ISTATE)= 0
              NTRANSIR(ISOT,ISTATE)= 0
              NTRANSMW(ISOT,ISTATE)= 0
              NBANDFS(ISOT,ISTATE)= 0
              NBANDVIS(ISOT,ISTATE)= 0
              NBANDIR(ISOT,ISTATE)= 0
              NBANDMW(ISOT,ISTATE)= 0
              NBVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              NEBPAS(ISOT,ISTATE)= 0
              DO  I= 1,NSTATES
                  NTRANSVIS(ISOT,ISTATE,I)= 0
                  NBANDEL(ISOT,ISTATE,I)= 0
                  ENDDO
              ENDDO
          NBANDS(ISOT)= 0
          ENDDO
      DO  ISTATE= 1,NSTATES
          VMX(ISTATE)= 0
          ENDDO
      NFSTOT= 0
      FSOMIT= 0
c========================================================================
c** Begin loop to read in data, band(or series)-by-band(or series).
c  STOP when run out of bands or when encounter a negative vibrational
c  quantum number.
c** Read all data for each isotopomer at one time.
      IBAND= 0
   10 CONTINUE
      IBAND= IBAND+1
      IF(IBAND.GT.NBANDMX) THEN
            IF(PRINP.GT.0) WRITE(6,609) IBAND,NBANDMX
            IBAND= IBAND-1
            GOTO 20
            ENDIF
c
c For each "band", read in:  (i) upper/lower vibrational quantum numbers
c   VP & VPP,  (ii) a two-character electronic-state alphameric label 
c   {enclosed in single quotes; e.g., 'X0' or 'A1'} for the upper
c   (LABLP) and lower (LABLP) state, and  (iii) integers NM1 & NM2 are
c   the mass numbers [corresponding to input atomic numbers AN(1) & 
c   AN(2)] identifying the particular isotopomer.  Note that LABLP also
c   identifies the type of data in the 'band' or data-group (see below).
c
c** LABLP = LABLPP  and  VP = VPP  for a microwave band
c   LABLP = LABLPP  and  VP.ne.VPP  for an infrared band 
c   LABLP = 'FLS'  identifies this data group/band as a fluorescence 
c           series from a single emitting level into vibrational levels
c           of electronic state LABLPP.  In this case: VP is the quantum
c           number v' for the emitting level, while VPP is actually the 
c           rotational quantum number J' for the emitting level and JP
c           [see below] the lower state vibrational quantum number v".
c   LABLP = 'PAS'  identifies this data group/band as a set of binding
c           energies [D-E(v,J,p)] for a given state.  Labels as for 'FLS'
c   LABLP = 'BVV'  identifies this data group/band as a set of Bv values
c           for electronic state LABLPP.  In this case, parameters  VP
c           & VPP are dummy variables, as are EFP, JPP and EFPP [see
c           below],  JP is actually the vibrational quantum number v",
c           FREQ the Bv value & UFREQ its uncertainty
c   LABLP = 'WID'  identifies this data group/band as a set of tunneling 
c           predissociation widths for electronic state LABLPP.  In this
c           case, parameters VP, VPP and EFP are dummy variables, while
c           the predissociating level is identified as: v"=JP, J"=JPP,
c           and parity p"=EFPP.
c   NOTE: !!!!!!!!!!! This last option is ignored by DSParFit !!!!!!!!!
c** STOP reading when run out of bands OR when read-in VPP is negative   
c-----------------------------------------------------------------------
      READ(4,*,END=20) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1,MN2
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 20
      IEP(IBAND)= -99
      IEPP(IBAND)= -99
      DO  I= -3,NSTATES
          IF(LABLP.EQ.SLABL(I)) IEP(IBAND)= I
          IF(LABLPP.EQ.SLABL(I)) IEPP(IBAND)= I
          ENDDO
c** Check that this isotopomer is one of those chosen to be fitted ...
      ISOT= 0
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      TOTUFREQ= 0.D0
      MAXUFREQ(IBAND)= 0
      JMAX(IBAND)= 0
      JMIN(IBAND)= 9999
      COUNT= COUNT+1
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      NTRANS= 0
      IFIRST(IBAND)= COUNT
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
      IF((ESPP.GT.0).AND.(ISOT.GT.0)) THEN
          VMAXespp= VMAX(ESPP)
          VMINespp= VMIN(ESPP)
          JTRUNCespp= JTRUNC(ESPP)
          IF(ISOT.GT.1) THEN
              VMAXespp= INT((VMAX(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
              VMINespp= INT((VMIN(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
              JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
              ENDIF
cc        VMAXesp= VMAX(ESPP)      ?????? huh?
          ENDIF
      IF((ESP.GT.0).AND.(ISOT.GT.0)) THEN
          VMAXesp= VMAX(ESP)
          VMINesp= VMIN(ESP)
          JTRUNCesp= JTRUNC(ESP)
          IF(ISOT.GT.1) THEN
              VMAXesp= INT((VMAX(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              VMINesp= INT((VMIN(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              JTRUNCesp= INT(JTRUNC(ESP)/RSQMU(ISOT))
              ENDIF
          ENDIF
c** For each of the lines in a given band/series, read upper level
c  rotational quantum number (JP) and e/f parity [EFP= +1 for e, -1 for
c  f, and  0 if e/f splitting unresolved and to  be ignored], and lower
c  level rotational quantum number (JPP) and parity [EFPP, as above],
c  the transition frequency  FREQ, and its uncertainty UFREQ.
c** For PAS or Tunneling Width data,  JP(COUNT)=v", JPP(COUNT)=J", 
c  EFPP(COUNT)=p", FREQ is the observable (a positive No.), while 
c  EFP(COUNT), VP(IBAND) & VPP(IBAND) are dummy variables.
c** For Bv values, JP(COUNT)=v" while JPP(COUNT), EFP(COUNT) and
c   EFPP(COUNT) as well as VP(IBAND) & VPP(IBAND) are dummy variables.
c-----------------------------------------------------------------------
   15 READ(4,*) JP(COUNT), EFP(COUNT), JPP(COUNT), EFPP(COUNT), 
     1                                       FREQ(COUNT), UFREQ(COUNT)
c-----------------------------------------------------------------------
c=======================================================================
c   Sample IR band data of HF for the '.4' file:
c   --------------------------------------------                          
c   1 0  'X0' 'X0'  1 19             % VP VPP LABLP LABLPP MN1 MN2
c   8 1   9 1  266.0131002  0.005    % JP EFP JPP EFPP FREQ UFREQ
c   9 1  10 1  265.8885896  0.003
c  10 1  11 1  265.7716591  0.002
c   .    .      .            .
c   .    .      .            .
c   [end of a band indicated by -ve JP and/or JPP value(s)]
c  -1 1  -1 1  -1.1         -1.1
c=======================================================================
      IF(EFP(COUNT).GT.1) EFP(COUNT)= 1
      IF(EFP(COUNT).LT.-1) EFP(COUNT)= -1
      IF(EFPP(COUNT).GT.1) EFPP(COUNT)= 1
      IF(EFPP(COUNT).LT.-1) EFPP(COUNT)= -1
c** At end of a band, exit from implicit loop
      IF((JPP(COUNT).LT.0).OR.(JP(COUNT).LT.0)) GOTO 18
c** If this band is not for one of the isotopomers chosen to be fitted,
c  omit its data from the fit
      IF(ISOT.EQ.0) GO TO 15
c** If this band involves electronic states other than those chosen to 
c   be treated, omit its data from the fit
      IF((ESP.EQ.-99).OR.(ESPP.EQ.-99)) GO TO 15
c** If a datum uncertainty of zero is accidentally read in, STOP
      IF(DABS(UFREQ(COUNT)).LE.0.d0) THEN
          WRITE(6,600) COUNT,FREQ(COUNT),IBAND
          STOP
          ENDIF
c** Omit data  with uncertainties outside specified limit UCUTOFF
      IF(UFREQ(COUNT).GT.UCUTOFF) GOTO 15
c** Require that datum lies within specified J & v ranges
      IF(ESP.GE.-2) THEN
          IF(((JTRUNCespp.GE.0).AND.(JPP(COUNT).GT.JTRUNCespp)).OR.
     1       ((JTRUNCespp.LT.0).AND.(JPP(COUNT).LT.-JTRUNCespp)))
     2                                                         GOTO 15
          IF((EFPP(COUNT)*EFSEL(ESPP)).LT.0) GOTO 15
          ENDIF
      IF(ESP.GT.0) THEN
          IF(VPP(IBAND).GT.VMAXespp) GOTO 15
          IF(VPP(IBAND).LT.VMINespp) GOTO 15
          IF(VP(IBAND).GT.VMAXesp) GOTO 15
          IF(VP(IBAND).LT.VMINesp) GOTO 15
          IF((JTRUNCesp.GE.0).AND.(JP(COUNT).GT.JTRUNCesp)) GOTO 15
          IF((JTRUNCesp.LT.0).AND.(JP(COUNT).LT.-JTRUNCesp)) GOTO 15
          IF((EFP(COUNT)*EFSEL(ESP)).LT.0) GOTO 15
        ELSE
          IF(JP(COUNT).GT.VMAXespp) GOTO 15
          IF(JP(COUNT).LT.VMINespp) GOTO 15
        ENDIF
c** If NOWIDTHS > 0  omit any tunneling width data from the fit.
      IF((ESP.EQ.-2).AND.(NOWIDTHS.GT.0)) GOTO 15
c
c** End of tests for datum inclusion.  Now count/sort data
c=======================================================================
      TVUP(COUNT)= 0
      TVLW(COUNT)= 0
      IF(ESP.GE.-1) UMIN= MIN(UMIN,UFREQ(COUNT))
c** Determine actual v & J range of data & count data for each v
c  JMIN & JMAX needed for printout summary & data-count for testing
c  no. parameters allowed in Band Constant fit.
c??? This segment imperfect & needs re-examination ?????????????
      IF(ESP.GT.0) THEN
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
          VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
c
c** Accumulate count of data associated with each vibrational level ...
          NDAT(VPP(IBAND),ISTP(IBAND),ESPP)=
     1                            NDAT(VPP(IBAND),ISTP(IBAND),ESPP)+ 1
          NDAT(VP(IBAND),ISTP(IBAND),ESP)=
     1                              NDAT(VP(IBAND),ISTP(IBAND),ESP)+ 1
        ELSEIF((ESP.LE.0).OR.(ESP.GE.-2)) THEN
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESPP)= MAX(VMX(ESPP),JP(COUNT))
          NDAT(JP(COUNT),ISTP(IBAND),ESPP)=
     1                      NDAT(JP(COUNT),ISTP(IBAND),ESPP)+ 1
        ELSEIF(ESP.EQ.-3) THEN
c... and for Bv data ...
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          NDAT(JPP(COUNT),ISTP(IBAND),ESPP)=
     1                     NDAT(JPP(COUNT),ISTP(IBAND),ESPP)+ 1
        ENDIF
      DFREQ(COUNT)= 0.d0
      IB(COUNT)= IBAND 
      TOTUFREQ= TOTUFREQ+UFREQ(COUNT) 
      IF(UFREQ(COUNT).GT.MAXUFREQ(IBAND)) MAXUFREQ(IBAND)= UFREQ(COUNT)
      COUNT= COUNT+1 
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      GOTO 15 
c** End of loop reading data for a given band/series 
c
c** Tidy up at end of reading for a given band
   18 COUNT= COUNT-1
      ILAST(IBAND)= COUNT 
      NTRANS= ILAST(IBAND)-IFIRST(IBAND)+1
      IF(NTRANS.GT.0) THEN
c** Treat PAS data as Fluorescence series unless  PASok > 0
          IF((IEP(IBAND).EQ.-1).AND.(PASok(IEPP(IBAND)).LE.0)) 
     1                                                    IEP(IBAND)=0
          IF((NTRANS.EQ.1).AND.(LABLP.EQ.'FLS')) THEN
c** Ignore any fluorescence series consisting of only one datum
              COUNT= COUNT-1
              IBAND= IBAND-1
              FSOMIT= FSOMIT+1
              GOTO 10
              ENDIF
          AVEUFREQ(IBAND)= TOTUFREQ/NTRANS
          NBANDS(ISTP(IBAND))= NBANDS(ISTP(IBAND))+1
        ELSE
          IBAND= IBAND-1
          GOTO 10
        ENDIF
c=======================================================================
c** Accumulate counters for bands/series of different types
      IF(ESP.EQ.0) THEN
c** For Fluorescence Series ... first enumerate the No. of bands & lines
          NFSTOT= NFSTOT+1
          FSBAND(NFSTOT)= IBAND
c** Define counter to label which f.s. is associated with band IBAND 
          NFS(IBAND)= NFSTOT
          NBANDFS(ISOT,ESPP)= NBANDFS(ISOT,ESPP)+1
          NBND= NBANDFS(ISOT,ESPP)
          NTRANSFS(ISOT,ESPP)= NTRANSFS(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each band
          YPR(ISOT,ESPP,1,1,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,1,2,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,1,3,NBND)= NTRANS
          YPR(ISOT,ESPP,1,4,NBND)= IBAND
          YPR(ISOT,ESPP,1,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,1,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.GT.0).AND.(ESP.NE.ESPP)) THEN
c** For vibrational band of a normal 2-state electronic transition
c ... count bands and transitions in visible (electronic) spectrum
          NBANDEL(ISOT,ESP,ESPP)= NBANDEL(ISOT,ESP,ESPP)+ 1
          NBANDVIS(ISOT,ESPP)= NBANDVIS(ISOT,ESPP)+ 1
          NBND= NBANDVIS(ISOT,ESPP)
          NTRANSVIS(ISOT,ESP,ESPP)= NTRANSVIS(ISOT,ESP,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,2,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,2,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,2,3,NBND)= NTRANS
          YPR(ISOT,ESPP,2,4,NBND)= IBAND
          YPR(ISOT,ESPP,2,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,2,6,NBND)= JMAX(IBAND)
          ENDIF 
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).NE.VPP(IBAND))) THEN
c** For an Infrared band of electronic state  s=ESPP=ESP
c** First cumulatively count the number of IR bands & transitions
          NBANDIR(ISOT,ESPP)= NBANDIR(ISOT,ESPP)+1
          NBND= NBANDIR(ISOT,ESPP)
          NTRANSIR(ISOT,ESPP)= NTRANSIR(ISOT,ESPP)+NTRANS 
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,3,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,3,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,3,3,NBND)= NTRANS
          YPR(ISOT,ESPP,3,4,NBND)= IBAND
          YPR(ISOT,ESPP,3,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,3,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).EQ.VPP(IBAND))) THEN
c** For Microwave transitions in electronic state  s=ESPP=ESP
c** First cumulatively count the number of MW bands & transitions
          NBANDMW(ISOT,ESPP)= NBANDMW(ISOT,ESPP)+1
          NBND= NBANDMW(ISOT,ESPP)
          NTRANSMW(ISOT,ESPP)= NTRANSMW(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,4,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,4,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,4,3,NBND)= NTRANS
          YPR(ISOT,ESPP,4,4,NBND)= IBAND
          YPR(ISOT,ESPP,4,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,4,6,NBND)= JMAX(IBAND)
          ENDIF
c
c** NOTE ... in YPR array a last index counts bands of this type for 
c  this isotopomer of this electronic state ... and put all Bv's, 
c  Tunneling Widths or PAS binding energies in one group.
      IF(ESP.EQ.-3) THEN
c** Data are not transition energies, but rather the values of Bv in
c  electronic state s=IEPP  [As in the published IBr(A-X) analysis].
ccc       IF((NBVPP(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
              WRITE(6,612) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NBVPP(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,5,3,1)= NTRANS
          YPR(ISOT,ESPP,5,4,1)= IBAND
          YPR(ISOT,ESPP,5,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,5,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-2) THEN
c** Data are tunneling predissociation linewidths (in cm-1) for levels
c  of electronic state IEPP=ESPP
ccc       IF((NWIDTH(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
ccc           WRITE(6,626) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NWIDTH(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,6,3,1)= NTRANS
          YPR(ISOT,ESPP,6,4,1)= IBAND
          YPR(ISOT,ESPP,6,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,6,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-1) THEN
c** Data are PhotoAssociation Binding Energies (in cm-1) for levels
c  of electronic state IEPP=ESPP
          WRITE(6,636) LABLPP,ISOT
          NEBPAS(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,7,3,1)= NTRANS
          YPR(ISOT,ESPP,7,4,1)= IBAND
          YPR(ISOT,ESPP,7,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,7,6,1)= JMAX(IBAND)
          ENDIF
c** Now return to read the next band
      GOTO 10
c========================================================================
c** Now, write a summary of the input data to the output file
   20 COUNTOT= COUNT
      NBANDTOT= 0
      DO  I= 1,NISTP
          NBANDTOT= NBANDTOT+ NBANDS(I)
          ENDDO
      ISOT= 1
      UCUTOFF= UMIN
      IF(FSOMIT.GT.0) WRITE(6,650) FSOMIT
      IF(PRINP.LE.0) RETURN
c** Print a summary of the data, one isotopomer at a time.
   26 WRITE(6,602) NBANDS(ISOT), (NAME(I),MN(I,ISOT),I=1,2)
c
      DO 50 ISTATE= 1,NSTATES
c ... For internal use, may wish to update VMAX(ISTATE) to the actual 
c  highest v in the data set for this state. ** Reactivate as needed.
c      VMAX(ISTATE)= VMX(ISTATE)
c ... and separately list data for each (lower) electronic state in turn
      IF(NTRANSMW(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDMW(ISOT,ISTATE)
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,4,2,I),
     1                     YPR(ISOT,ISTATE,4,1,I),
     2                  YPR(ISOT,ISTATE,4,3,I),YPR(ISOT,ISTATE,4,5,I),
     3                  YPR(ISOT,ISTATE,4,6,I), 
     3                  AVEUFREQ(YPR(ISOT,ISTATE,4,4,I)),
     4                  MAXUFREQ(YPR(ISOT,ISTATE,4,4,I))
              ENDDO
          ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDIR(ISOT,ISTATE)
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                     YPR(ISOT,ISTATE,3,1,I),
     2                  YPR(ISOT,ISTATE,3,3,I),YPR(ISOT,ISTATE,3,5,I),
     3                  YPR(ISOT,ISTATE,3,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,3,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,3,4,I))
              ENDDO
          ENDIF
c
c** Book-keeping for electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1         (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),SLABL(ISTATE),
     2                                    NBANDEL(ISOT,ISTATEE,ISTATE)
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                            YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3                  YPR(ISOT,ISTATE,2,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,2,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,2,4,I))
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,614) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              WRITE(6,616) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,1,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,1,4,I))
              ENDDO
          ENDIF
      IF(NBVPP(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for  Bv  data
          WRITE(6,618) NBVPP(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,5,4,1)
          WRITE(6,620) YPR(ISOT,ISTATE,5,3,1),YPR(ISOT,ISTATE,5,5,1),
     1       YPR(ISOT,ISTATE,5,6,1),AVEUFREQ(YPR(ISOT,ISTATE,5,4,1)),
     2                               MAXUFREQ(YPR(ISOT,ISTATE,5,4,1))
          ENDIF
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          WRITE(6,628) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,6,3,1),
     1                  YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,6,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,6,4,1))
          ENDIF
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS Binding Energy  data
          WRITE(6,632) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,7,3,1),
     1                  YPR(ISOT,ISTATE,7,5,1),YPR(ISOT,ISTATE,7,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,7,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,7,4,1))
          ENDIF
   50 CONTINUE
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 26
          ENDIF 
      WRITE(6,622)
      RETURN
  600 FORMAT(/' *** INPUT ERROR ***  Datum   FREQ(',i5,')=',f12.4,
     1 '  in   IBAND=',i4,'   has zero uncertainty!!!')
  601 FORMAT(23x,'or with',A3,'-parity.')
  603 FORMAT(/' Neglect data with:  Uncertainties > UCUTOFF=',G12.3,
     1  ' (cm-1)')
  605 FORMAT(7x,'and State ',A3,' data with  J < JTRUNC=',I4,
     1  '  or  v  outside range',i3,'  to',i4)
  607 FORMAT(7x,'and State ',A3,' data with  J > JTRUNC=',I4,
     2  '  or  v  outside range',i3,'  to',i4)
  602 FORMAT(/1x,20('===')/'  *** Input data for',i5,' bands/series of '
     1  ,A2,'(',I3,')-',A2,'(',I3,') ***'/1x,20('==='))
  604 FORMAT(1x,28('--')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' sets'/1x,28('--')/"   v'  ",
     1 'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/1x,25('--'))
  606 FORMAT(I4,I4,3I7,1x,1P2D10.1)
  608 FORMAT(1x,32('--')/I6,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands'/1x,32('--')/
     2 "   v'  ",'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A3,'}--{State ',A3,'} Transitions in',i4,' Bands'/1x,35('--')/
     2 "   v'",'  v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  612 FORMAT(/" NOTE that all read-in Bv's for   ISTATE=",i2,'   ISOT=',
     1  i2/32x,' must be input as a single "band" or data group')
cc612 FORMAT(/" *** STOP INPUT *** and put all read-in Bv's for   ISTATE
cc   1=",i2,'   ISOT=',i2/ 10x,'into one "band" or data group.')
  614 FORMAT(1x,38('==')/I5,' Fluorescence transitions into State ',
     1 A3,2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series'/
     2 1x,38('==')/"   v'  j' p' ",'#data  v"min  v"max  Avge.Unc.  Max.
     3Unc.'/1x,51('-'))
  616 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1)
  618 FORMAT(1x,65('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Bv values treated as independent data'/1x,24('--')/
     2 '  #values   v(min)  v(max)  Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  620 FORMAT(I7,I9,I8,3x,1P2D11.1)
  622 FORMAT(1x,25('===')/1x,25('==='))
  626 FORMAT(/" NOTE that all read-in Tunneling Widths for   ISTATE=",
     1 i2,'   ISOT=',i2/10x,' must be in a single "band" or data group')
cc626 FORMAT(/" *** STOP INPUT *** and put all read-in Tunneling Widths'
cc   1  '  for   ISTATE=",i2,'   ISOT=',i2/ 
cc   2  10x,'into one "band" or data group.')
  628 FORMAT(1x,61('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Tunneling Widths included as data'/
     2 1x,61('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  630 FORMAT(I7,I9,I8,2x,1P2D11.1)
  632 FORMAT(1x,70('=')/I4,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') PAS Binding Energies included in data set'/
     2 1x,70('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  636 FORMAT(/' NOTE that all read-in PAS Binding Energies for   ISTATE=
     1 ',a2,'  ISOT=',i2/10x,' must be in a single "band" or data group'
     2 )
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
  650 FORMAT(/' Data input IGNORES',i4,' fluorescence series consisting'
     1 ,' of only  onee  line!')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE TVSORT(ISTATE,NPARM,VMAX,NTVALL)
c***********************************************************************
c** Subroutine to sort through global data file, and for each isotopomer
c  in state ISTATE:  (1) find the number of transitions coupled to each
c  level (v,J,p),  (2) for levels in order (v,J,p), add a free parameter
c  for each level involved in one or more transitions, and  (3) label each
c  transition involving one of these levels by the index/counter of the
c  parameter associated with that term value.
c             ********* Version of 27 August 2004 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On Entry:
c------------
c  ISTATE is the electronic state being considered.
c  NPARM  enters as the cumulative count of parameters prior to entry
c  TVUP(i) and TVLW(i) in COMMON equal zero for all data
c
c** On Return:
c-------------
c  NPARM  is updated to include the number of term values for this state
c  TVUP(i) & TVLW(i): if the upper and/or lower level of transition-i is
c      to be represented by a term value, TVUP and TVLW (respectively)
c      is the associated parameter index; otherwise they = 0.
c
c** Internally
c-------------
c  NLV(v,J.p) * initially, counts transitions for level {v,J,p} of a 
c                          given isotopologue
c           * later reset it as the parameter index for that term value
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      INTEGER I,J,P,IBAND,ISOT,ISTATE,NPARM,LOWEST,VMAX(NSTATEMX),
     1 NLV(0:NVIBMX,0:NROTMX,-1:1),NTVS(NSTATEMX,NISTPMX),
     2 NTVALL(0:NSTATEMX)
c
c** Type statements & common block for data
cc
cc    REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
cc   1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
cc   2  RMUP(0:9,NISTPMX)
cc    INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
cc   1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
cc   2  EFP(NDATAMX),EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),
cc   3  FSBAND(NBANDMX),NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),
cc   4  ISTP(NBANDMX),IFIRST(NBANDMX),ILAST(NBANDMX),
cc   5  NTV(NSTATEMX,NISTPMX)
cc    CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX)
cc    COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,COUNTOT,
cc   1 NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,TVUP,TVLW,VP,VPP,
cc   2 FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
      WRITE(6,600) SLABL(ISTATE) 
      LOWEST= 1
      IF(ISTATE.GT.1) LOWEST= 0
      NTVALL(ISTATE)= 0
      DO  ISOT= 1, NISTP
c** First ... zero transition counter array for this isotopomer
          DO  I= 0, VMAX(ISTATE)
              DO  J= 0, NROTMX 
                  DO  P= -1,1
                      NLV(I,J,P)= 0
                      ENDDO
                  ENDDO
              ENDDO
          DO  IBAND= 1, NBANDTOT
c** Then ... search for bands involving isotopomer ISOT in this state
              IF(((IEP(IBAND).EQ.ISTATE).OR.(IEPP(IBAND).EQ.ISTATE))
     1          .AND.(ISTP(IBAND).EQ.ISOT).AND.(IEP(IBAND).GE.0)) THEN
                  DO  I= IFIRST(IBAND), ILAST(IBAND)
c ... for each such band, loop over all transitions, and increment NLV 
c     for each {v,J,p} level encountered in a transision
                      IF(IEP(IBAND).EQ.ISTATE) THEN
                          IF(JP(I).GT.NROTMX) THEN
c ... check for array dimension overruns
                              WRITE(6,602) ISTATE,ISOT,JP(I),NROTMX
                              STOP
                              ENDIF
                          NLV(VP(IBAND),JP(I),EFP(I))= 
     1                                  NLV(VP(IBAND),JP(I),EFP(I))+ 1
                          ENDIF
                      IF(IEPP(ISTATE).EQ.ISTATE) THEN
                          IF(JPP(I).GT.NROTMX) THEN
                              WRITE(6,604) ISTATE,ISOT,JPP(I),NROTMX
                              STOP
                              ENDIF
                          NLV(VPP(IBAND),JPP(I),EFPP(I))
     1                             = NLV(VPP(IBAND),JPP(I),EFPP(I))+ 1
                          ENDIF
                      ENDDO
                  ENDIF
c** Finished scan over all data set for this isotopologue
              ENDDO
c
c** Now ... count a free parameter for each level in a transition
c** NTV  is the total number of term values for case (ISTATE,ISOT) 
c   NTVS is the no. of them involved in only a single transition
          NTV(ISTATE,ISOT)= 0
          NTVS(ISTATE,ISOT)= 0
          DO  I= 0, VMAX(ISTATE)
              DO  J= 0, NROTMX
                  DO  P= -1,1
                      IF(NLV(I,J,P).GT.0) THEN
                          IF(LOWEST.EQ.1) THEN
c** If using term values for `lowest' state (defined as the first state
c  considered), its lowest observed level for isotopologue-1 defines the
c   absolute energy zero
                              WRITE(6,606) I,J,P,ISOT,SLABL(ISTATE)
                              LOWEST= 0
                              NLV(I,J,P)= 0
                              GOTO 20
                              ENDIF
                          NPARM= NPARM+ 1
                          NTV(ISTATE,ISOT)= NTV(ISTATE,ISOT)+ 1
                          IF(NLV(I,J,P).EQ.1) NTVS(ISTATE,ISOT)=
     1                                           NTVS(ISTATE,ISOT) +1 
                          WRITE(7,700) SLABL(ISTATE),I,J,P,ISOT
c ... reset NLV(v,J,p) as the parameter index for that term value
                          NLV(I,J,P)= NPARM
                          ENDIF
   20                 CONTINUE
                      ENDDO
                  ENDDO
              ENDDO
c
c** Finally - label each transition with term-value parameter index for
c   (as appropriate) upper & lower level of each transition
          DO  IBAND= 1, NBANDTOT
              IF(((IEP(IBAND).EQ.ISTATE).OR.(IEPP(IBAND).EQ.ISTATE))
     1          .AND.(ISTP(IBAND).EQ.ISOT).AND.(IEP(IBAND).GE.0)) THEN
c ... for each band involving state ISTATE of this isotopologue, label 
c     each transition with the term value parameter index (which is zero
c     if the state is not represented by term values!).
                  DO  I= IFIRST(IBAND), ILAST(IBAND)
                      IF(IEP(IBAND).EQ.ISTATE) 
     1                            TVUP(I)= NLV(VP(IBAND),JP(I),EFP(I))
                      IF(IEPP(IBAND).EQ.ISTATE) 
     1                          TVLW(I)= NLV(VPP(IBAND),JPP(I),EFP(I))
                      ENDDO
                  ENDIF
              ENDDO
          WRITE(6,608) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                              NTV(ISTATE,ISOT),NTVS(ISTATE,ISOT)
          NTVALL(ISTATE)= NTVALL(ISTATE)+ NTV(ISTATE,ISOT)
          ENDDO
      RETURN
  600 FORMAT(/' For State ',A3,'  fit to individual term values for each
     1  {v,J,p,isot}'/1x,6('******'))
  602 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JP(ISTATE=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NROTMX=',i4)
  604 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JPP(ISTATE=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NROTMX=',i4)
  606 FORMAT(/'  Absolute zero of energy is fixed at level {v=',i3,
     1 ', J=',i3,', p=',i2,'}'/1x,12('**'),10x,'of isotopomer ',i2,
     2 ' of  State ',A3)
  608 FORMAT(' For ',A2,'(',i3,')-',A2,'(',I3,')  fit to',i4,
     1 ' T(v,J,p) term values,'/20x,'of which',i4,' are involved in only
     2 one transition')
  700 FORMAT('Tv(',A2,':',i3,',',i3,',',SP,i2,';',SS,i2,')=')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MKPREDICT(NSTATES,NDAT)
c***********************************************************************
c** Subroutine to prepare fake input data array which will cause ParFit
c  to make transition energy predictions for electronic or infrared band
c  or microwave transitions.  On entry:
c  NSTATES  is the number of states involved in the data set.  
c    NSTATES= 1 generates infrared or microwave bands for state SLABL(1)
c    NSTATES= 2 generates electronic bands from lower state SLABL(1) 
c               into upper state SLABL(2)
c  VMIN(s) and VMAX(s) are the bounds on the vibrational energy range 
c      for state 's' specified in the main input file.
c** On return:
c  NDAT(v,i,s)  is the number of transitions associated with
c    vibrational level-v of isotopomer-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in band specifications on Channel-5 and writes
c   the transition energy specifications to channel-4
c-----------------------------------------------------------------------
c                         Version of 21 August 2004
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      CHARACTER*2 LABLP,LABLPP
      INTEGER I,J,J2,JD,J2DL,J2DU,J2DD,JMAXX,PP,PPP,NTRANS,COUNT,
     1  IBAND,JMAX(NBANDMX),JMIN(NBANDMX),
     1  VMX(NSTATEMX),ISOT,ESP,ESPP,ISTATE,MN1,MN2
      INTEGER NSTATES,NDAT(0:NVIBMX,NISTPMX,NSTATEMX)
c
c** Type statements & common block for data
cc
cc    REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
cc   1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
cc   2  RMUP(0:9,NISTPMX)
cc    INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
cc   1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
cc   2  EFP(NDATAMX),EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),
cc   3  FSBAND(NBANDMX),NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),
cc   4  ISTP(NBANDMX),IFIRST(NBANDMX),ILAST(NBANDMX),
cc   5  NTV(NSTATEMX,NISTPMX)
cc    CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX),LABLP,LABLPP
cc    COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,COUNTOT,
cc   1 NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,TVUP,TVLW,VP,VPP,
cc   2 FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV, NAME,SLABL
c
c** Type statements & common blocks for characterizing transitions
c
      REAL*8  AVEUFREQ(NBANDMX),MAXUFREQ(NBANDMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NBVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NBANDS(NISTPMX),
     7  YPR(NISTPMX,NSTATEMX,7,6,NBANDMX)
c
      COMMON /TYPEBLK/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NBVPP,NWIDTH,
     2  NEBPAS,NBANDS,YPR
c-----------------------------------------------------------------------
c** Initialize counters for book-keeping on input data
      COUNT= 0
      DO  ISOT= 1,NISTP
          DO  ISTATE= 1,NSTATES
              NTRANSFS(ISOT,ISTATE)= 0
              NTRANSIR(ISOT,ISTATE)= 0
              NTRANSMW(ISOT,ISTATE)= 0
              NBANDFS(ISOT,ISTATE)= 0
              NBANDVIS(ISOT,ISTATE)= 0
              NBANDIR(ISOT,ISTATE)= 0
              NBANDMW(ISOT,ISTATE)= 0
              NBVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              DO  I= 1,NSTATES
                  NTRANSVIS(ISOT,ISTATE,I)= 0
                  NBANDEL(ISOT,ISTATE,I)= 0
                  ENDDO
              ENDDO
          NBANDS(ISOT)= 0
          ENDDO
      DO  ISTATE= 1,NSTATES
          VMX(ISTATE)= 0
          ENDDO
      NFSTOT= 0
      IBAND= 0
   10 IBAND= IBAND+ 1
      IF(IBAND.GT.NBANDMX) THEN
            WRITE(6,609) IBAND,NBANDMX
            IBAND= IBAND-1
            GOTO 99
            ENDIF
c** Generate "empty" band data sets to allow ParFit to make predictions
c  for those sets of transitions.  
c** LABLP & LABLPP are the two-character variables identifying the upper
c     and lower electronic states, respectively.  LABLP=LABLPP for IR or
c     MW transitions within a given electronic state
c** VP & VPP are the v' & v" values identifying the band;
c** PP & PPP specify the rotational parity of the upper and lower levels
c** MN1 & MN2 identify the isotopomer
c** Generate 'lines' for  J"= 0 to JMAXX subject to selection rule that
c  Delta(J) runs from J2DL to J2DU in steps of J2DD
c-----------------------------------------------------------------------
      READ(5,*,end=99) VP(IBAND),VPP(IBAND),LABLP,LABLPP,MN1,MN2,PP,PPP,
     1                                            JMAXX,J2DL,J2DU,J2DD
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 99
c** Set electronic state number for upper & lower levels.  
c* Always set lower state as 1'st state considered in input [SLABL(1)]
c* For NSTATES= 1, upper state is the same one.  For NSTATES= 2 the 
c  upper state is 2'nd one considered [SLABL(2)]
      IEPP(IBAND)= 1
      IEP(IBAND)= NSTATES
      WRITE(4,400) VP(IBAND),VPP(IBAND),LABLP,LABLPP,MN1,MN2
      ISOT= 0
c** Determine the correct isotopomer-number for this band.
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      MAXUFREQ(IBAND)= 0
      JMAX(IBAND)= JMAXX
      JMIN(IBAND)= 0
      NTRANS= 0
      IFIRST(IBAND)= COUNT+ 1
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
c** Now - loop over J to generate all possible transitions ...
      DO  J= 0, JMAXX
          DO  JD= J2DL, J2DU, J2DD
              J2= J+ JD
              IF((J2.GE.0).AND.((J.NE.0).OR.(J2.NE.0))) THEN
                  COUNT= COUNT+1
                  IF(COUNT.GT.NDATAMX) THEN
                      WRITE(6,640) COUNT,NDATAMX
                      STOP
                      ENDIF
                  WRITE(4,402) J2,PP,J,PPP
                  JP(COUNT)= J2
                  EFP(COUNT)= PP
                  JPP(COUNT)= J
                  EFPP(COUNT)= PPP
                  FREQ(COUNT)= 0.d0
                  UFREQ(COUNT)= 0.001d0
                  DFREQ(COUNT)= 0.d0
                  IB(COUNT)= IBAND
c** Accumulate count of data associated with each vibrational level ...
                  NDAT(VPP(IBAND),ISTP(IBAND),ESPP)=
     1                            NDAT(VPP(IBAND),ISTP(IBAND),ESPP)+ 1
                  NDAT(VP(IBAND),ISTP(IBAND),ESP)=
     1                              NDAT(VP(IBAND),ISTP(IBAND),ESP)+ 1
                  ENDIF
              ENDDO
          ENDDO
      WRITE(4,404)
  400 FORMAT(2I4,"   '",A2,"'  '",A2,"'   ",2I4)
  402 FORMAT(I4,I3,I5,I3,'    0.d0     1.0d-3')
  404 FORMAT('  -1 -1   -1 -1    -1.d0    -1.d-3'/)
      VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
      VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
      ILAST(IBAND)= COUNT
      NTRANS= ILAST(IBAND)-IFIRST(IBAND)+1
      GOTO 10
   99 RETURN
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NDEXM(ISTATE,SLABL,NLR,DLIMIT,VD,CN,NSIG,ZMU,NDEORD,
     1                                            NISTP,RSQMU,PNDE,XM)
c** If using NDE functions to represent vib. energies & rot. and/or CDC
c  constants, use NLR and initial trial CN value to calculate limiting
c  ND-theory scaling factor XM(IORDR,ISTATE,ISOT) of orders IORDR up to
c  NDEORD, with values for ISOT=1 rounded off to NSIG signif. digits.
c  If  NSIG.le.0  do NO rounding.                 Version date: 21/08/04
c***********************************************************************
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      CHARACTER*2 SLABL
      INTEGER I,M,ISTATE,NLR,NSIG,NISTP,NDEORD,ISOT
      REAL*8 XBAR(3:6,0:9),PNDE(0:9,NSTATEMX),XM(0:9,NSTATEMX,NISTPMX),
     1  RSQMU(NISTPMX),ZMU,DLIMIT,VD,CN,ZN,PWNDE,FCT,FCT1,FCT2
c
      DATA XBAR/36410.D0,   13433.D0,   9170.9D0,   7932.0D0,
     1      60221.03D0, 4275.748D0, 1178.287D0, 546.6391D0,
     2       -24901.D0,  -205.65D0,  -15.377D0,  -2.7539D0,
     3        1348.0D0,  -3.7691D0, -1.7742D-1, -1.9186D-2,
     4       -73.367D0, -2.9067D-1, -5.3435D-3, -2.9942D-4,
     5       -60.674D0, -.29837D-1, -.21680D-3, -.62117D-5,
     6       -38.694D0, -.35686D-2, -.10281D-4, -.15002D-6,
c ** RJL's guesses for limiting N-D theory constants defining Ov & Pv
c    7          -20.D0,    -4.6D-4,    -5.2D-7,   -0.39D-8,
c Fudged  Ov for n=5
     7          -20.D0,    -4.6D-4,    -4.151888D-7,   -0.39D-8,
     8          -10.d0,    -5.5d-5,    -2.4d-8,  -0.94d-10, 
     9          -0.d0,      0.d0,       0.d0,      0.d0/ 
c-----------------------------------------------------------------------
      ZN= DFLOAT(NLR)
      PWNDE= 2.d0*ZN/(ZN- 2.d0)
      FCT= 1.d0/(ZMU**NLR *CN**2)**(1.d0/(ZN-2.d0))
      DO  M= 0,NDEORD
          PNDE(M,ISTATE)= PWNDE
          PWNDE= PWNDE- 2.d0
          FCT2= XBAR(NLR,M)*FCT
          XM(M,ISTATE,1)= FCT2
          IF(NSIG.GT.0) THEN
c** If desired, round off the XM constant for Isotopomer-1 to NSIG digits
              FCT1= 1.d0
              IF(FCT2.LT.0.d0) FCT1= -FCT1
              FCT2= DABS(FCT2)
              I= IDINT(DLOG10(FCT2))+1
              IF(FCT2.LT.1.d0) I= I-1
              XM(M,ISTATE,1)= FCT1*DFLOAT(IDINT(FCT2*
     1                         10.d0**(NSIG-I)+0.5d0))*10.d0**(I-NSIG)
              ENDIF
          IF(NISTP.GT.1) THEN
              FCT2= 2.d0*ZN/(ZN- 2.d0)
              DO  ISOT= 2,NISTP
                  XM(M,ISTATE,ISOT)= XM(M,ISTATE,1)*RSQMU(ISOT)**FCT2
                  ENDDO
              ENDIF
          ENDDO
      WRITE(6,610) SLABL,CN,NLR,DLIMIT,VD
      DO  M= 0, NDEORD
          WRITE(6,612) (M,ISOT,XM(M,ISTATE,ISOT),ISOT=1,NISTP)
          ENDDO
      RETURN
  610 FORMAT(/' NDE function(s) for State ',A3,' which has potential tai
     1l:','   D -',1PD13.6,'/R**',i1/4x,'& initial parameters:   D(limit
     2) =',0PF12.4,'   and   '/4x,'vD(ISOT=1)=',f13.8,5x,
     3  'based on factors  XM(m,ISOT):')
  612 FORMAT(4x,3('   XM(',i1,',',i1,')=',1pD14.6:)/
     1     (31x,2('   XM(',i1,',',i1,')=',1pD14.6:)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789

c***********************************************************************
      SUBROUTINE NDEDGB(ISTATE,NISTP,NEWGv,NEWBv,RSQMU,VMAXX)
c** Subroutine to prepare various contributions to NDE partial
c  derivatives for later use in DYIDPJ;  If  NDEBv.le.0  only for Gv.
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8  VDMV,VDMVP,RSQMU(NISTPMX),GFCT,BFCT,SNUM,DNUM,SDEN,DDEN,
     1  EB,BV
      INTEGER I,ISOT,ISTATE,IV,NISTP,NEWGv,NEWBv,VMAXX
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
c
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c
c** Type statement and common block for NDE partial derivative stuff
c
      REAL*8 DGPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     1  DGQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     2  DGVD(-1:NVIBMX,NSTATEMX,NISTPMX),
     3  DBPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     4  DBQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     5  DBVD(-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /DERVBLK/DGPM,DGQM,DGVD,DBPM,DBQM,DBVD
c-----------------------------------------------------------------------
      DO  ISOT= 1,NISTP
          GFCT= 1.d0
          BFCT= 1.d0
          IF(ITYPE(ISTATE).EQ.2) GFCT= PNDE(0,ISTATE)
          IF(ITYPB(ISTATE).EQ.2) BFCT= PNDE(1,ISTATE)
          DO  IV= -1,VMAXX
              IF(IV.LT.0) THEN
                  VDMV= (VD(ISTATE)+ 0.5d0)*RSQMU(ISOT)
                ELSE
                  VDMV= (VD(ISTATE) - IV)*RSQMU(ISOT)
                ENDIF
              IF(NEWGv.GT.0) THEN
c** First get vib. energy numerator/exponent polynomial & its vD deriv.
                  SNUM= 1.d0
                  DNUM= 0.d0
                  SDEN= 1.d0     
                  DDEN= 0.d0
                  IF(NP0(ISTATE).GT.0) THEN
                      VDMVP= VDMV**IP0(ISTATE)
                      DO  I= 1,NP0(ISTATE)
                          DNUM= DNUM+(IP0(ISTATE)+I)*PM0(I,ISTATE)*VDMVP
                          VDMVP= VDMVP*VDMV
                          SNUM= SNUM+ PM0(I,ISTATE)*VDMVP
                          ENDDO
                      ENDIF
                  IF(NQ0(ISTATE).GT.0) THEN
c ... then get vib. energy denominator polynomial & its vD derivative
                      VDMVP= VDMV**IQ0(ISTATE)
                      DO  I= 1,NQ0(ISTATE)        
                          DDEN= DDEN+(IQ0(ISTATE)+I)*QM0(I,ISTATE)*VDMVP
                          VDMVP= VDMVP*VDMV      
                          SDEN= SDEN+ QM0(I,ISTATE)*VDMVP
                          ENDDO
                      ENDIF
c ... and now, store values & derivative components for use in DYIDPJ
                  EB= XM(0,ISTATE,ISOT)*VDMV**PNDE(0,ISTATE)
                  IF(ITYPE(ISTATE).EQ.3) THEN
                      EB= EB*DEXP(SNUM- 1.d0)
                      DGPM(IV,ISTATE,ISOT)= -EB
                      DGVD(IV,ISTATE,ISOT)= -EB*(PNDE(0,ISTATE)/VDMV
     1                                                         + DNUM)
                    ELSE
                      EB= EB*(SNUM/SDEN)**GFCT
                      DGPM(IV,ISTATE,ISOT)= -EB*GFCT/SNUM
                      DGQM(IV,ISTATE,ISOT)= EB*GFCT/SDEN
                      DGVD(IV,ISTATE,ISOT)= -EB*(PNDE(0,ISTATE)/VDMV
     1                               + GFCT*DNUM/SNUM- GFCT*DDEN/SDEN)
                    ENDIF
                  ZK(0,IV,ISTATE,ISOT)= DLIMIT(ISTATE)- EB 
                  ENDIF
c====== End of vibrational energy/derivative calculations =============
c
              IF((NDEBv(ISTATE).GT.0).AND.(NEWBv.GT.0)) THEN
c*** Now ... derivatives of Bv w.r.t. expansion coefficients & vD
c   First get Rotational NDE numerator/exponent polynomial & its vD deriv.
                  SNUM= 1.d0
                  DNUM= 0.d0
                  SDEN= 1.d0
                  DDEN= 0.d0
                  IF(NP1(ISTATE).GT.0) THEN
                      VDMVP= VDMV**IP1(ISTATE)
                      DO  I= 1,NP1(ISTATE)
                          DNUM= DNUM+(IP1(ISTATE)+I)*PM1(I,ISTATE)*VDMVP
                          VDMVP= VDMVP*VDMV
                          SNUM= SNUM+ PM1(I,ISTATE)*VDMVP
                          ENDDO
                      ENDIF
                  IF(NQ1(ISTATE).GT.0) THEN
c ... then get Rotational NDE denominator polynomial & its vD derivative
                      VDMVP= VDMV**IQ1(ISTATE)
                      DO  I= 1,NQ1(ISTATE)
                          DDEN= DDEN+(IQ1(ISTATE)+I)*QM1(I,ISTATE)*VDMVP
                          VDMVP= VDMVP*VDMV
                          SDEN= SDEN+ QM1(I,ISTATE)*VDMVP
                          ENDDO
                      ENDIF
c ... and now, store values & derivative components for use in DYIDPJ
                  BV= XM(1,ISTATE,ISOT)*VDMV**PNDE(1,ISTATE)
                  IF(ITYPB(ISTATE).EQ.3) THEN
                      BV= BV*DEXP(SNUM- 1.d0)
                      DBPM(IV,ISTATE,ISOT)= BV
                      DBVD(IV,ISTATE,ISOT)=BV*(PNDE(1,ISTATE)/VDMV+DNUM)
                    ELSE
                      BV= BV*(SNUM/SDEN)**BFCT
                      DBPM(IV,ISTATE,ISOT)= BV*BFCT/SNUM
                      DBQM(IV,ISTATE,ISOT)= -BV*BFCT/SDEN
                      DBVD(IV,ISTATE,ISOT)= BV*(PNDE(1,ISTATE)/VDMV
     1                                  + BFCT*(DNUM/SNUM- DDEN/SDEN))
                    ENDIF
                  ZK(1,IV,ISTATE,ISOT)= BV
                  IF(NEWBv.GT.0) ZK(1,IV,ISTATE,ISOT)= BV
                  ENDIF
              ENDDO
          ENDDO
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NDEDUN(ISTATE,NISTP,ZMASS,RSQMU,PU,CM)
c** Subroutine to calculate conventional Dunham parameters we, wexe, Be
c  and alpha_e (=AE) and their derivatives w.r.t. NDE parameters, from
c  input NDE-expansion functions.
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8  VDPH,VDPHI,VDPHP,FCT,FCT1,FCT2,FCT3,
     1  PWG,PWB,PU(NPARMX),CM(NPARMX,NPARMX),ZMASS(3,NISTPMX),
     2  RSQMU(NISTPMX),PFCT,SNUM,DNUM,D2NUM,D3NUM,SDEN,DDEN,D2DEN,
     3  D3DEN,EB,S1,S2,S3,S4
      REAL*8  WE(NISTPMX),XE(NISTPMX),BE(NISTPMX),AE(NISTPMX),
     1  RE(NISTPMX),UWE(NISTPMX),UXE(NISTPMX),UBE(NISTPMX),
     2  UAE(NISTPMX),DWE(2*NDUNMX,NISTPMX),DXE(2*NDUNMX,NISTPMX),
     3  URE(NISTPMX),DBE(2*NDUNMX,NISTPMX),DAE(2*NDUNMX,NISTPMX)
      INTEGER I,I1,ISOT,ISTATE,ICMS,ICMF,J,J1,NISTP,IP0I,IQ0I
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
c
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c-----------------------------------------------------------------------
c** First ... zero uncertainties and partial derivatives
      DO  ISOT= 1,NISTP
          UWE(ISOT)= 0.d0
          UXE(ISOT)= 0.d0
          UBE(ISOT)= 0.d0
          UAE(ISOT)= 0.d0
          DO  I= 1,2*NDUNMX
              DWE(I,ISOT)= 0.d0
              DXE(I,ISOT)= 0.d0
              DBE(I,ISOT)= 0.d0
              DAE(I,ISOT)= 0.d0
              ENDDO
          ENDDO
      PWG= PNDE(0,ISTATE)
      PWB= PNDE(1,ISTATE)
      IP0I= IP0(ISTATE)
      IQ0I= IQ0(ISTATE)
      ICMS= IPSTATE(ISTATE)+ 1
      IF(IFXD(ISTATE).LE.0) ICMS= ICMS+ 1
      ICMF= ICMS- 1
      IF(IFXVD(ISTATE).LE.0) ICMF= ICMF+1
      DO  ISOT= 1,NISTP
c** Loop over all isotopomers
          PFCT= 1.d0
          IF(ITYPE(ISTATE).EQ.2) PFCT= PWG
          VDPH= (VD(ISTATE)+0.5D0)*RSQMU(ISOT)
          VDPHI= RSQMU(ISOT)/VDPH
c** First get vib. energy numerator/exponent polynomial & its first 
c  three derivatives w.r.t. vD
          SNUM= 1.d0
          DNUM= 0.d0
          D2NUM= 0.d0
          D3NUM= 0.d0
          SDEN= 1.d0     
          DDEN= 0.d0
          D2DEN= 0.d0
          D3DEN= 0.d0
          IF(NP0(ISTATE).GT.0) THEN
              VDPHP= VDPH**IP0I
              DO  I= 1,NP0(ISTATE)
                  ICMF= ICMF+ 1
                  VDPHP= VDPHP*VDPH
                  FCT= PM0(I,ISTATE)*VDPHP
                  SNUM= SNUM+ FCT
                  FCT=(IP0I+ I)*FCT*VDPHI
                  DNUM= DNUM+ FCT
                  FCT= (IP0I+ I- 1)*FCT*VDPHI
                  D2NUM= D2NUM+ FCT
                  D3NUM= D2NUM+ (IP0I+ I- 2)*FCT*VDPHI
                  ENDDO
              ENDIF
          IF(NQ0(ISTATE).GT.0) THEN
c ... then get vib. energy denominator polynomial & its first three
c  derivatives w.r.t. vD
              VDPHP= VDPH**IQ0(ISTATE)
              DO  I= 1,NQ0(ISTATE)        
                  ICMF= ICMF+ 1
                  VDPHP= VDPHP*VDPH      
                  FCT= QM0(I,ISTATE)*VDPHP
                  SDEN= SDEN+ FCT
                  FCT=(IQ0I+ I)*FCT*VDPHI
                  DDEN= DDEN+ FCT
                  FCT= (IQ0I+ I- 1)*FCT*VDPHI
                  D2DEN= D2DEN+ FCT
                  D3DEN= D3DEN+ (IQ0I+ I- 2)*FCT*VDPHI
                  ENDDO
              ENDIF
c** Now generate  we, wexe & their deriv's w.r.t. vD & the p_i's
c  [derivative (0) w.r.t. vD;  (i) w.r.t. p_i ]
          EB= XM(0,ISTATE,ISOT)*VDPH**PWG
          J= 0
          IF(ITYPE(ISTATE).EQ.3) THEN
c** First, for exponential vibrational NDE
              EB= EB*DEXP(SNUM- 1.d0)
              WE(ISOT)= EB*(PWG*VDPHI + DNUM)
              XE(ISOT)= 0.5d0*EB*(PWG*(PWG-1.d0)*VDPHI**2 + 
     1                             DNUM*(2.d0*PWG*VDPHI+ DNUM) +D2NUM)
c ... first ... derivative w.r.t. vD
              IF(IFXVD(ISTATE).LE.0) THEN
                  J= J+1
                  DWE(J,ISOT)= XE(ISOT)
                  DXE(J,ISOT)= 0.5d0*EB*(PWG*VDPHI*(3.d0*(DNUM**2+ 
     1         D2NUM) + (PWG-1.d0)*VDPHI*(3.d0*DNUM+(PWG-2.d0)*VDPHI))
     2                              + DNUM**3+ 3.d0*DNUM*D2NUM+ D3NUM)
                  ENDIF
              IF(NP0(ISTATE).GT.0) THEN
c ... and then derivatives w.r.t. p_i's
                  FCT= VDPH**IP0I
                  DO  I= 1,NP0(ISTATE)
                      J= J+1
                      FCT= FCT*VDPH
                      DWE(J,ISOT)= EB*FCT*((PWG+ I+ IP0I)*VDPHI + DNUM)
                      DXE(J,ISOT)= FCT*(XE(ISOT)+ 0.5d0*EB*(I+IP0I)*
     1                  VDPHI*(VDPHI*(2.d0*PWG+i+IP0I-1) + 2.d0*DNUM))
                      ENDDO
                  ENDIF
            ELSE
c** Now for rational polynomial vibrational NDE's: ITYPE= 1 & 2
              EB= EB*(SNUM/SDEN)**PFCT
              FCT= PWG*VDPHI + PFCT*(DNUM/SNUM - DDEN/SDEN)
              FCT1= -PWG*VDPHI**2 + PFCT*(D2NUM/SNUM - (DNUM/SNUM)**2
     1                                  - D2DEN/SDEN + (DDEN/SDEN)**2)
              WE(ISOT)= EB*FCT
              XE(ISOT)= 0.5d0*EB*(FCT**2 + FCT1)
c ... now ... derivative w.r.t. vD
              IF(IFXVD(ISTATE).LE.0) THEN
                  J= J+1
                  DWE(J,ISOT)= 2.d0*XE(ISOT)
                  DXE(J,ISOT)= 0.5d0*WE(ISOT)*(FCT**2 + FCT1) + 
     1                     0.5d0*EB*(2.d0*FCT*FCT1 + 2.d0*PWG*VDPHI**3
     2                                  + PFCT*(D3NUM/SNUM -D3DEN/SDEN
     3                - 3.d0*(D2NUM*DNUM/SNUM**2 - D2DEN*DDEN/SDEN**2)
     4                      + 2.d0*((DNUM/SNUM)**3 - (DDEN/SDEN)**3)))
                  ENDIF
              IF(NP0(ISTATE).GT.0) THEN
c ... and then derivatives w.r.t. p_i's
                  FCT2= PFCT*EB*VDPH**IP0I/SNUM
                  DO  I= 1, NP0(ISTATE)
                      J= J+1
                      FCT2= FCT2*VDPH
                      FCT3= (I+IP0I)*VDPHI - DNUM/SNUM
                      DWE(J,ISOT)= FCT2*(FCT + FCT3)
                      DXE(J,ISOT)= 0.5d0*FCT2*(FCT**2 + FCT1+ 
     1       2.d0*FCT*FCT3 + (I+IP0I)*(I+IP0I-1)*VDPHI**2 - D2NUM/SNUM
     2                                          - 2.d0*FCT3*DNUM/SNUM)
                      ENDDO
                  ENDIF
              IF(NQ0(ISTATE).GT.0) THEN
c ... and then derivativesw.r.t. q_i's
                  FCT2= PFCT*EB*VDPH**IQ0I/SDEN
                  DO  I= 1, NQ0(ISTATE)
                      J= J+1
                      FCT2= FCT2*VDPH
                      FCT3= (I+IQ0I)*VDPHI - DDEN/SDEN
                      DWE(J,ISOT)= -FCT2*(FCT + FCT3)
                      DXE(J,ISOT)= -0.5d0*FCT2*(FCT**2+ FCT1+ 
     1       2.d0*FCT*FCT3 + (I+IQ0I)*(I+IQ0I-1)*VDPHI**2 - D2DEN/SDEN
     2                                          - 2.d0*FCT3*DDEN/SDEN)
                      ENDDO
                  ENDIF
     1                 
            ENDIF
c
c** Now get Bv exponent polynomial & its derivative w.r.t. vD
          IF(NDEBv(ISTATE).GT.0) THEN
              VDPHP= 1.d0
              SNUM= 0.d0
              DNUM= 0.d0
              D2NUM= 0.d0
              IF(NP1(ISTATE).GT.0) THEN
                  DO  I= 1,NP1(ISTATE)
                      ICMF= ICMF+ 1
                      VDPHP= VDPHP*VDPH                 
                      FCT= PM1(I,ISTATE)*VDPHP
                      SNUM= SNUM+ FCT
                      FCT= FCT*I*VDPHI
                      DNUM= DNUM+ FCT
                      D2NUM= D2NUM+ (I-1)*FCT*VDPHI
                      ENDDO
                  ENDIF
              IF(NQ1(ISTATE).GT.0) THEN
                  DO  I= 1,NQ1(ISTATE)
c???????? unfinished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                      ENDDO
                  ENDIF
              BE(ISOT)= XM(1,ISTATE,ISOT)*VDPH**PWB*DEXP(SNUM)
              AE(ISOT)= BE(ISOT)*(PWB*VDPHI+ DNUM)
c ... first, derivatives w.r.t. vD
              IF(IFXVD(ISTATE).LE.0) THEN
                  DBE(1,ISOT)= AE(ISOT)
                  DAE(1,ISOT)= AE(ISOT)*(PWB*VDPHI+ DNUM) 
     1                                + BE(ISOT)*(D2NUM- PWB*VDPHI**2)
                  ENDIF
              IF(NP1(ISTATE).GT.0) THEN
                  FCT= 1.d0
                  DO  I= 1,NP1(ISTATE)
                      J= J+1
                      FCT= FCT*VDPH
                      DBE(J,ISOT)= BE(ISOT)*FCT
                      DAE(J,ISOT)= DBE(J,ISOT)*(VDPHI*(PWB+I) + DNUM)
                      ENDDO
                  ENDIF
              IF(NQ1(ISTATE).GT.0) THEN
c???????? unfinished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                  ENDIF
              ENDIF
c** Finally ... use correlation matrix and uncertainties in fitted 
c  parameter to calculate uncertainties in our Dunham constants.
          I1= 0
          DO  I= ICMS,ICMF
              I1= I1+1
              J1= 0
              S1= 0.d0
              S2= 0.d0
              S3= 0.d0
              S4= 0.d0
              DO  J= ICMS,ICMF
                  J1= J1+1
                  S1= S1+ CM(I,J)*PU(J)*DWE(J1,ISOT)
                  S2= S2+ CM(I,J)*PU(J)*DXE(J1,ISOT)
                  IF(NDEBv(ISTATE).GT.0) THEN
                      S3= S3+ CM(I,J)*PU(J)*DBE(J1,ISOT)
                      S4= S4+ CM(I,J)*PU(J)*DAE(J1,ISOT)
                      ENDIF
                  ENDDO
              UWE(ISOT)= UWE(ISOT)+ PU(I)*DWE(I1,ISOT)*S1
              UXE(ISOT)= UXE(ISOT)+ PU(I)*DXE(I1,ISOT)*S2
              IF(NDEBv(ISTATE).GT.0) THEN
                  UBE(ISOT)= UBE(ISOT)+ PU(I)*DBE(I1,ISOT)*S3
                  UAE(ISOT)= UAE(ISOT)+ PU(I)*DAE(I1,ISOT)*S4
                  ENDIF
              ENDDO
          UWE(ISOT)= DSQRT(UWE(ISOT))
          UXE(ISOT)= DSQRT(UXE(ISOT))
          IF(NDEBv(ISTATE).GT.0) THEN
              UBE(ISOT)= DSQRT(UBE(ISOT))
              UAE(ISOT)= DSQRT(UAE(ISOT))
              RE(ISOT)= DSQRT(16.85762908D0/(BE(ISOT)*ZMASS(3,ISOT)))
              URE(ISOT)= RE(ISOT)*0.5d0*UBE(ISOT)/BE(ISOT)
              ENDIF
          ENDDO
      WRITE(6,600) ISTATE
  600 FORMAT(/' State-',i1,' isotopic Dunham parameters generated from N
     1DE functions:'/1x,64('-'))
      WRITE(6,602) (WE(ISOT),UWE(ISOT),ISOT= 1,NISTP)
      WRITE(6,604) (XE(ISOT),UXE(ISOT),ISOT= 1,NISTP)
      WRITE(6,606) (ZK(0,0,ISTATE,ISOT),ISOT= 1,NISTP)
      WRITE(6,608) (ZK(0,-1,ISTATE,ISOT),ISOT= 1,NISTP)
      WRITE(6,610) (BE(ISOT),UBE(ISOT),ISOT= 1,NISTP)
      WRITE(6,612) (AE(ISOT),UAE(ISOT),ISOT= 1,NISTP)
      WRITE(6,614) (RE(ISOT),URE(ISOT),ISOT= 1,NISTP)
  602 FORMAT('         we =',3(F11.5,' (',F9.5,')':)/
     1                                 (10x,3(F11.5,' (',F9.5,')':)))
  604 FORMAT('       wexe =',3(F11.5,' (',F9.5,')':)/
     1                                 (10x,3(F11.5,' (',F9.5,')':)))
  606 FORMAT('     T(v= 0)=',3(F14.5,9x:)/(10x,3(F14.5,9x:)))
  608 FORMAT('  T(v= -1/2)=',3(F14.5,9x:)/(10x,3(F14.5,9x:)))
  610 FORMAT('  B(v= -1/2)=',3(1PD12.5,' (',D8.1,')':)/
     1                                 (10x,3(D12.5,' (',D8.1,')':)))
  612 FORMAT('    alpha_e =',3(1PD12.5,' (',D8.1,')':)/
     1                                 (10x,3(D12.5,' (',D8.1,')':)))
  614 FORMAT(' Re(v= -1/2)=',3(F11.6,' (',F9.6,')':)/
     1                                 (10x,3(F11.6,' (',F9.6,')':)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPTOT,YC,PV,PD,PS)
c** Partial derivative subroutine called by general least-squares fitting
c  subroutine NLLSSRR for use in parameter-fits to diatomic molecule
c  spectroscopic data by R.J. Le Roy's program  DParFit. 
c                     Version of  02 April 2016
c   0701/13    Removed IFXP & RMSR as unused input parameters
c-----------------------------------------------------------------------
c** Input parameters PS (& RMSR) are used in cases when determine partial
c  derivatives by first differences (e.g., BCONT);  DUMMY variables here.
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      INTEGER  IVPP,IVP,M,NDATA,NPTOT,IDAT,IPAR,IPX,I,IBB,ESP,ESPP,ISOT,
     1  ISTATE,IPARVD,MQ0,MQM,ATOM,ATOM2,L,LAMIN
      REAL*8 ZATOM,JJPW,JJPPW,DER,PDVD,VDMVP,VDMVPP,VDMVPW,VDMVPPW,
     1  JJP,JJPI,JJPQ,JJPP,JJPPI,JJPPQ,YC,PV(NPTOT),PD(NPTOT),PS(NPTOT),
     2  SwP,SwPLR,SwPP,SwPPLR,dSwPVS,dSwPdVS,dSwPPVS,dSwPPdVS,YYDun
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c
c** Type statement and common block for NDE partial derivative stuff
c
      REAL*8 DGPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     1  DGQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     2  DGVD(-1:NVIBMX,NSTATEMX,NISTPMX),
     3  DBPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     4  DBQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     5  DBVD(-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /DERVBLK/DGPM,DGQM,DGVD,DBPM,DBQM,DBVD
c
c** For the first datum, call subroutine MAPPAR to map current NLLSSRR
c  values of fitted parameters PV onto "internal" molecular parameters
      IF(IDAT.EQ.1) CALL MAPPAR(NPTOT,PV)
c** Call subroutine to return current predicted value of datum-IDAT
      CALL PREDICT(IDAT,YC,NCDC,NLDMX,IOMEG,efREF,DLIMIT,ORIGIN,PV,
     1                                                          ZK,ZQ)
c
      IBB= IB(IDAT)
      ISOT= ISTP(IBB)
      ESP= IEP(IBB)
      ESPP= IEPP(IBB)
      IVP= VP(IBB)
      IVPP= VPP(IBB)
      IF(ESP.LE.0) IVPP= JP(IDAT)
      JJPP= JPP(IDAT)*(JPP(IDAT)+ 1.d0)
      JJPPQ= JJPP
      IF(IOMEG(ESPP).GT.0) JJPP= JJPP - IOMEG(ESPP)**2
      JJPPI= JJPP
c** Lower level isotope scaling ... as required ...
      IF((ISOT.GT.1).AND.(NDEGv(ESPP).GE.0)) JJPPI= JJPP*RMUP(1,ISOT)
      IF(ESP.GT.0) THEN
          JJP= JP(IDAT)*(JP(IDAT)+ 1.d0)
          JJPQ= JJP
          IF(IOMEG(ESP).GT.0) JJP= JJP - IOMEG(ESP)**2
          JJPI= JJP
c** Upper level isotope scaling ... as required ...
          IF((ISOT.GT.1).AND.(NDEGv(ESP).GE.0)) JJPI= JJP*RMUP(1,ISOT)
          ENDIF
c** Zero all partial derivatives ...
      DO  I= 1,NPTOT
          PD(I)= 0.d0
          ENDDO
      IF((ESP.EQ.-3).AND.(NDEBv(ESPP).EQ.-1).AND.(IFXBv(ESPP).LE.0))THEN
c** If Datum is a Bv value for a state being fitted by band constants ..
          IPX= NPAR(IVPP,ESPP,ISOT)
          IF(NRC(IVPP,ESPP,ISOT).GT.0) THEN
              IF(FITGV(IVPP,ESPP,ISOT).GT.0) IPX= IPX+ 1
              PD(IPX+ 1)= 1.d0
              ENDIF
          RETURN
          ENDIF
c
      DO 90 ISTATE= 1,NSTATES
c** Begin loop over states, accumulating parameter count & partial derivs.
          IF((ISTATE.NE.ESP).AND.(ISTATE.NE.ESPP)) GO TO 90
          IPAR= IPSTATE(ISTATE)
          IF((NDEGv(ISTATE).EQ.2).OR.(NDEBv(ISTATE).EQ.2)) THEN
c* If use MXS for Gv or Bv, generate switching functions & partial der.
              IF(ISTATE.EQ.ESP) THEN
                  SwPLR= dexp((IVP- VSISO(ESP,ISOT))/DVSISO(ESP,ISOT))
                  SwP= 1.d0/(1.d0+ SwPLR)
                  SwPLR= SwPLR*SwP
                  IF(IFXVS(ISTATE).LE.0) THEN
                      dSwPVS= SwPLR/DVS(ISTATE)
                      dSwPdVS= dSwPVS*(IVP- VSISO(ESP,ISOT))/
     1                                                DVSISO(ESP,ISOT)
                      ENDIF
                  ENDIF
              IF(ISTATE.EQ.ESPP) THEN
                  SwPPLR= dexp((IVPP- VSISO(ESPP,ISOT))/
     1                                              DVSISO(ESPP,ISOT))
                  SwPP= 1.d0/(1.d0+ SwPPLR)
                  SwPPLR= SwPPLR*SwPP
                  IF(IFXVS(ISTATE).LE.0) THEN
                      dSwPPVS= SwPPLR/DVS(ISTATE)
                      dSwPPdVS= dSwPPVS*(IVP- VSISO(ESP,ISOT))/
     1                                                DVSISO(ESP,ISOT)
                      ENDIF
                  ENDIF
              ENDIF
c
c** If fitting to parameters defining Gv's for this state ..............
c=======================================================================
          IF(IFXGv(ISTATE).LE.0) THEN
              IF(NDEGv(ISTATE).EQ.-2) THEN
c** If fitting to individual term values for this state ...
c==========================================================
                  IF((ISTATE.EQ.ESP).AND.(TVUP(IDAT).GT.0))
     1                                            PD(TVUP(IDAT))= 1.d0
                  IF((ISTATE.EQ.ESPP).AND.(TVLW(IDAT).GT.0))
     1                                           PD(TVLW(IDAT))= -1.d0
                  GOTO 90
                  ENDIF
              IF(NDEGv(ISTATE).EQ.-1) THEN
c** If fitting Gv (& rotational constants) as band constants ...
c================================================================
                  IF(ISTATE.EQ.ESPP) THEN
c ... derivatives for ISTATE being the lower state
                      IPX= NPAR(IVPP,ESPP,ISOT)
                      IF(FITGV(IVPP,ESPP,ISOT).GT.0) THEN
                          IPX= IPX+ 1
                          PD(IPX)= PD(IPX) - 1.d0
                          ENDIF
                      IF(NRC(IVPP,ESPP,ISOT).GT.0) THEN
                          JJPPW= 1.d0
                          DO  M= 1,NRC(IVPP,ESPP,ISOT)
                              IPX= IPX+ 1
                              JJPPW= JJPPW*JJPP
                              PD(IPX)= PD(IPX) - JJPPW
                              ENDDO
                          ENDIF
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
c ... then ... derivatives for ISTATE being the upper state
                      IPX= NPAR(IVP,ESP,ISOT)
                      IF(FITGV(IVP,ESP,ISOT).GT.0) THEN
                          IPX= IPX+ 1
                          PD(IPX)= PD(IPX)+ 1.d0
                          ENDIF
                      IF(NRC(IVP,ESP,ISOT).GT.0) THEN
                          JJPW= 1.d0
                          DO  M= 1,NRC(IVP,ESP,ISOT)
                              IPX= IPX+ 1
                              JJPW= JJPW*JJP
                              PD(IPX)= PD(IPX)+ JJPW
                              ENDDO
                          ENDIF
                      ENDIF
                  IPAR= IPX
                  GO TO 60
                  ENDIF
c
c** If Gv's for this state fitted to Dunham or mixed (MXS) function .. 
c=======================================================================
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).EQ.2)) THEN
                  IF(((ESP.EQ.ESPP).AND.(IVP.EQ.IVPP)).OR.
     1                                               (ESP.EQ.-3)) THEN
c** For MW data or "experimental" Bv's ... skip vib. derivatives
                      IF(ISTATE.GT.1) IPAR= IPAR+1
                      IF(LMAX(0,ISTATE).GT.0) IPAR= IPAR+ LMAX(0,ISTATE)
                      GO TO 30
                      ENDIF
                  IF(ISTATE.GT.1) THEN
c** First ... derivative w.r.t. Te  for this state 
                      IPAR= IPAR+ 1
                      IF(ESP.NE.ESPP) THEN
                          IF(ISTATE.EQ.ESP) THEN
                              IF(NDEGv(ISTATE).EQ.0) PD(IPAR)= 1.d0
                              IF(NDEGv(ISTATE).GE.2) PD(IPAR)= SwP
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              IF(NDEGv(ISTATE).EQ.0) PD(IPAR)=-1.d0
                              IF(NDEGv(ISTATE).GE.2) PD(IPAR)=-SwPP
                              ENDIF
                          ENDIF
                      ENDIF
c ... Now ... derivatives w.r.t. Dunham vibrational coefficients
                  IF(LMAX(0,ISTATE).GT.0) THEN
                      YYDun= YLM(0,0,ISTATE)
                      DO  L= 1,LMAX(0,ISTATE)
                          IPAR= IPAR+ 1
                          IF(ISTATE.EQ.ESP) THEN
                              PD(IPAR)= VPHPW(IVP,L)*RSQMUP(L,ISOT)
                              IF(NDEGv(ISTATE).GE.2) THEN
                                  PD(IPAR)= PD(IPAR)*SwP
                                  YYDun= YYDun+ PD(IPAR)*YLM(L,0,ISTATE)
                                  ENDIF
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              IF(NDEGv(ISTATE).EQ.0) PD(IPAR)= PD(IPAR)-
     1                                    VPHPW(IVPP,L)*RSQMUP(L,ISOT)
                              IF(NDEGv(ISTATE).GE.2) PD(IPAR)= PD(IPAR)-
     1                               VPHPW(IVPP,L)*RSQMUP(L,ISOT)*SwPP
                              ENDIF
                          ENDDO
                      ENDIF
                  ENDIF
c
c** If Gv's for this state fitted to NDE or MXS functions ... 
c=============================================================
              IF(NDEGv(ISTATE).GE.1) THEN
                  IF(((ESP.EQ.ESPP).AND.(IVP.EQ.IVPP)).OR.
     1                                               (ESP.EQ.-3)) THEN
c** For MW data or "experimental" Bv's ... skip vib. derivatives
                      IF(IFXD(ISTATE).LE.0) IPAR= IPAR+ 1
                      IF(IFXVD(ISTATE).LE.0) THEN
                          IPAR= IPAR+ 1
                          IPARVD= IPAR
                          PDVD= 0.d0
                          IPAR= IPAR+ NP0(ISTATE)+ NQ0(ISTATE)
                          ENDIF
                      GO TO 20
                      ENDIF
                  IF(IFXD(ISTATE).LE.0) THEN
c** First ... derivative w.r.t. DLIMIT for this state (if appropriate)
c** Note that DLIMIT's for multiple NDE-represented states are coupled
                      IPAR= IPAR+ 1
                      IF(ESP.NE.ESPP) THEN
                          IF(ISTATE.EQ.ESP) THEN
                              IF(NDEGv(ISTATE).EQ.1) PD(IPAR)= 1.d0
                              IF(NDEGv(ISTATE).GE.2) PD(IPAR)= SwPLR
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              IF(NDEGv(ISTATE).EQ.1) PD(IPAR)= -1.d0
                              IF(NDEGv(ISTATE).GE.2) PD(IPAR)= -SwPPLR
c ... for PAS data - no dependence on DLIMIT ...
                              IF(ESP.EQ.-1) PD(IPAR)= 0.d0
                              ENDIF
                          ENDIF
                      ENDIF
                  IF(IFXVD(ISTATE).LE.0) THEN
c** Now ... vibrational derivative w.r.t. vD
                      IPAR= IPAR+ 1
                      IPARVD= IPAR
                      PDVD= 0.d0
                      IF(ISTATE.EQ.ESP) THEN
                          PDVD= DGVD(IVP,ISTATE,ISOT)
                          IF(NDEGv(ISTATE).GE.2) PDVD= PDVD*SwPLR
                          ENDIF
                      IF(ISTATE.EQ.ESPP) THEN
                          IF(NDEGv(ESPP).EQ.1) 
     1                             PDVD= PDVD - DGVD(IVPP,ISTATE,ISOT)
                          IF(NDEGv(ESPP).GE.2) 
     1                      PDVD= PDVD - DGVD(IVPP,ISTATE,ISOT)*SwPPLR
                          ENDIF
                      PD(IPAR)= PDVD
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
c ... prepare for numerator/denominator term sums ...
                      VDMVP= (VD(ISTATE)-IVP)*RSQMU(ISOT)
                      VDMVPW= VDMVP**IP0(ISTATE)
                      ENDIF
                  IF(ISTATE.EQ.ESPP) THEN
                      VDMVPP= (VD(ISTATE)-IVPP)*RSQMU(ISOT)
                      VDMVPPW= VDMVPP**IP0(ISTATE)
                      ENDIF
                  IF(NP0(ISTATE).GT.0) THEN
c ... then w.r.t. NDE Gv numerator polynomial coefficients
                      DO  I= 1,NP0(ISTATE)
                          IPAR= IPAR+ 1
                          DER= 0.d0
                          IF(ISTATE.EQ.ESP) THEN
                              VDMVPW= VDMVPW* VDMVP
                              DER= DGPM(IVP,ISTATE,ISOT)*VDMVPW
                              IF(NDEGv(ESP).GE.2) DER= DER*SwPLR
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              VDMVPPW= VDMVPPW* VDMVPP 
                              IF(NDEGv(ESPP).EQ.1) DER= DER - 
     1                                  DGPM(IVPP,ISTATE,ISOT)*VDMVPPW
                              IF(NDEGv(ESPP).GE.2) DER= DER - 
     1                           DGPM(IVPP,ISTATE,ISOT)*VDMVPPW*SwPPLR
                              ENDIF
                          PD(IPAR)= DER
                          ENDDO
                      ENDIF
                  IF(NQ0(ISTATE).GT.0) THEN
c ... then w.r.t. NDE Gv denominator polynomial coefficients
                      IF(ISTATE.EQ.ESP) VDMVPW= VDMVP**IQ0(ISTATE)
                      IF(ISTATE.EQ.ESPP) VDMVPPW= VDMVPP**IQ0(ISTATE)
                      DO  I= 1,NQ0(ISTATE)
                          IPAR= IPAR+ 1
                          DER= 0.d0
                          IF(ISTATE.EQ.ESP) THEN
                              VDMVPW= VDMVPW* VDMVP
                              DER= DGQM(IVP,ISTATE,ISOT)*VDMVPW
                              IF(NDEGv(ESP).GE.2) DER= DER*SwPPLR
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              VDMVPPW= VDMVPPW* VDMVPP
                              IF(NDEGv(ESPP).EQ.1) DER= DER-
     1                                  DGQM(IVPP,ISTATE,ISOT)*VDMVPPW
                              IF(NDEGv(ESPP).GE.2) DER= DER-
     1                           DGQM(IVPP,ISTATE,ISOT)*VDMVPPW*SwPPLR
                              ENDIF
                          PD(IPAR)= DER
                          ENDDO
                      ENDIF
c ... now include Bv contribution in the derivative w.r.t.  vD 
   20             IF((IFXVD(ISTATE).LE.0).AND.(NDEBv(ISTATE).GT.0)) THEN
                      IF(ISTATE.EQ.ESP) THEN
                          IF(NDEBv(ESP).EQ.1) 
     1                             PDVD= PDVD+ DBVD(IVP,ESP,ISOT)*JJPI
                          IF(NDEBv(ESP).GE.2) 
     1                       PDVD= PDVD+ DBVD(IVP,ESP,ISOT)*JJPI*SwPLR
                          ENDIF
                      IF(ISTATE.EQ.ESPP) THEN
                          IF(NDEBv(ESPP).EQ.1) PDVD= PDVD-
     1                                      DBVD(IVPP,ESPP,ISOT)*JJPPI
                          IF(NDEBv(ESPP).GE.2) PDVD= PDVD-
     1                               DBVD(IVPP,ESPP,ISOT)*JJPPI*SwPPLR
                          ENDIF
                      PD(IPARVD)= PDVD
                      ENDIF
                  ENDIF
c ... end of fit to vibrational NDE
              ENDIF
c ... end of fit to vibrational parameters
c=======================================================================
c** Begin treatment of fitted Bv parameters ...........................
c=======================================================================
   30     IF(IFXBv(ISTATE).LE.0) THEN
              IF(NDEBv(ISTATE).EQ.-1) THEN
c-----------------------------------------------------------------------
c** If fitting to band-constant rotational (and hence CDC) constants
c  while Gv's represented by a Dunham, NDE or MXS analytic function ...
c-----------------------------------------------------------------------
                  IF(ISTATE.EQ.ESPP) THEN
c ... derivatives for ISTATE being the lower state
                      IPX= NPAR(IVPP,ESPP,ISOT)
                      IF(NRC(IVPP,ESPP,ISOT).GT.0) THEN
                          JJPPW= 1.d0
                          DO  M= 1,NRC(IVPP,ESPP,ISOT)
                              IPX= IPX+ 1
                              JJPPW= JJPPW*JJPP
                              PD(IPX)= PD(IPX) - JJPPW
                              ENDDO
                          ENDIF
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
c ... then ... derivatives for ISTATE being the upper state
                      IPX= NPAR(IVP,ESP,ISOT)
                      IF(NRC(IVP,ESP,ISOT).GT.0) THEN
                          JJPW= 1.d0
                          DO  M= 1,NRC(IVP,ESP,ISOT)
                              IPX= IPX+ 1
                              JJPW= JJPW*JJP
                              PD(IPX)= PD(IPX)+ JJPW
                              ENDDO
                          ENDIF
                      ENDIF
c** Set counter in case wishing to fit to BOB corrn. for Gv ...
                  IPAR= NEBC(ISTATE)
                  GO TO 60
                  ENDIF
c** If fitting Bv's to a Dunham or NDE or MXS function ...
c---------------------------------------------------------
              IF((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).EQ.2)) THEN
c* If Bv's for this state fitted to a Dunham or MXS function ...
c---------------------------------------------------------------
                  DO  L= 0,LMAX(1,ISTATE)
                      IPAR= IPAR+ 1
                      DER= 0.d0
                      IF(ESP.EQ.-3) THEN
                          DER= VPHPW(IVPP,L)
                          IF(NDEBv(ISTATE).EQ.2) DER= DER*SwPP
                        ELSE
                          IF(ISTATE.EQ.ESP) THEN
                              DER= VPHPW(IVP,L)*JJPI
                              IF(NDEBv(ISTATE).GE.2) DER= DER*SwP
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              IF(NDEBv(ISTATE).EQ.0) DER= DER-
     1                                             VPHPW(IVPP,L)*JJPPI
                              IF(NDEBv(ISTATE).GE.2) DER= DER-
     1                                        VPHPW(IVPP,L)*JJPPI*SwPP
                              ENDIF
                        ENDIF
                      PD(IPAR)= DER*RSQMUP(L,ISOT)
                      ENDDO
                  IF((ESP.EQ.-3).AND.(NDEBv(ISTATE).EQ.0)) RETURN
                  ENDIF
              IF(NDEBv(ISTATE).GT.0) THEN
c
c** If Bv's for this state fitted to pure NDE or MXS functions ...
c-----------------------------------------------------------------
                  IF(IFXGv(ISTATE).GT.0) THEN
c ... if NDE Gv's not being fitted, define expansions variables here ...
                      IF(ISTATE.EQ.ESP)
     1                             VDMVP= (VD(ISTATE)-IVP)*RSQMU(ISOT)
                      IF(ISTATE.EQ.ESPP)
     1                           VDMVPP= (VD(ISTATE)-IVPP)*RSQMU(ISOT)
                      ENDIF
                  IF(NP1(ISTATE).GT.0) THEN
c ... for numerator/exponent polynomial coefficients ...
                      IF(ISTATE.EQ.ESP) VDMVPW= JJPI*VDMVP**IP1(ISTATE)
                      IF(ISTATE.EQ.ESPP) VDMVPPW=
     1                                       JJPPI*VDMVPP**IP1(ISTATE)
                      IF(ESP.EQ.-3) VDMVPPW= VDMVP**IP1(ISTATE)
                      DO  I= 1,NP1(ISTATE)
                          IPAR= IPAR+ 1
                          DER= 0.d0
                          IF(ISTATE.EQ.ESP) THEN
                              VDMVPW= VDMVPW* VDMVP
                              DER= DBPM(IVP,ISTATE,ISOT)*VDMVPW
                              IF(NDEBv(ISTATE).GE.2) DER= DER*SwPLR
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              VDMVPPW= VDMVPPW* VDMVPP
                              IF(NDEBv(ISTATE).EQ.1) DER= DER-
     1                                  DBPM(IVPP,ISTATE,ISOT)*VDMVPPW
                              IF(NDEBv(ISTATE).GE.2) DER= DER-
     1                           DBPM(IVPP,ISTATE,ISOT)*VDMVPPW*SwPPLR
                              ENDIF
                          PD(IPAR)= DER
                          ENDDO
                      ENDIF
                  IF(NQ1(ISTATE).GT.0) THEN
                      IF(ISTATE.EQ.ESP) VDMVPW= JJPI*VDMVP**IQ1(ISTATE)
                      IF(ISTATE.EQ.ESPP) VDMVPPW= 
     1                                       JJPPI*VDMVPP**IQ1(ISTATE)
                      IF(ESP.LT.-3) VDMVPPW= VDMVPP**IQ1(ISTATE)
                      DO  I= 1,NQ1(ISTATE)
                          IPAR= IPAR+ 1
                          DER= 0.d0
                          IF(ISTATE.EQ.ESP) THEN
                              VDMVPW= VDMVPW* VDMVP
                              DER= DBQM(IVP,ISTATE,ISOT)*VDMVPW
                              IF(NDEBv(ISTATE).GE.2) DER= DER*SwPLR
                              ENDIF
                          IF(ISTATE.EQ.ESPP) THEN
                              VDMVPPW= VDMVPPW* VDMVPP
                              IF(NDEBv(ISTATE).EQ.1) DER= DER-
     1                                  DBQM(IVPP,ISTATE,ISOT)*VDMVPPW
                              IF(NDEBv(ISTATE).GE.2) DER= DER-
     1                           DBQM(IVPP,ISTATE,ISOT)*VDMVPPW*SwPPLR
                              ENDIF
                          PD(IPAR)= DER
                          ENDDO
                      ENDIF
                  IF(ESP.EQ.-3) RETURN
                  ENDIF
              ENDIF
c
c=======================================================================
c** Begin treatment of fitted CDC parameters ...........................
c=======================================================================
          IF(IFXCDC(ISTATE).LE.0) THEN
              IF(NDECDC(ISTATE).EQ.-1) THEN
c----------------------------------------------------------------------
c** If fitting CDCs to band constants (while Gv & Bv are Dunham or ...)
c----------------------------------------------------------------------
                  IF(ISTATE.EQ.ESPP) THEN
c ... derivatives for ISTATE being the lower state
                      IPX= NPAR(IVPP,ESPP,ISOT)
                      IF(NRC(IVPP,ESPP,ISOT).GT.0) THEN
                          JJPPW= JJPP
                          DO  M= 2,NRC(IVPP,ESPP,ISOT)
                              IPX= IPX+ 1
                              JJPPW= JJPPW*JJPP
                              PD(IPX)= PD(IPX) - JJPPW
                              ENDDO
                          ENDIF
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
c ... then ... derivatives for ISTATE being the upper state
                      IPX= NPAR(IVP,ESP,ISOT)
                      IF(NRC(IVP,ESP,ISOT).GT.0) THEN
                          JJPW= JJP
                          DO  M= 2,NRC(IVP,ESP,ISOT)
                              IPX= IPX+ 1
                              JJPW= JJPW*JJP
                              PD(IPX)= PD(IPX)+ JJPW
                              ENDDO
                          ENDIF
                      ENDIF
c** Set counter in case wishing to fit to BOB corrn. for Gv or Bv ...
                  IPAR= NEBC(ISTATE)
                  GO TO 60
                  ENDIF
c
              IF(NDECDC(ISTATE).EQ.0) THEN
c-------------------------------------------------------
c** If fitting to Dunham expansions for the CDC's ......
c-------------------------------------------------------
                  JJPW= JJPI
                  JJPPW= JJPPI
                  DO 50 M= 2,NCDC(ISTATE)+1
                      IF(ESP.EQ.-3) THEN
                          IPAR= IPAR+ LMAX(M,ISTATE)+ 1
                          GO TO 50
                          ENDIF
                      JJPW= JJPW*JJPI
                      JJPPW= JJPPW*JJPPI
                      DO  L= 0,LMAX(M,ISTATE)
                          IPAR= IPAR+ 1
                          DER= 0.d0
                          IF(ISTATE.EQ.ESP) DER= DER+ VPHPW(IVP,L)
     1                                            *RSQMUP(L,ISOT)*JJPW
                          IF(ISTATE.EQ.ESPP) DER= DER- VPHPW(IVPP,L)
     1                                           *RSQMUP(L,ISOT)*JJPPW
                          PD(IPAR)= DER                               
                          ENDDO
   50                 CONTINUE 
                  ENDIF 
              ENDIF
c=======================================================================
c** Now ... derivatives w.r.t. any Lambda or ^2\Sigma doubling constants
c=======================================================================
   60     IF((IOMEG(ISTATE).NE.0).AND.(IFXLD(ISTATE).LE.0)
     1                                 .AND.(NLDMX(ISTATE).GE.1)) THEN
c ... if fitting to Band-Constant values for doubling constants ...
c------------------------------------------------------------------
              MQ0= MAX0(0,IOMEG(ISTATE)-1)
              IF(NDELD(ISTATE).LT.0) THEN
                  IF(ISTATE.EQ.ESPP) THEN
c ... derivatives for ISTATE being the lower state .....
                      IF(NQC(IVPP,ESPP,ISOT).GT.0) THEN
                          IF(IOMEG(ISTATE).LT.0) THEN
                              JJPPW= JPP(IDAT)-MIN(0,EFPP(IDAT))
                            ELSE
                              JJPPW= JJPPQ**IOMEG(ISTATE)
                            ENDIF
                          IPAR= NQPAR(IVPP,ESPP,ISOT)
                          DO  M= 1,NQC(IVPP,ESPP,ISOT)
                              IPAR= IPAR+ 1
c ... NOTE ... neglect derivative if  EFPP(idat)= 0
                              IF(EFPP(IDAT).NE.0) PD(IPAR)= PD(IPAR)  
     1                     - (0.5D0*(EFPP(IDAT)- efREF(ISTATE))*JJPPW)
                              JJPPW= JJPPW*JJPP
                              ENDDO
                          ENDIF
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
c ... derivatives for ISTATE being the upper state .....
                      IF(NQC(IVP,ISTATE,ISOT).GT.0) THEN
                          IF(IOMEG(ISTATE).LT.0) THEN
                              JJPW= JP(IDAT)-MIN(0,EFP(IDAT))
                            ELSE
                              JJPW= JJPQ**IOMEG(ISTATE)
                            ENDIF
                          IPAR= NQPAR(IVP,ISTATE,ISOT)
                          DO  M= 1,NQC(IVP,ESP,ISOT)
                              IPAR= IPAR+ 1
                              IF(EFP(IDAT).NE.0) PD(IPAR)= PD(IPAR)
     1                         + 0.5D0*(EFP(IDAT)- efREF(ISTATE))*JJPW
                              JJPW= JJPW*JJP
                              ENDDO
                          ENDIF
                      ENDIF
                  IPAR= NQPAR(VMAX(ISTATE)+1,ISTATE,NISTP)
                  ENDIF
c
c ... if using Dunham-type expansions for doubling constants ...
c------------------------------------------------------------------
              IF(NDELD(ISTATE).GE.0) THEN
c ... First, define rotational q.No. factors ...
                  IF(ISTATE.EQ.ESPP) THEN
                      IF(IOMEG(ISTATE).LT.0)
     1               JJPPW= (JPP(IDAT)-MIN(0,EFPP(IDAT)))*RMUP(1,ISOT)
                      IF(IOMEG(ISTATE).GT.0)
     1                    JJPPW= (JJPPQ * RMUP(2,ISOT))**IOMEG(ISTATE)
                      ENDIF
                  IF(ISTATE.EQ.ESP) THEN
                      IF(IOMEG(ISTATE).LT.0)
     1                  JJPW= (JP(IDAT)-MIN(0,EFP(IDAT)))*RMUP(1,ISOT)
                      IF(IOMEG(ISTATE).GT.0)
     1                      JJPW= (JJPQ * RMUP(2,ISOT))**IOMEG(ISTATE)
                      ENDIF
                  IPAR= NQPAR(0,ISTATE,1)
                  DO  M= 1,NLDMX(ISTATE)
                      MQM= MQ0+ M
                      IF(LDMAX(MQM,ISTATE).GE.0) THEN
                          DO  L= 0,LDMAX(MQM,ISTATE)
                              IPAR= IPAR+ 1
                              DER= 0.d0
c ... Note - neglect if  EFP(idat)= 0
                              IF((ISTATE.EQ.ESP).AND.(EFP(IDAT).NE.0))
     1                          DER= DER + VPHPW(IVP,L)*RSQMUP(L,ISOT)
     2                                *(EFP(IDAT)- efREF(ISTATE))*JJPW
                              IF((ISTATE.EQ.ESPP).AND.(EFPP(IDAT).NE.0))
     1                          DER= DER- VPHPW(IVPP,L)*RSQMUP(L,ISOT)
     2                              *(EFPP(IDAT)- efREF(ISTATE))*JJPPW
                              PD(IPAR)= 0.5d0*DER
                              ENDDO
                          ENDIF
                      JJPW= JJPW*JJPI
                      JJPPW= JJPPW*JJPPI
                      ENDDO
                  ENDIF
              ENDIF
c======================================================
c** Partial derivatives w.r.t. B-O-B  delta  parameters
c======================================================
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).LE.0)) THEN
c ... loop over atoms ... then over M ... & then over L
              ZATOM= 1.d0 - ZMASS(1,1)/ZMASS(1,ISOT)
              ATOM2= 2
              IF(AN(1).EQ.AN(2)) THEN
c ... and for a homonuclear species, add contributions for the two atoms
                  ATOM2= 1
                  ZATOM= ZATOM+ 1.d0 - ZMASS(2,1)/ZMASS(2,ISOT)
                  ENDIF 
              DO  ATOM= 1,ATOM2
                  JJPW= 1.d0
                  JJPPW= 1.d0
                  LAMIN= 0
                  IF((ISTATE.EQ.1).AND.(BOB00.LE.0)) LAMIN= 1
                  DO 80 M= 0,BOBORD(ISTATE)
                      IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                          IF(ESP.EQ.-3) THEN
                              IPAR= IPAR+ LAMAX(ATOM,M,ISTATE)+ 1-LAMIN
                              GO TO 80
                              ENDIF
                          DO  L= LAMIN,LAMAX(ATOM,M,ISTATE)
                              IPAR= IPAR+ 1
                              DER= 0.d0
                              IF(ISTATE.EQ.ESP) THEN
                                  DER= ZATOM*VPHPW(IVP,L)*
     1                                         RSQMUP(L,ISOT)*JJPW
                                  ENDIF
                              IF(ISTATE.EQ.ESPP) THEN
                                  DER= DER- ZATOM*VPHPW(IVPP,L)*
     1                                        RSQMUP(L,ISOT)*JJPPW
                                  ENDIF
                              PD(IPAR)= DER
                              ENDDO
                          ENDIF
                      JJPW= JJPW*JJPI
                      JJPPW= JJPPW*JJPPI
   80                 LAMIN= 0
                  ZATOM= 1.d0 - ZMASS(2,1)/ZMASS(2,ISOT)
                  ENDDO
              ENDIF
   90     CONTINUE
c
c** And finally ... derivatives w.r.t. to Fluorescence series origin
c===================================================================
      IF(ESP.EQ.0) PD(NPTOT- NFSTOT+ NFS(IBB))= 1.d0
c     if(idat.eq.300) write(9,901) nparm,(pd(i),i=1,nparm)
c 901 format(/i6/(4(1Pd20.11)))
      RETURN
c** End of Partial Derivatives for fit to Dunham and/or ND expansions
c=======================================================================
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MAPPAR(NPTOT,PV)
c** Subroutine to MAP PARameters PV(i) worked with inside NLLSSRR onto 
c  the NPTOT free parameters of the actual model.
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      INTEGER  ATOM,ATOM2,I,IV,IPARDLIM,IFS,IPAR,ISTATE,ISOT,
     1  L,LAMIN,M,MQ0,MQM,MMIN,MMAX,NPTOT, NEWGv,NEWBv
      REAL*8  PV(NPTOT),XX,XXP,YY, ZATOM, Sw, SwLR
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
c
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c
c** Type statement and common block for NDE partial derivative stuff
c
      REAL*8 DGPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     1  DGQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     2  DGVD(-1:NVIBMX,NSTATEMX,NISTPMX),
     3  DBPM(-1:NVIBMX,NSTATEMX,NISTPMX),
     4  DBQM(-1:NVIBMX,NSTATEMX,NISTPMX),
     5  DBVD(-1:NVIBMX,NSTATEMX,NISTPMX)
      COMMON /DERVBLK/DGPM,DGQM,DGVD,DBPM,DBQM,DBVD
c
      IPAR= 0
      IPARDLIM= -1
      NEWGv= 0
      NEWBv= 0
      DO 50 ISTATE= 1,NSTATES
c** Now - identify and update parameters from the fit
c
c** If use all Term Value representation for this state ...
c==========================================================
          IF((NDEGv(ISTATE).EQ.-2).AND.(IFXGv(ISTATE).LE.0)) THEN
c ... need to loop over isotopomers ...
              DO  ISOT= 1,NISTP
c ... and cumulatively count parameters ...
                  IPAR= IPAR+ NTV(ISTATE,ISOT)
                  ENDDO
              GOTO 50
              ENDIF
c
c** If use band-constants to represent vib/rot term values of this state
c=======================================================================
          IF((NDEGv(ISTATE).EQ.-1).AND.(IFXGv(ISTATE).LE.0)) THEN
c ... first do outer loop over isotopomers ...
              DO  ISOT= 1,NISTP
c ... and then inner loop over vibrational levels.
                  DO  IV= VMIN(ISTATE), VMAX(ISTATE)
                      IF(FITGV(IV,ISTATE,ISOT).GT.0) THEN
                          IPAR= IPAR+ 1
                          ZK(0,IV,ISTATE,ISOT)=  PV(IPAR)
                        ENDIF
                      IF(NRC(IV,ISTATE,ISOT).GT.0) THEN
                          DO  M= 1,NRC(IV,ISTATE,ISOT)
                              IPAR= IPAR+ 1
                              ZK(M,IV,ISTATE,ISOT)= PV(IPAR)
                              ENDDO
                          ENDIF
                      ENDDO
                  ENDDO
                  GO TO 30
              ENDIF
c====end of section for fit to term values or vib-rot band constants====
c
c*** First ... update Gv expansion parameters
c============================================
c*** If using Dunham or NDE or MXS Gv function for this state
          MMIN= -1
          IF(IFXGv(ISTATE).LE.0) THEN
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).GE.2)) THEN
c ... First ... for upper (ISTATE > 1) electronic state, update Te for
c     pure Dunham or MXS case
                  IF(ISTATE.GT.1) THEN
                      IPAR= IPAR+ 1
                      Te(ISTATE)= PV(IPAR)
                      YLM(0,0,ISTATE)= Te(ISTATE)
                      ENDIF
c ... next update Dunham Gv parameters for MXS or pure Dunham cases ...
                  IF(LMAX(0,ISTATE).GE.1) THEN
                      DO  L= 1, LMAX(0,ISTATE)
                          IPAR= IPAR+ 1
                          YLM(L,0,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  MMIN= 0
                  MMAX= 0
                  ENDIF
c
c*** If using MXS function for Gv - update VS & DVS if they are fitted
              IF((NDEGv(ISTATE).GE.2).AND.(IFXVS(ISTATE).LE.0)) THEN
                  IPAR= IPAR+ 1
                  VS(ISTATE)= PV(IPAR)
                  IPAR= IPAR+ 1
                  DVS(ISTATE)= PV(IPAR)
                  XX= VS(ISTATE)+ 0.5d0
                  DO  ISOT= 1, NISTP
                      VSISO(ISTATE,ISOT)= XX/RSQMU(ISOT)- 0.5d0
                      DVSISO(ISTATE,ISOT)= DVS(ISTATE)/RSQMU(ISOT)
                      ENDDO 
                  ENDIF
c
c*** If using NDE or MXS functions, update NDE parameters for Gv
              IF(NDEGv(ISTATE).GE.1) THEN
                  NEWGv= 1
                  IF(IFXD(ISTATE).LE.0) THEN
c ... First the value of DLIM (if it was floated)
                      IPAR= IPAR+1
                      DLIMIT(ISTATE)= PV(IPAR)
                      ENDIF
                  IF(IFXVD(ISTATE).LE.0) THEN
c ... then the value of  vD  for this state (if free)      
                      IPAR= IPAR+ 1
                      VD(ISTATE)= PV(IPAR)   
                      ENDIF
c ... then the vibrational numerator polynomial coefficients
                  IF(NP0(ISTATE).GT.0) THEN
                      DO  I= 1,NP0(ISTATE)
                          IPAR= IPAR+ 1
                          PM0(I,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  IF(NQ0(ISTATE).GT.0) THEN   
c ... then the vibrational denominator polynomial coefficients
                      DO  I= 1,NQ0(ISTATE)
                          IPAR= IPAR+ 1
                          QM0(I,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF
c=========================================
c** Now ... update Bv expansion parameters
c=========================================
          IF(IFXBv(ISTATE).LE.0) THEN
              IF(NDEBv(ISTATE).EQ.-1) THEN
c** If use band constants for Rotational (including CDC) constants, but 
c  NOT for Gv ... first do outer loop over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... and then inner loop over vibrational levels.
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
                          IF(NRC(IV,ISTATE,ISOT).GT.0) THEN
                              DO  M= 1,NRC(IV,ISTATE,ISOT)
                                  IPAR= IPAR+ 1
                                  ZK(M,IV,ISTATE,ISOT)= PV(IPAR)
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).GE.2)) THEN
c ... First Dunham Bv parameters for MXS or pure Dunham cases ...
                  IF(LMAX(1,ISTATE).GE.0) THEN
                      DO  L= 0, LMAX(1,ISTATE)
                          IPAR= IPAR+ 1
                          YLM(L,1,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  MMAX= 1
                  MMIN= 0
                  IF((NDEGv(ISTATE).EQ.1).OR.(NDEGv(ISTATE).LT.0).OR.
     1                                    (IFXGv(ISTATE).GT.0)) MMIN=1
                  ENDIF
              IF(NDEBv(ISTATE).GT.0) THEN
                  NEWBv= 1
                  IF(NP1(ISTATE).GT.0) THEN
c ... if appropriate, update the Bv NDE numerator polynomial coeffts.
                      DO  I= 1, NP1(ISTATE)
                          IPAR= IPAR+ 1
                          PM1(I,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  IF(NQ1(ISTATE).GT.0) THEN
c... if appropriate, update the Bv NDE denominator polynomial coeffts.
                      DO  I= 1, NQ1(ISTATE)
                          IPAR= IPAR+ 1
                          QM1(I,ISTATE)= PV(IPAR)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF
c** Call subroutine to generate updated NDE-type Gv & Bv values AND to
c   generate core of partial derivatives for next cycle.
          IF(((NDEGv(ISTATE).GT.0).AND.(IFXGv(ISTATE).LE.0)).OR.
     1           ((NDEBv(ISTATE).GT.0).AND.(IFXBv(ISTATE).LE.0))) THEN
              CALL NDEDGB(ISTATE,NISTP,NEWGv,NEWBv,RSQMU,VMAX(ISTATE))
              IF(NDEGv(ISTATE).EQ.1) Te(ISTATE)= ZK(0,-1,ISTATE,1)
              ENDIF
c=============================================
c** If fitting to CDC's, update parameters ...
c=============================================
          IF(IFXCDC(ISTATE).LE.0) THEN
              IF((NDECDC(ISTATE).EQ.-1).AND.(NDEBv(ISTATE).GE.0)) THEN
c** If use band constants for CDCs, but NOT for Bv (or Gv)
c ... first do outer loop over isotopomers ...
                  DO  ISOT= 1,NISTP
c ... and then inner loop over vibrational levels.
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
                          IF(NRC(IV,ISTATE,ISOT).GT.1) THEN
                              DO  M= 2,NRC(IV,ISTATE,ISOT)
                                  IPAR= IPAR+ 1
                                  ZK(M,IV,ISTATE,ISOT)= PV(IPAR)
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDDO
c=====end of section for fit to band constants for CDCs & doubling======
                  ENDIF
              IF(NDECDC(ISTATE).EQ.0) THEN
c** If fitting to Dunham expansions for the CDCs
                  IF(MMIN.LT.0) MMIN= 2
                  MMAX= NCDC(ISTATE)+ 1
                  DO  M= 2,MMAX
                      IF(LMAX(M,ISTATE).GE.0) THEN
                          DO  L= 0, LMAX(M,ISTATE)
                              IPAR= IPAR+ 1
                              YLM(L,M,ISTATE)= PV(IPAR)
                              ENDDO
                          ENDIF
                      ENDDO
                  ENDIF
              ENDIF
c
          IF(MMIN.GE.0) THEN
c=======================================================================
c** Now generate all relevant Dunham band constants  ZK(M,v,ISTATE,ISOT)
c=======================================================================
              DO  IV=VMIN(ISTATE), VMAX(ISTATE)
                  DO  ISOT= 1,NISTP
                      XX= (IV+ 0.5d0)*RSQMU(ISOT)
                      IF(NDEGv(ISTATE).GE.2) THEN
                          SwLR= dexp((IV- VSISO(ISTATE,ISOT))/
     1                                            DVSISO(ISTATE,ISOT))
                          Sw= 1.d0/(1.d0+ SwLR)
                          SwLR= SwLR*Sw
                          ENDIF
                      DO  M= MMIN,MMAX
                          IF(LMAX(M,ISTATE).GE.0) THEN
                              YY= 0.d0
                              DO  L= LMAX(M,ISTATE),0,-1
                                  YY= YY*XX + YLM(L,M,ISTATE)
                                  ENDDO
                              YY= YY*RMUP(M,ISOT)
                              IF(((M.EQ.0).AND.(NDEGv(ISTATE).GE.2))
     1                   .OR.((M.EQ.1).AND.(NDEBv(ISTATE).GE.2))) THEN
                                  ZK(M,IV,ISTATE,ISOT)= YY*Sw +
     1                                       SwLR*ZK(M,IV,ISTATE,ISOT)
                                ELSE
                                  ZK(M,IV,ISTATE,ISOT)= YY
                                ENDIF
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDDO
              ENDIF
c=====================================================================
c*** If appropriate, update Lambda/Gamma doubling expansion parameters
c=====================================================================
   30     IF((IOMEG(ISTATE).NE.0).AND.(IFXLD(ISTATE).LE.0)
     1                                 .AND.(NLDMX(ISTATE).GE.1)) THEN
              MQ0= MAX0(0,IOMEG(ISTATE)-1)
              IF(NDELD(ISTATE).EQ.-1) THEN
c ... if using Band Constant form for doubling constants 
                  DO  ISOT= 1,NISTP
c ... first do outer loop over isotopomers ...
                      DO  IV= VMIN(ISTATE), VMAX(ISTATE)
c ... and then inner loop over vibrational levels.
                          IF(NQC(IV,ISTATE,ISOT).GT.0) THEN
                              DO  M= 1,NQC(IV,ISTATE,ISOT)
                                  IPAR= IPAR+ 1
                                  ZQ(M+MQ0,IV,ISTATE,ISOT)= PV(IPAR)
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF(NDELD(ISTATE).GE.0) THEN
c ... if using Dunham-type representation for doubling constants ...
                  DO  M= 1,NLDMX(ISTATE)
                      MQM= MQ0+ M
                      IF(LDMAX(MQM,ISTATE).GE.0) THEN
                          DO  L= 0, LDMAX(MQM,ISTATE)
                              IPAR= IPAR+ 1
                              QLM(L,MQM,ISTATE)= PV(IPAR)
                              ENDDO
c ... then generate values of the resulting Lambda doubling coeffts. for
c   each vibrational level of each isotopomer.
                          DO  IV= VMIN(ISTATE),VMAX(ISTATE)
                              DO   ISOT= 1,NISTP
                                  XX= (IV+ 0.5d0)*RSQMU(ISOT)
                                  YY= QLM(0,MQM,ISTATE)
                                  IF(LDMAX(MQM,ISTATE).GE.1) THEN
                                      XXP= 1.d0
                                      DO   L= 1,LDMAX(MQM,ISTATE)
                                          XXP= XXP*XX
                                          YY= YY+ QLM(L,MQM,ISTATE)*XXP
                                          ENDDO
                                      ENDIF
                                  IF(IOMEG(ISTATE).LT.0)
     1                           ZQ(M,IV,ISTATE,ISOT)= YY*RMUP(M,ISOT)
                                  IF(IOMEG(ISTATE).GT.0)
     1            ZQ(MQM,IV,ISTATE,ISOT)= YY*RMUP(1,ISOT)**(MQ0+MQM+1)
                                  ENDDO
                              ENDDO
                          ENDIF
                      ENDDO
                  ENDIF
              ENDIF
c
c** If appropriate, update B-O-B  delta  expansion parameters
c=======================================================================
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).LE.0)) THEN
c ... first ... for atom-A
              ATOM2= 2
              IF(AN(1).EQ.AN(2)) ATOM2= 1
              DO  ATOM= 1,ATOM2
                  LAMIN= 0
                  IF((ISTATE.EQ.1).AND.(BOB00.LE.0)) LAMIN= 1
                  DO  M= 0,BOBORD(ISTATE)
                      IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                          DO  L= LAMIN, LAMAX(ATOM,M,ISTATE)
                              IPAR= IPAR+ 1
                              DELTA(ATOM,L,M,ISTATE)= PV(IPAR)
                              ENDDO
                          ENDIF
                      LAMIN= 0
                      ENDDO
                  ENDDO
c ... & then use these  delta's  to update the isotopomeric band constants
              CONTINUE
              DO  ISOT= 1,NISTP
                  DO  IV= 0, VMAX(ISTATE)
                      LAMIN= 0
                      IF((ISTATE.EQ.1).AND.(BOB00.LE.0)) LAMIN= 1
                      DO  M= 0,BOBORD(ISTATE)
                          YY= 0.d0
                          ZATOM= 1.d0 - ZMASS(1,1)/ZMASS(1,ISOT)
                          IF(AN(1).EQ.AN(2)) ZATOM= ZATOM+ 1.d0 -
     1                                        ZMASS(2,1)/ZMASS(2,ISOT)
                          DO  ATOM= 1,ATOM2
                              IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                                  DO  L= LAMIN,LAMAX(ATOM,M,ISTATE)
                                      YY= YY+VPHPW(IV,L)*RSQMUP(L,ISOT)*
     1                                    DELTA(ATOM,L,M,ISTATE)*ZATOM
                                      ENDDO
                                  ENDIF
                              ZATOM= 1.d0 - ZMASS(2,1)/ZMASS(2,ISOT)
                              ENDDO
                          ZK(M,IV,ISTATE,ISOT)= ZK(M,IV,ISTATE,ISOT) +
     1                                                 YY*RMUP(M,ISOT)
                          LAMIN= 0
                          ENDDO
                      ENDDO
                  ENDDO
              ENDIF
   50     CONTINUE
      IF(NFSTOT.GT.0) THEN
c** If appropriate, map onto values of fluorescence series band origins.
          DO  IFS= 1,NFSTOT
              IPAR= IPAR+ 1
              ORIGIN(IFS)= PV(IPAR)
              ENDDO
          ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PREDICT(IDAT,YC,NCDC,NLDMX,IOMEG,efREF,DLIMIT,ORIGIN,
     1                                                       PV,ZK,ZQ)
c** Subroutine using existing band constants ZK, lambda doubling
c  constants ZQ to calculate the value of datum-IDAT.           19/04/05
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      REAL*8 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     1  ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX),DLIMIT(NSTATEMX),PV(NPARMX),
     2  YC,JJP,JJPQ,JJPP,JJPPQ,JJPW,JJPPW
      INTEGER  IDAT,IBB,ESP,ESPP,IVP,IVPP,ISOT,M,MMAX,MQ0,
     1  NCDC(NSTATEMX),NLDMX(NSTATEMX),IOMEG(NSTATEMX),efREF(NSTATEMX) 
c
c** Type statements & common block for data
cc
cc    REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
cc   1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
cc   2  RMUP(0:9,NISTPMX)
cc    INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
cc   1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
cc   2  EFP(NDATAMX),EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),
cc   3  FSBAND(NBANDMX),NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),
cc   4  ISTP(NBANDMX),IFIRST(NBANDMX),ILAST(NBANDMX),
cc   5  NTV(NSTATEMX,NISTPMX)
cc    CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX)
cc    COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,COUNTOT,
cc   1 NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,TVUP,TVLW,VP,VPP,
cc   2 FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
      IBB= IB(IDAT)
      ESP= IEP(IBB)
      ESPP= IEPP(IBB)
      IVP= VP(IBB)
      IVPP= VPP(IBB)
      IF(ESP.LE.0) IVPP= JP(IDAT)
      ISOT= ISTP(IBB)
      IF(ESP.EQ.-3) THEN
c** For input state s=ESPP Bv value, return current predicted value ...
          YC = ZK(1,IVPP,ESPP,ISOT)
          RETURN
          ENDIF 
c
c** First treat lower level of the transition .........
c=======================================================================
      IF(TVLW(IDAT).GT.0) THEN
c ... if it is represented by an individual fitted term value
c [TVLW is parameter counter: Lower level term value for transition IDAT]
c------------------------------------------------------------
          YC= - PV(TVLW(IDAT))
          IF(ESP.EQ.0)  YC= YC+ ORIGIN(NFS(IBB))
          IF(ESP.EQ.-1) YC= YC+ DLIMIT(ESPP)
          IF(ESP.LE.0) RETURN
        ELSE
c Otherwise - first define the lower state centrifugal factor 
          JJPPQ= JPP(IDAT)*(JPP(IDAT)+ 1)
          IF(IOMEG(ESPP).LE.0) THEN
              JJPP= JJPPQ
            ELSE
              JJPP= JJPPQ - IOMEG(ESPP)**2
            ENDIF
          YC= 0.d0
c ... then, for FS or PAS datum or term-value upper level, sum
c------------------------------- lower-level mechanical rotation terms
          IF((ESP.LE.0).OR.(TVUP(IDAT).GT.0)) THEN
              MMAX= NCDC(ESPP)+ 1
              JJPPW= 1.d0   
              YC= - ZK(0,IVPP,ESPP,ISOT)
              DO  M= 1,MMAX
                  JJPPW= JJPPW*JJPP
                  YC= YC - ZK(M,IVPP,ESPP,ISOT)*JJPPW     
                  ENDDO
              ENDIF
          ENDIF
c
c** Now ... treat upper level or origin of transition ...
c=======================================================================
      IF(TVUP(IDAT).GT.0) THEN
c ... if it is represented by an individual fitted term value
          YC= YC+ PV(TVUP(IDAT))
          IF(IOMEG(ESPP).NE.0) GOTO 20
          RETURN
          ENDIF
      IF((ESP.EQ.0).OR.(ESP.EQ.-1)) THEN
c ... if it is a fluorescence series or PAS datum ...
          IF(ESP.EQ.0)  YC= YC + ORIGIN(NFS(IBB))
          IF(ESP.EQ.-1) YC= YC + DLIMIT(ESPP)
          ENDIF
      IF(ESP.GT.0) THEN
c** For mechanical upper-state level, first define rotational factor ...
          JJPQ= JP(IDAT)*(JP(IDAT)+ 1)
          IF(IOMEG(ESP).LE.0) THEN
              JJP=JJPQ
            ELSE
              JJP= JJPQ - IOMEG(ESP)**2
            ENDIF
          IF(TVLW(IDAT).GT.0) THEN
c ... sum mechanical rotation terms if lower level given by term value
              MMAX= NCDC(ESP)+ 1
              JJPW= 1.d0
              YC= YC+ ZK(0,IVP,ESP,ISOT)
              DO  M= 1, MMAX
                  JJPW= JJPW* JJP
                  YC= YC+ ZK(M,IVP,ESP,ISOT)*JJPW
                  ENDDO
            ELSE
c ...  else minimize truncation errors using joint upper/lower rot. sums
              MMAX= MAX(NCDC(ESP),NCDC(ESPP))+ 1
              JJPW= 1.d0
              JJPPW= 1.d0
              YC= ZK(0,IVP,ESP,ISOT)- ZK(0,IVPP,ESPP,ISOT)
              IF(MMAX.GE.1) THEN
                  DO  M= 1,MMAX
                      JJPW= JJPW*JJP
                      JJPPW= JJPPW*JJPP
                      YC= YC+ ZK(M,IVP,ESP,ISOT)*(JJPW- JJPPW) -
     1                (ZK(M,IVPP,ESPP,ISOT)- ZK(M,IVP,ESP,ISOT))*JJPPW
                      ENDDO
                  ENDIF
            ENDIF
          ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c** If appropriate include lower-level Lambda or ^2\Sigma doubling shift
c-----------------------------------------------------------------------
c If considering doublet Sigma e/f splitting (gamma-doubling)
c   e par is   +1/2  N    { g1 + g2 N(N+1) + g3 [N(N+1)]^2 + ...} 
c   f par is   -1/2 (N+1) { g1 + g2 N(N+1) + g3 [N(N+1)]^2 + ...} 
c while for Lambda doubling, for both parities (EFPP= +/-1 for e/f)
c  Delta(E)= (1/2)(EFPP - efREF) [J(J+1)] { q1 + q2[J(J+1) - OMEGA^2] 
c                                       + q3 [J(J+1) - OMEGA^2]^2 + ...}
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   20 IF((IOMEG(ESPP).NE.0).AND.(EFPP(IDAT).NE.0).AND.
     1                                        (NLDMX(ESPP).GE.1)) THEN
          MQ0= MAX0(0,IOMEG(ESPP)-1)
          IF(IOMEG(ESPP).LT.0) THEN
              JJPPW= JPP(IDAT)- MIN(0,EFPP(IDAT))
            ELSE
              JJPPW= JJPPQ**IOMEG(ESPP) 
            ENDIF
          DO  M= 1,NLDMX(ESPP)
              YC= YC - 0.5d0*(EFPP(IDAT) - efREF(ESPP))*
     1                                  ZQ(M+MQ0,IVPP,ESPP,ISOT)*JJPPW
              JJPPW= JJPPW*JJPP
              ENDDO
          ENDIF
c ... then deal with upper state Lambda or ^2\Sigma doubling, if present
      IF(ESP.GT.0) THEN
          IF((IOMEG(ESP).NE.0).AND.(EFP(IDAT).NE.0)
     1                                    .AND.(NLDMX(ESP).GE.1)) THEN
              MQ0= MAX0(0,IOMEG(ESP)-1)
              IF(IOMEG(ESP).LT.0) THEN
                  JJPW= JP(IDAT)- MIN(0,EFP(IDAT)) 
                ELSE
                  JJPW= JJPQ**IOMEG(ESP)
                ENDIF
              DO  M= 1,NLDMX(ESP)
                  YC= YC + 0.5d0*(EFP(IDAT)- efREF(ESP))*
     1                                     ZQ(M+MQ0,IVP,ESP,ISOT)*JJPW
                  JJPW= JJPW*JJP
                  ENDDO
              ENDIF
          ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PPISOT(NISTP,AN,MN,PV,PU,PS,CM,ZMASS,RSQMUP,RMUP,NAME,
     1                                                 SLABL,NAMEPARM)
c** Subroutine to create and print maximally rounded (to Sensitivity)
c  Dunham Ylm or NDE parameters for minority isotopomers (ISOT>1) .
c** While this subroutine only (currently) works with Dunham-type YLM
c  (or qLM) parameters, the looping must consider ALL parameters, to 
c  get the count/labelling correct.              Version date: 15/05/05
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      REAL*8 ZMASS(3,NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),RMUP(0:9,NISTPMX)
      REAL*8  PYLM(0:NDUNMX,0:9,NISTPMX),UYLM(0:NDUNMX,0:9,NISTPMX),
     1 SYLM(0:NDUNMX,0:9,NISTPMX),PqLM(0:NDUNMX,0:9,NISTPMX),
     2 UqLM(0:NDUNMX,0:9,NISTPMX),SqLM(0:NDUNMX,0:9,NISTPMX),PV(NPARMX),
     3 PU(NPARMX),PS(NPARMX),CM(NPARMX,NPARMX),PFCT,FMU,FMU1,D10,D20,
     4 D01,D11, Y00SC(NISTPMX),U00SC(NISTPMX),UZPE(NISTPMX)
      real*8  tst
      INTEGER  I,ILM,IPAR,ISTATE,ISOT,IVIB,IROT,L,LMIN,LAMIN,M,MQ0,MQM,
     1  MMIN,MMAX,MN(2,NISTPMX),NISTP,ATOM,ATOM2,AN(2)
      INTEGER  IYLM(0:NDUNMX,0:9),IqLM(0:NDUNMX,0:9),IDLM(0:NDUNMX,0:9)
      CHARACTER*20 NAMEPARM(NPARMX),NAMEY00
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
c** Type statements and common block for case (type of representation)
c
      REAL*8  XM(0:9,NSTATEMX,NISTPMX),PNDE(0:9,NSTATEMX)
c
      INTEGER  NSTATES,IBAND,VMIN(NSTATEMX),VMAX(NSTATEMX),
     1 NCDC(NSTATEMX),IOMEG(NSTATEMX),NLDMX(NSTATEMX),efREF(NSTATEMX),
     2 MULTPLT(NSTATEMX),NDEGv(NSTATEMX),NDEBv(NSTATEMX),
     3 NDECDC(NSTATEMX),NDELD(NSTATEMX),IFXGv(NSTATEMX),IFXBv(NSTATEMX),
     4 IFXCDC(NSTATEMX),IFXLD(NSTATEMX),BOBORD(NSTATEMX),
     5 NUMNDE(NSTATEMX),IFXD(NSTATEMX),IFXVD(NSTATEMX),ITYPE(NSTATEMX),
     6 NP0(NSTATEMX),NQ0(NSTATEMX),IP0(NSTATEMX),IQ0(NSTATEMX),
     7 ITYPB(NSTATEMX),NP1(NSTATEMX),NQ1(NSTATEMX),IP1(NSTATEMX),
     8 IQ1(NSTATEMX),LMAX(0:9,NSTATEMX),LDMAX(9,NSTATEMX),
     9 IFXVS(NSTATEMX),IFXDVS(NSTATEMX),BOB00,LAMAX(2,0:9,NSTATEMX),
     a IPSTATE(NSTATEMX),NPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     b NQPAR(0:NVIBMX,NSTATEMX,NISTPMX),
     c FITGV(0:NVIBMX,NSTATEMX,NISTPMX),NRC(0:NVIBMX,NSTATEMX,NISTPMX),
     d NQC(0:NVIBMX,NSTATEMX,NISTPMX),NEBC(NSTATEMX)
c
      COMMON /CASEBLK/XM,PNDE, NSTATES,IBAND,VMIN,VMAX,NCDC,IOMEG,NLDMX,
     1 efREF,MULTPLT,NDEGv,NDEBv,NDECDC,NDELD,IFXGv,IFXBv,IFXCDC,IFXLD,
     2 IFXVS,IFXDVS,BOBORD,NUMNDE,IFXD,IFXVD,ITYPE,NP0,NQ0,IP0,IQ0,
     3 ITYPB,NP1,NQ1,IP1,IQ1,LMAX,LDMAX,BOB00,LAMAX,IPSTATE,NPAR,NQPAR,
     4 FITGV,NRC,NQC,NEBC
c
c** Type statements and common block for actual parameter values
c
      REAL*8  Te(NSTATEMX),VPHPW(0:NVIBMX,0:NDUNMX),
     1 YLM(0:NDUNMX,0:9,NSTATEMX),DELTA(2,0:NDUNMX,0:9,NSTATEMX),
     2 QLM(0:NDUNMX,9,NSTATEMX),DLIMIT(NSTATEMX),VD(NSTATEMX),
     3 PM0(NDUNMX,NSTATEMX),QM0(NDUNMX,NSTATEMX),PM1(NDUNMX,NSTATEMX),
     4 QM1(NDUNMX,NSTATEMX),VS(NSTATEMX),DVS(NSTATEMX),
     5 VSISO(NSTATEMX,NISTPMX),DVSISO(NSTATEMX,NISTPMX),ORIGIN(NBANDMX),
     6 ZK(0:9,-1:NVIBMX,NSTATEMX,NISTPMX),
     6 ZQ(9,-1:NVIBMX,NSTATEMX,NISTPMX)
c
      COMMON /PARMBLK/Te,VPHPW,YLM,DELTA,QLM,DLIMIT,VD,PM0,QM0,PM1,QM1,
     1                                VS,DVS,VSISO,DVSISO,ORIGIN,ZK,ZQ
c
      DATA NAMEY00/'    Delta{T(v=-1/2)}'/
c-----------------------------------------------------------------------
      IPAR= 0
ccc   IPARDLIM= -1
      DO 50 ISTATE= 1,NSTATES
c** Treat one electronic state at a time.  Skip over band constant fit
c  cases or cases when Gv & Bv held fixed. 
          DO  ISOT= 1, NISTP
              PYLM(0,0,ISOT)= 0.d0
              UYLM(0,0,ISOT)= 0.d0
              SYLM(0,0,ISOT)= 0.d0
              ENDDO
          IF(NDEGv(ISTATE).LT.0) GO TO 50
          IPAR= IPSTATE(ISTATE)
          DO  M= 0,9
              DO  L= 0,NDUNMX
                  IYLM(L,M)= 0
                  ENDDO
              ENDDO
c** Alternately, if use Dunham or NDE functions for this state ...
          MMIN= -1
          IF(IFXGv(ISTATE).LE.0) THEN
              IF((NDEGv(ISTATE).EQ.0).OR.(NDEGv(ISTATE).GE.2)) THEN
c*** If Dunham parameters are used in pure Dunham or MXS for Gv ...
c ... first - for upper (ISTATE > 1) electronic states, count Te 
                  IF(ISTATE.GT.1) IPAR= IPAR+ 1
                  IVIB= IPAR+ 1
                  DO  L= 1, LMAX(0,ISTATE)
c ... then generate minority isotopomer YLM's
                      IPAR= IPAR+ 1
                      IYLM(L,0)= IPAR
                      DO  ISOT= 1,NISTP
                          PFCT= RMUP(0,ISOT)*RSQMUP(L,ISOT)
                          PYLM(L,0,ISOT)= PV(IPAR)*PFCT
                          UYLM(L,0,ISOT)= (PU(IPAR)*PFCT)**2
                          SYLM(L,0,ISOT)= (PS(IPAR)*PFCT)**2
                          ENDDO
                      ENDDO
                  MMIN= 0
                  MMAX= 0
                  ENDIF
c
              IF(NDEGv(ISTATE).GT.0) THEN
c** If NDE expressions used (in MXS or pure NDE) functions for Gv ...
                  IF(IFXD(ISTATE).LE.0) IPAR= IPAR+1
c ... first, count free  D  &  vD, if appropriate
                  IF(IFXVD(ISTATE).LE.0) IPAR= IPAR+ 1
c ... then count any vibrational numerator polynomial coefficients
                  IF(NP0(ISTATE).GT.0) THEN
                      DO  I= 1,NP0(ISTATE)
                          IPAR= IPAR+ 1
                          ENDDO
                      ENDIF
                  IF(NQ0(ISTATE).GT.0) THEN   
c ... then count any vibrational denominator polynomial coefficients
                      DO  I= 1,NQ0(ISTATE)
                          IPAR= IPAR+ 1
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF
c
          IF(IFXBv(ISTATE).LE.0) THEN
c** Now ... consider/count Bv parameters ...
              IF(NDEBv(ISTATE).EQ.-1) THEN
c** If using band constants for Bv's but not Gv, count to get BOB 
c   corrections right.
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE), VMAX(ISTATE)
                          IF(NRC(I,ISTATE,ISOT).GE.1) IPAR= IPAR + 1
                          ENDDO
                      ENDDO
                  ENDIF
              IF((NDEBv(ISTATE).EQ.0).OR.(NDEBv(ISTATE).GE.2)) THEN
c*** If Dunham parameters are used in pure Dunham or MXS for Bv ...
                  IROT= IPAR+ 1
                  DO  L= 0, LMAX(1,ISTATE)
                      IPAR= IPAR+ 1
                      IYLM(L,1)= IPAR
                      DO  ISOT= 1,NISTP
                          PFCT= RMUP(1,ISOT)*RSQMUP(L,ISOT)
                          PYLM(L,1,ISOT)= PV(IPAR)*PFCT
                          UYLM(L,1,ISOT)= (PU(IPAR)*PFCT)**2
                          SYLM(L,1,ISOT)= (PS(IPAR)*PFCT)**2
                          ENDDO
                      ENDDO
                  MMIN= 0
                  IF((NDEGv(ISTATE).EQ.1).OR.(IFXGv(ISTATE).GT.0))MMIN=1
                  MMAX= 1
                  ENDIF
              IF(NDEBv(ISTATE).GT.0) THEN
c ... and if NDE function is used for Bv's ...
                  IF(NP1(ISTATE).GT.0) THEN
c ... then count any rotational NDE numerator polynomial coefficients
                      DO  I= 1, NP1(ISTATE)
                          IPAR= IPAR+ 1
                          ENDDO
                      ENDIF
                  IF(NQ1(ISTATE).GT.0) THEN
c ... then count any rotational NDE denominator polynomial coefficients
                      DO  I= 1, NQ1(ISTATE)
                          IPAR= IPAR+ 1
                          ENDDO
                      ENDIF
                  ENDIF
              ENDIF
c
          IF(IFXCDC(ISTATE).LE.0) THEN
              IF(NDECDC(ISTATE).EQ.-1) THEN
c** If using band constants for CDC's but not Bv and/or Gv, need 
c  count to get BOB corrections right. 
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE),VMAX(ISTATE)
                          IF(NRC(I,ISTATE,ISOT).GE.2) 
     1                              IPAR= IPAR+ NRC(I,ISTATE,ISOT) - 1
                          ENDDO
                      ENDDO 
                  ENDIF
              IF(NDECDC(ISTATE).EQ.0) THEN
c** If fitting to CDC's using Dunham form, count parameters & prepare ..
                  IF(MMIN.LT.0) MMIN= 2
                  MMAX= NCDC(ISTATE)+ 1
                  DO  M= 2,MMAX
                      DO  L= 0, LMAX(M,ISTATE)
                          IPAR= IPAR+ 1
                          IYLM(L,M)= IPAR
                          DO  ISOT= 1,NISTP
                              PFCT= RMUP(M,ISOT)*RSQMUP(L,ISOT)
                              PYLM(L,M,ISOT)= PV(IPAR)*PFCT
                              UYLM(L,M,ISOT)= (PU(IPAR)*PFCT)**2
                              SYLM(L,M,ISOT)= (PS(IPAR)*PFCT)**2
                              ENDDO
                          ENDDO
                      ENDDO
                  ENDIF
              ENDIF
c
          IF((IOMEG(ISTATE).NE.0).AND.(NDELD(ISTATE).GE.0).AND.
     1             (NLDMX(ISTATE).GT.0).AND.(IFXLD(ISTATE).LE.0)) THEN
c** If fitting to Dunham-type Lambda/Gamma doubling expansion parameters
              MQ0= MAX0(0,IOMEG(ISTATE)-1)
              DO  M= 1,NLDMX(ISTATE)
                  MQM= MQ0+ M
                  IF(LDMAX(MQM,ISTATE).GE.0) THEN
                      DO  L= 0, LDMAX(MQM,ISTATE)
                          IPAR= IPAR+ 1
                          IqLM(L,MQM)= IPAR
                          DO  ISOT= 1, NISTP
                              IF(IOMEG(ISTATE).LT.0) 
     1                              PFCT= RSQMUP(L,ISOT)* RMUP(M,ISOT)
                              IF(IOMEG(ISTATE).GT.0)
     1                     PFCT= PFCT* RMUP(1,ISOT)* RMUP(MQ0,ISOT)**2
                              PqLM(L,MQM,ISOT)= PV(IPAR)*PFCT
                              UqLM(L,MQM,ISOT)= PU(IPAR)*PFCT
                              SqLM(L,MQM,ISOT)= PS(IPAR)*PFCT
                              ENDDO
                          ENDDO
                      ENDIF
                  ENDDO
              ENDIF
c
          IF((BOBORD(ISTATE).GE.0).AND.(IFXGv(ISTATE).LE.0)) THEN
c** If appropriate, count B-O-B  delta  expansion parameters, correct 
c  mass-scaled YLM's and generate appropriate fully correlated uncert.
              DO  M= 0,9
                  DO  L= 0,NDUNMX
                      IDLM(L,M)= 0
                      ENDDO
                  ENDDO
              ATOM2= 2
c** If both atoms are the same chemical species ...
              IF(AN(1).EQ.AN(2)) ATOM2= 1
              DO  ATOM= 1,ATOM2
                  LAMIN= 0
                  IF((ISTATE.EQ.1).AND.(BOB00.LE.0)) LAMIN= 1
                  DO  M= 0,BOBORD(ISTATE)
                      IF(LAMAX(ATOM,M,ISTATE).GE.LAMIN) THEN
                          DO  L= LAMIN, LAMAX(ATOM,M,ISTATE)
                              IPAR= IPAR+ 1
                              ILM= IYLM(L,M)
                              IF(ATOM.EQ.1) IDLM(L,M)= IPAR
c ... correct isotopomeric YLM's for delta contributions
                              DO  ISOT= 1, NISTP
                                  PFCT= RMUP(M,ISOT)*RSQMUP(L,ISOT)
                                  FMU= PFCT*(1.d0- ZMASS(ATOM,1)/
     1                                               ZMASS(ATOM,ISOT))
c ... using a combined mass scaling factor for the homonuclear case
                                  IF(ATOM2.EQ.1) FMU= FMU+ PFCT*
     1                                (1.d0- ZMASS(2,1)/ZMASS(2,ISOT))
                                  PYLM(L,M,ISOT)= PYLM(L,M,ISOT)+ 
     1                                                    FMU*PV(IPAR)
c ... and generate corrected correlated uncertainties & sensitivities
                                  UYLM(L,M,ISOT)= UYLM(L,M,ISOT)+ 
     1                                               (FMU*PU(IPAR))**2
                                  SYLM(L,M,ISOT)= SYLM(L,M,ISOT)+ 
     1                                               (FMU*PS(IPAR))**2
                                  IF(ILM.GT.0) THEN
                                      UYLM(L,M,ISOT)= UYLM(L,M,ISOT)+ 
     1                    2.d0*FMU*PU(IPAR)*PFCT*PU(ILM)*CM(IPAR,ILM)
                                      ENDIF
                                  IF((ATOM.EQ.2).AND.(IDLM(L,M).GT.0))
     1                                                            THEN
c ... including, if appropriate, the Atom-1/Atom-2 cross term
                                      FMU1= (1.d0- ZMASS(1,1)/
     1                                                  ZMASS(1,ISOT))
                                      UYLM(L,M,ISOT)= UYLM(L,M,ISOT)+ 
     1         2.d0*FMU*PU(IPAR)*FMU1*PU(IDLM(L,M))*CM(IPAR,IDLM(L,M))
                                      ENDIF                              
                                  ENDDO
                              ENDDO
                          ENDIF
                      LAMIN= 0
                      ENDDO
                  ENDDO
              ENDIF
c** Now ... round and print the resulting isotopomeric YLM's
          IF((IFXGv(ISTATE).GT.0).OR.(IFXBv(ISTATE).GT.0)) GO TO 50
          IF(NISTP.GT.1) WRITE(6,600) SLABL(ISTATE),((NAME(ATOM),
     1                          MN(ATOM,ISOT),ATOM=1,2),ISOT= 2,NISTP)
          IF(LMAX(0,ISTATE).GT.0) WRITE(6,604) SLABL(1),
     1                             (ZK(0,0,ISTATE,ISOT),ISOT= 1,NISTP)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c** Calculate & print uncertainties in [G(v=0)-G(v=-1/2)]
ccc  RJL should re-think what this is supposed to be doing, and why!
ccc
ccc       IF(LMAX(0,ISTATE).GT.0) THEN
ccc           DO  ISOT= 1,NISTP
ccc               UZPE(ISOT)= 0.d0
ccc               DO  L= 1, LMAX(0,ISTATE)
ccc                   PFCT= 0.d0
ccc                   DO  I= 1, LMAX(0,ISTATE)
ccc                       PFCT= PFCT + VPHPW(0,I)*RSQMUP(I,ISOT)*
ccc  1                            UYLM(I,0,ISOT)*CM(IVIB+L-1,IVIB+I-1)
ccc                       ENDDO
ccc                   UZPE(ISOT)= UZPE(ISOT)
ccc  1                  + PFCT*VPHPW(0,L)*RSQMUP(L,ISOT)*UYLM(I,0,ISOT)
ccc                   ENDDO
ccc               UZPE(ISOT)= DSQRT(UZPE(ISOT))
ccc               ENDDO
ccc           WRITE(6,606) (UZPE(ISOT), ISOT= 1,NISTP)
ccc  606 FORMAT(/' Uncertainty in  [G(v=0)-G(v=-1/2)]  for reference isotop
ccc     1omer:',F10.6:/5x,'and others:',6F10.6:/(16x,6F10.6:))
ccc           ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c** Calculate & print Y00(semiclassical) & its uncertainties
          IF((LMAX(0,ISTATE).GE.2).AND.(LMAX(1,ISTATE).GE.1)) THEN
              DO  ISOT= 1, NISTP
                  PFCT= PYLM(1,0,ISOT)*PYLM(1,1,ISOT)/
     1                                          (12.d0*PYLM(0,1,ISOT))
                  Y00SC(ISOT)= 0.25D0*(PYLM(0,1,ISOT)+ PYLM(2,0,ISOT))
     1                                 - PFCT + PFCT**2/PYLM(0,1,ISOT)
                  D10= (-PFCT + 2.d0*PFCT**2/PYLM(0,1,ISOT))
                  D11= DSQRT(UYLM(1,1,ISOT))*(D10/PYLM(1,1,ISOT))
     1                                                 *RSQMUP(3,ISOT)
                  D10= DSQRT(UYLM(1,0,ISOT))*D10/PYLM(1,0,ISOT)
     1                                                 *RSQMUP(1,ISOT)
                  D20= 0.25d0*RMUP(1,ISOT)*DSQRT(UYLM(2,0,ISOT)) 
                  D01= (0.25d0 + PFCT/PYLM(0,1,ISOT) - 3.d0*(PFCT/
     1          PYLM(0,1,ISOT))**2)*RMUP(1,ISOT)*DSQRT(UYLM(0,1,ISOT))
                  PFCT= D10*(D10+ D20*CM(IVIB,IVIB+1)+ D01*CM(IVIB,IROT)
     1         + D11*CM(IVIB,IROT+1)) + D20*(D10*CM(IVIB,IVIB+1) + D20
     2         + D01*CM(IVIB+1,IROT) + D11*CM(IVIB+1,IROT+1)) 
     3         + D01*(D10*CM(IVIB,IROT) + D20*CM(IVIB+1,IROT) + D01 
     4         + D11*CM(IROT,IROT+1)) + D11*(D10*CM(IVIB,IROT+1) 
     5         + D20*CM(IVIB+1,IROT+1) + D01*CM(IROT,IROT+1) + D11)
                  U00SC(ISOT)= DSQRT(PFCT)
                  ENDDO
              WRITE(6,608) (Y00SC(ISOT),U00SC(ISOT),ISOT= 1,NISTP)
              ENDIF
          WRITE(6,603) 
c
c** Now print overall isotopomer Ylm's and their uncertainties!
          IF((MMIN.GE.0).AND.(NISTP.GE.2)) THEN
              DO  M= MMIN,MMAX
                  LMIN= 0
                  IF((M.EQ.0).AND.(ISTATE.EQ.1)) LMIN= 1  
                  IF((M.EQ.0).AND.(ISTATE.GT.1).AND.
     1                             (BOBORD(ISTATE).GE.0)) WRITE(6,612)
                  DO  L= LMIN,LMAX(M,ISTATE)
                      DO  ISOT= 1,NISTP
                          UYLM(L,M,ISOT)= SQRT(UYLM(L,M,ISOT))
                          SYLM(L,M,ISOT)= SQRT(SYLM(L,M,ISOT))
                          IF(ISOT.GT.1) THEN
                             IF(DABS(PYLM(L,M,ISOT)).GT.0.d0) 
     1                      CALL ROUNDSEN(PYLM(L,M,ISOT),SYLM(L,M,ISOT))
                             ENDIF
                          ENDDO
                      IF((M.EQ.0).AND.(L.EQ.0)) THEN
                          WRITE(6,602) NAMEY00,(PYLM(L,M,ISOT),
     1                                   UYLM(L,M,ISOT),ISOT= 2,NISTP)
                          WRITE(9,902) NAMEY00,PYLM(L,M,1),
     1                      UYLM(L,M,1),(PYLM(L,M,ISOT),ISOT= 2,NISTP)
                        ELSE 
                          WRITE(6,602) NAMEPARM(IYLM(L,M)),
     1                   (PYLM(L,M,ISOT),UYLM(L,M,ISOT),ISOT= 2,NISTP)
                          WRITE(9,902) NAMEPARM(IYLM(L,M)),PYLM(L,M,1),
     1                      UYLM(L,M,1),(PYLM(L,M,ISOT),ISOT= 2,NISTP)
                        ENDIF
                      ENDDO
                  WRITE(6,603)
                  ENDDO
              ENDIF
c
c** If appropriate ... round and print the resulting isotopomeric qLM's
          IF((IOMEG(ISTATE).NE.0).AND.(NDELD(ISTATE).GE.0).AND.
     1              (NLDMX(ISTATE).GT.0).AND.(IFXLD(ISTATE).LE.0).AND.
     2                                              (NISTP.GE.2)) THEN
              WRITE(6,610) SLABL(ISTATE),((NAME(ATOM),MN(ATOM,ISOT),
     1                                        ATOM=1,2),ISOT= 2,NISTP)
              WRITE(6,603)
              DO  M= 1,NLDMX(ISTATE)
                  MQM= MQ0+ M
                  IF(LDMAX(MQM,ISTATE).GE.0) THEN
                      DO  L= 0,LDMAX(MQM,ISTATE)
                          DO  ISOT= 1,NISTP
                              IF(ISOT.GT.1) THEN
                                 IF(DABS(PqLM(L,MQM,ISOT)).GT.0.d0) 
     1                CALL ROUNDSEN(PqLM(L,MQM,ISOT),SqLM(L,MQM,ISOT))
                                 ENDIF
                              ENDDO
                          WRITE(6,602) NAMEPARM(IqLM(L,MQM)),
     1               (PqLM(L,MQM,ISOT),UqLM(L,MQM,ISOT),ISOT= 2,NISTP)
                          WRITE(9,902) NAMEPARM(IqLM(L,MQM)),
     1    PqLM(L,MQM,1),UqLM(L,MQM,1),(PqLM(L,MQM,ISOT),ISOT= 2,NISTP)
                          ENDDO
                      ENDIF
                  ENDDO
              WRITE(6,603)
              ENDIF
   50     CONTINUE
      RETURN
  600 FORMAT(/" State-",A3," Sensitivity-Rounded parameters's for Minori
     1ty Isotopomers:"/1x,32('==')/(4(4x,a2,"(",i3,")-",a2,"(",i3,")":)
     2 ))
  602 FORMAT(a20,2(1PD19.11,' (',D7.1,')':)/
     1                               (20x,2(1PD19.11,' (',D7.1,')':)))
  603 FORMAT(' ')
  604 FORMAT(/' Zero point level T(v=0) relative to v= -1/2 of the first
     1 state considered (',A2,')'/4x,'for the reference isotopomer is:',
     2  F15.6:/4x,'and for the others:',4F14.6:/(16x,4F14.6:))
  608 FORMAT(/' Semiclassical  Y00  of the reference isotopomer is:',
     1 F10.6,'(',F9.6,')':/4x,'& of others:',3(F10.6,'(',F9.6,')':)/
     2 (16x,3(F10.6,'(',F9.6,')':)))
  610 FORMAT(/" State-",A3," Sensitivity-Rounded qLM's for Minority Isot
     1opomers:"/(4(4x,a2,"(",i3,")-",a2,"(",i3,")":)))
  612 FORMAT(2x,'Delta{T(v=-1/2)} = [T(v=-1/2;{this isotopomer}) - T(v=-
     11/2;{ref.isotopomer})]')
c 901 FORMAT(A20,10(' &',1PD19.11,' (',D7.1,')'))
  902 FORMAT('10^{0} &',A20,' &',1PD19.11,' (',D7.1,')',9(' &',
     1                                                    1PD19.11:))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE ROUNDSEN(PV,PS)
c** Subroutine to round off parameter with value  PV  at the
c  IROUND'th significant digit of the quantity [its sensitivity] PS . 
c** On return, the rounded value replaced the initial value  PV.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER    IROUND,IRND,KRND 
      REAL*8  PV,PS,CRND,XRND,FCT,UU,CNST
      XRND= 0.d0
      IROUND= 1
      CNST= PV
      UU= CNST
      XRND= DLOG10(PS)
c** First ... fiddle with log's to perform the rounding
      IRND= INT(XRND)            
      IF(XRND.GT.0) IRND=IRND+1
      IRND= IRND- IROUND
      FCT= 10.D0**IRND
      CRND= PV/FCT
      XRND= 0.d0                                                         
      IF(DABS(CRND).GE.1.D+8) THEN
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.2.D+17) THEN
              WRITE(6,601) IROUND
  601 FORMAT(1x,39('==')/' Caution:',i3,'-digit rounding would exceed (a
     1assumed) REAL*8'/' ********   precision overflow at 1.D+16, so kee
     2p unrounded constant')
               RETURN
               ENDIF
c ... to avoid problems from overflow of I*4 integers ...
          KRND= NINT(CRND/1.D+8)
          XRND= KRND*1.D+8
          CRND= CRND-XRND
          XRND= XRND*FCT
          END IF
      IRND= NINT(CRND)
      CNST= IRND*FCT+ XRND
      PV= CNST
c     WRITE(6,614) UU,PS,PV
c 614 FORMAT(1x,30('==')/' Round Off   PV=',1PD22.14,'  with   PS=',
c    1    d9.2/3x,'fixing it as ',D22.14)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

***********************************************************************
      SUBROUTINE DIFFSTATS(NSTATES,ROBUST,MKPRED)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of 15 May 2005
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
c
      INTEGER I,IBB,ISOT,ISTATE,ISTATEE,J,NSTATES,MKPRED,ROBUST
      REAL*8 AVE,AVETOT,DIV,RMSR,RMSTOT,SSQTOT
      CHARACTER*3 MARKER,NEF(-1:1)
c
c** Type statements & common blocks for characterizing transitions
c
      REAL*8  AVEUFREQ(NBANDMX),MAXUFREQ(NBANDMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NBVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NBANDS(NISTPMX),
     7  YPR(NISTPMX,NSTATEMX,7,6,NBANDMX)
c
      COMMON /TYPEBLK/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NBVPP,NWIDTH,
     2  NEBPAS,NBANDS,YPR
c
      DATA  NEF/'  f',' ef','  e'/
c========================================================================
      ISOT= 1
      SSQTOT= 0.d0
      IF(MKPRED.GT.0) WRITE(6,600)
c** Summarize data discrepancies for one isotopomer at a time.
   10 WRITE(6,602) NBANDS(ISOT),(NAME(I),MN(I,ISOT),I= 1,2)
c
c** Loop over bands for each (lower) electronic state, in turm
      DO 90 ISTATE= 1,NSTATES
      IF(NTRANSMW(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          WRITE(6,605)
          WRITE(8,604) NTRANSMW(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,4,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,4,3,I)*AVE
                  WRITE(6,606)YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I)
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSMW(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSMW(ISOT,ISTATE)
          IF(MKPRED.LE.0) WRITE(6,630) NTRANSMW(ISOT,ISTATE),AVETOT,
     1                                                          RMSTOT
	    ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          WRITE(6,605)
          WRITE(8,608) NTRANSIR(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,3,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,3,3,I)*AVE
                  WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I)
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSIR(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSIR(ISOT,ISTATE)
          IF(MKPRED.LE.0) WRITE(6,630) NTRANSIR(ISOT,ISTATE),AVETOT,
     1                                                          RMSTOT
          ENDIF
c
c** Book-keeping for Electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              WRITE(6,605)
              WRITE(8,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              RMSTOT= 0.d0
              AVETOT= 0.d0
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      IF(MKPRED.LE.0) THEN
                          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,
     1                                        RMSR,SSQTOT,DFREQ,UFREQ)
                          RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,2,3,I)*RMSR**2
                          AVETOT= AVETOT+ YPR(ISOT,ISTATE,2,3,I)*AVE
                          WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                  YPR(ISOT,ISTATE,2,1,I),YPR(ISOT,ISTATE,2,3,I),
     2                  YPR(ISOT,ISTATE,2,5,I),YPR(ISOT,ISTATE,2,6,I),
     3                            AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                          ENDIF
                      WRITE(8,605)
                      IF(MKPRED.LE.0) WRITE(8,606) 
     1                  YPR(ISOT,ISTATE,2,2,I),YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3             YPR(ISOT,ISTATE,2,6,I),AVEUFREQ(IBB),MAXUFREQ(IBB),
     4                                                        AVE,RMSR
                      IF(MKPRED.GT.0) WRITE(8,606) 
     1                  YPR(ISOT,ISTATE,2,2,I),YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3                                          YPR(ISOT,ISTATE,2,6,I)
                      CALL PBNDERR(IBB,MKPRED,NEF)
                      ENDIF
                  ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSVIS(ISOT,ISTATEE,ISTATE))
              AVETOT= AVETOT/NTRANSVIS(ISOT,ISTATEE,ISTATE)
              IF(MKPRED.LE.0) WRITE(6,630) 
     1                    NTRANSVIS(ISOT,ISTATEE,ISTATE),AVETOT,RMSTOT
              ENDIF
          ENDDO
c
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,612) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          WRITE(6,617)
          WRITE(8,612) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),
     1                 (NAME(I),MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,1,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,1,3,I)*AVE
              WRITE(6,614) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              WRITE(8,617)
              WRITE(8,614) YPR(ISOT,ISTATE,1,1,I),
     1                  YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSFS(ISOT,ISTATE))
              AVETOT= AVETOT/NTRANSFS(ISOT,ISTATE)
              WRITE(6,632) NTRANSFS(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
c
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS  data
          IBB= YPR(ISOT,ISTATE,7,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,627)
          DO  I= IFIRST(IBB),ILAST(IBB)
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.5.d0) ) marker='*  '
              IF( (DIV.GE.5.d0).AND.(DIV.LT.10.d0) ) marker='** '
              IF( (DIV.GE.10.d0) ) marker='***'
              WRITE(8,628) JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     1                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),MARKER
              ENDDO
          WRITE(6,629)
          WRITE(8,629)
          ENDIF
c
      IF(NBVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Bv  data
          IBB= YPR(ISOT,ISTATE,5,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,618) JP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,618) JP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          IBB= YPR(ISOT,ISTATE,6,4,1)    
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),
     1                             UFREQ(J),DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c** End of loop over the various (lower) electronic states
   90 CONTINUE
c=======================================================================
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 10
          ENDIF 
      RMSR= DSQRT(SSQTOT/COUNTOT)
      WRITE(6,624) COUNTOT,RMSR
      RETURN
  600 FORMAT(/1x,36('**')/'  Write to Channel-8 Predictions From Complet
     1e Set of Input Parameters!'/1x,36('**'))
  602 FORMAT(/1x,21('===')/'  *** Discrepancies for',I5,' bands/series o
     1f ',A2,'(',I3,')-',A2,'(',I3,') ***'/1x,21('==='))
  604 FORMAT(/1x,21('===')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' vib. levels')
  605 FORMAT(1x,16('==='),'== Avge. ========'/"   v' ",
     2  ' v" #data  J"min  J"max  Av.Unc.  Max.Unc.   Err/Unc   DRMSD'/
     1  1x,13('-----'))
  606 FORMAT(2I4,I6,3x,I4,3x,I4,1x,1P2D9.1,0PF11.5,F8.3)
  608 FORMAT(/1x,63('=')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A3,'}--{State ',A3,'} Transitions in',i4,' bands')
  612 FORMAT(/1x,75('=')/I5,' Fluorescence transitions into State ',A3,
     1 2x,A2,'(',I3,')-',A2,'(',I3,') in',i5,' series')
  617 FORMAT(1x,52('='),'= Avge. ',15('=')/"   v'  j' p' ",
     2  '#data  v"min  v"max','  AvgeUnc  Max.Unc.   Err/Unc   DRMSD'/
     3  1x,25('---'))
  614 FORMAT(2I4,A3,I6,2I7,1x,1P2D9.1,0PF11.5,F8.3)
  616 FORMAT(/1x,66('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Bv values treated as independent data'/1x,20('=='),
     2 '  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  618 FORMAT(I5,A3,2x,F12.8,1PD9.1,0PF13.8,F12.4)
  620 FORMAT(/1x,73('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,59('=')/'   v   J  p     Width',7x,'u(Width)  [calc-obs]  [cal
     5c-obs]/unc'/1x,59('-'))
  622 FORMAT(2I4,A3,1PD14.6,D10.1,2D13.2)
  624 FORMAT(/1x,29('==')/' For overall fit to',i6,' data,  DRMS(deviati
     1ons)=',G11.4/1x,30('==')) 
  626 FORMAT(/1x,29('==')/I5,' PAS Binding Energies for State ',A3,2x,
     1 A2,'(',I3,')-',A2,'(',I3,')'/1x,50('='),' Avge. ',('=')/
     2 ' #data  v_min  v_max   AvgeUnc  Max.Unc.  Err/Unc  DRMSD'/
     3  1x,29('--')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3)
  627 FORMAT(1x,48('='),'  calc-obs'/'   v   j  p      PAS(Eb)     u(Eb)
     1      calc-obs   /u(FREQ)'/1x,29('--'))
  628 FORMAT(2I4,A3,F14.6,1PD10.1,D13.2,0PF11.4,1X,A3)
  629 FORMAT(1x,29('=='))
  630 FORMAT(1x,7('--'),' For these',i6,' lines, overall:',F11.5,F8.3)
  632 FORMAT(1x,17('-'),' For these',i6,' lines, overall:',F11.5,F8.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE BNDERR(FIRST,LAST,ROBUST,AVEDD,RMSDD,SSQTOT,DFREQ,
     1                                                          UFREQ)
c** Calculate the average (AVEDD) & the root mean square dimensionless 
c  deviation (RSMDD) for the band running from datum # FIRST to LAST.
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      REAL*8 DFREQ(NDATAMX),UFREQ(NDATAMX),AVEDD,RMSDD,SSQTOT
      INTEGER FIRST,LAST,NDAT,I,ROBUST
c
      AVEDD= 0.d0
      RMSDD= 0.d0
      DO  I= FIRST,LAST
          AVEDD= AVEDD+ DFREQ(I)/UFREQ(I)
          IF(ROBUST.LE.0) RMSDD= RMSDD+ (DFREQ(I)/UFREQ(I))**2
          IF(ROBUST.GT.0) RMSDD= RMSDD+ DFREQ(I)**2/
     1                                (UFREQ(I)**2 + DFREQ(I)**2/3.d0)
          ENDDO
      SSQTOT= SSQTOT+ RMSDD
      NDAT= LAST-FIRST+1
      AVEDD= AVEDD/NDAT
      RMSDD= DSQRT(RMSDD/NDAT)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PBNDERR(IBB,MKPRED,NEF)
c** Print to channel-8 a listing of the [obs.-calc.] values for the band
c  running from datum # FIRST to LAST.                           
cc    INCLUDE 'arrsizes.h'             
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** BLOCK DATE Utility routine 'arrsizes.h' governing array dimensioning
c   in dParFiT that MUST be installed under this name in the same
c   (sub)directory containing the folowing FORTRAN file for Program
c    dParFit16 when it is being compiled,
c-----------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NDATAMX,NBANDMX,NVIBMX,NSTATEMX,NDUNMX,
     1   NROTMX
c*  NISTPMX  is the maximum number of isotopomers allowed for
      PARAMETER (NISTPMX = 10)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for
      PARAMETER (NPARMX  = 3000)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 22000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2700)
c*  NDUNMX  is the maximum number of Dunham/NDE power series coeffts.
      PARAMETER (NDUNMX   = 20)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX = 155)
c** NROTMX  is the maximum number of rotational (J or N) values for a
c         given vib level.  Required for term-value fit data counting
      PARAMETER (NROTMX = 200)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
cc    PARAMETER (NDUNMX=0)    % when used wity DPotFit
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX)
c** Differs from PotFit version because these factors not needed.
cc   2 ,ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
cc   3 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,AN,MN,NISTP
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX)
c
      INTEGER  COUNTOT,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),EFP(NDATAMX),
     2 EFPP(NDATAMX),TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NBANDMX),
     3 NFS(NBANDMX),IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),
     4 IFIRST(NBANDMX),ILAST(NBANDMX),NTV(NSTATEMX,NISTPMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*3 SLABL(-6:NSTATEMX)
c
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,COUNTOT,NFSTOT,NBANDTOT,
     1 IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,NFS,IEP,IEPP,ISTP,
     2 IFIRST,ILAST,NTV, NAME,SLABL
c=======================================================================
      REAL*8 DIV
      INTEGER IBB,I,MKPRED
      CHARACTER*3 marker, NEF(-1:1)
c----------------------------------------------------------------------- 
      IF(MKPRED.LE.0) WRITE(8,600)
      IF(MKPRED.GT.0) WRITE(8,601)
      DO  I= IFIRST(IBB),ILAST(IBB)
          IF(MKPRED.LE.0) THEN
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.4.d0) ) marker='*  '
              IF( (DIV.GE.4.d0).AND.(DIV.LT.8.d0) ) marker='** '
              IF( (DIV.GE.8.d0) ) marker='***'
              IF(IEP(IBB).GT.0) WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),
     1         VPP(IBB),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                      DFREQ(I)/UFREQ(I),marker
              IF(IEP(IBB).EQ.0) WRITE(8,602) VP(IBB),VPP(IBB),
     1                  NEF(EFP(I)),JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     2                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),marker
            ELSE
              WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),VPP(IBB),JPP(I),
     1                                           NEF(EFPP(I)),DFREQ(I)
c* Print predictions in alternate (Lyon) format
c             WRITE(11,606)VP(IBB),VPP(IBB),JPP(I),JP(I)-JPP(I),DFREQ(I)
c 606 FORMAT(2I4,I5,I4,f13.4)
            ENDIF
          ENDDO
      WRITE(8,604)
      RETURN                       
  600 FORMAT(1x,59('='),'  calc-obs'/ "   v'  J' p'",
     1  '  v"  J" p"    FREQ(obs)     u(FREQ)    calc-obs  /u(FREQ)'/
     2  1x,69('-'))
  601 FORMAT(1x,36('=')/ "   v'  J' p'",'  v"  J" p"   FREQ(calc)'/
     1  1x,36('-'))
  602 FORMAT(2(2I4,A3),f14.6,2f12.6,f10.4,1x,A3)
  604 FORMAT(1x,69('-'))
      END   
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

C***********************************************************************
      SUBROUTINE GPROUND(IROUND,NPTOT,NPMAX,NPAR1,NPAR2,LPRINT,IFXP,
     1                                                          PV,PU)
c** Subroutine to round off parameters PV(i), i= NPAR1 to NPAR2, at the
c  |IROUND|'th significant digit of the smallest of their uncertainties
c  min{U(i)}.  This procedure does NOT attempt to correct the remaining
c  parameters to compensate for these changes (as ROUND does), so this
c  procedure is not appropriate for nonlinear parameters.
c** On return, the rounded values replaces the initial values of  PV(i).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000-2004  by  Robert J. Le Roy             +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c                      Version of 27 January 2004                      +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  IROUND,NPMAX,NPTOT,NPAR1,NPAR2,NPARM,IRND,KRND,LPRINT
      INTEGER  IFXP(NPTOT)
      REAL*8  PV(NPMAX),PU(NPMAX),CNST,CRND,XRND,FCT,XX,YY,UNC
c
c** Loop over & round off the parameters # NPAR1 to NPAR2 
      IF(LPRINT.GE.2) WRITE(6,602)  NPAR2-NPAR1+1,NPTOT,NPAR1,NPAR2
      UNC= 99.d99
      DO  NPARM= NPAR1, NPAR2
          IF(PU(NPARM).LT.UNC) UNC= PU(NPARM)
          ENDDO
      DO  NPARM= NPAR1, NPAR2
c** First ... fiddle with log's to perform the rounding
          XRND= DLOG10(UNC)
          IRND= INT(XRND)
          IF(XRND.GT.0) IRND=IRND+1
          IRND= IRND- IROUND
          FCT= 10.D0**IRND
          CNST= PV(NPARM)
          YY= CNST
          CRND= PV(NPARM)/FCT
          XRND= 0.d0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.1.D+16) THEN
              WRITE(6,600) IROUND,NPARM
               RETURN
               ENDIF
          IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
              KRND= NINT(CRND/1.D+8)
              XRND= KRND*1.D+8
              CRND= CRND-XRND
              XRND= XRND*FCT
              END IF
          IRND= NINT(CRND)
          CNST= IRND*FCT+ XRND
          PV(NPARM) = CNST
          IFXP(NPARM)= 1
          IF(LPRINT.GE.2) WRITE(6,604) NPARM,YY,PV(NPARM)
  604 FORMAT(5x,'Round parameter #',i4,' from',G20.12,'  to',G20.12)
          ENDDO
          NPARM= NPARM- 1
      RETURN
  600 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
  602 FORMAT(' Rounding off ',i5,' of the ',i5,' parameters #:',i5,
     1 ' to',i5)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,CYCMAX,IROUND,ROBUST,LPRINT,
     1                      IFXP,YO,YU,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting 
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].         25/03/16
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2016  by  Robert J. Le Roy                +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program uses orthogonal decomposition of the "design" (partial 
c  derivative) matrix for the core locally linear (steepest descent) 
c  step, following a method introduced (to me) by Dr. Michael Dulick. 
c** If no parameters are free (NPTOT=0), simply return RMS(residuals) as
c  calculated from the input parameter values {PV(j)}.
c** A user MUST SUPPLY subroutine  DYIDPJ  to generate the predicted
c  value of each datum and the partial derivatives of each datum w.r.t.
c  each parameter (see below) from the current trial parameters.
c  
c** On entry: 
c    NDATA  is the number of data to be fitted 
c    NPTOT  the total number of parameters in the model (.le.NPMAX).
c           If NPTOT.le.0 , assume  YD(i)=YO(i)  and calculate the (RMS 
c           dimensionless deviation)=DSE  from them & YU(i) 
c    NPMAX is the maximum number of model parameters allowed by current
c          external array sizes.  Should set internal NPINTMX = NPMAX 
c          (may be freely changed by the user).
c    CYCMAX is the upper bound on the allowed number of iterative cycles 
c    IROUND .ne. 0  causes Sequential Rounding & Refitting to be 
c             performed, with each parameter being rounded at the 
c            |IROUND|'th sig. digit of its local incertainty.
c        > 0  rounding selects in turn remaining parameter with largest
c             relative uncertainy
c        < 0  round parameters sequentially from last to first
c        = 0  simply stops after full convergence (without rounding).
c    ROBUST > 0  causes fits to use Watson's ``robust'' weighting  
c        1/[u^2 +{(c-o)^2}/3].  ROBUST > 1 uses normal 1/u^2 on first
c        fit cycle and  'robust' on later cycles.
c    LPRINT  specifies the level of printing inside NLLSSRR
c          if: =  0, no print except for failed convergence.
c               < 0  only converged, unrounded parameters, PU & PS's
c              >= 1  print converged parameters, PU & PS's
c              >= 2  also print parameter change each rounding step
c              >= 3  also indicate nature of convergence
c              >= 4  also print convergence tests on each cycle
c              >= 5  also parameters changes & uncertainties, each cycle
c    IFXP(j)  specifies whether parameter  j  is to be held fixed
c           [IFXP > 0] or to be freely varied in the fit [IFXP= 0]
c    YO(i)  are the NDATA 'observed' data to be fitted  
c    YU(i)  are the uncertainties in these YO(i) values
c    PV(j)  are initial trial parameter values (for non-linear fits);  
c           should be set at zero for initially undefined parameters.
c
c** On Exit:   
c    YD(i)  is the array of differences  [Ycalc(i) - YO(i)]
c    PV(j)  are the final converged parameter values
c    PU(j)  are 95% confidence limit uncertainties in the PV(j)'s
c    PS(j)  are 'parameter sensitivities' for the PV(j)'s, defined such 
c           that the RMS displacement of predicted data  due to rounding
c           off parameter-j by PS(j) is .le. DSE/10*NPTOT
c    CM(j,k)  is the correlation matrix obtained by normalizing variance
c           /covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c    TSTPS = max{|delta[PV(j)]/PS(j)|}  is the parameter sensitivity 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    TSTPU = max{|delta[PV(j)]/PU(j)|}  is the parameter uncertainty 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    DSE    is the predicted (dimensionless) standard error of the fit
c
c  NOTE that the squared 95% confidence limit uncertainty in a property 
c  F({PV(j)}) defined in terms of the fitted parameters {PV(j)} (where
c  the L.H.S. involves  [row]*[matrix]*[column]  multiplication) is:
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c
c** Externally dimension:  YO, YU and YD  .ge. NDATA 
c             PV, PU  and  PS  .ge.  NPTOT (say as NPMAX), 
c             CM   as a square matrix with column & row length  NPMAX
c***********************************************************************
      INTEGER MXPdim             !! internal limit on max # parameters
      PARAMETER (MXPdim=3000)    !! must be .GE. external max # NPMAX
      INTEGER I,J,K,L,IDF,ITER,NITER,CYCMAX,IROUND,ISCAL,JROUND,LPRINT,
     1 NDATA,NPTOT,NPMAX,NPARM,NPFIT,JFIX,QUIT,ROBUST,
     2 IFXP(NPMAX),JFXP(MXPdim)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT), 
     1 PS(NPTOT),PSS(MXPdim),PC(MXPdim),PCS(MXPdim),PX(MXPdim),
     2 PY(MXPdim),CM(NPMAX,NPMAX), F95(10),
     3 RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, YC, Zthrd
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.MXPdim).OR.(NPTOT.GT.NDATA)
     1                                     .OR.(NPMAX.GT.MXPdim)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,602) NPTOT,MXPdim,NPMAX,NDATA
          STOP
          ENDIF
      Zthrd= 0.d0
      IF(ROBUST.GE.2) Zthrd= 1.d0/3.d0
      TSTPS= 0.d0
      RMSR= 0.d0
      NITER= 0
      QUIT= 0
      NPARM= NPTOT
      DO J= 1, NPTOT
          PS(J)= 0.d0
          JFXP(J)= IFXP(J)
          IF(IFXP(J).GT.0) NPARM= NPARM- 1
          ENDDO
      NPFIT= NPARM
      JROUND= IABS(IROUND)
c=======================================================================
c** Beginning of loop to perform rounding (if desired).  NOTE that in 
c  sequential rounding, NPARM is the current (iteratively shrinking) 
c  number of free parameters. 
    6 IF(NPARM.GT.0) TSTPS= 9.d99
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DFLOAT(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF) 
        ELSE
          TFACT= 0.D0
        ENDIF
c======================================================================
c** Begin iterative convergence loop:  try for up to  CYCMAX  cycles
      DO 50 ITER= 1, CYCMAX
          ISCAL= 0
          NITER= NITER+ 1
          DSE= 0.d0 
          TSTPSB= TSTPS
          RMSRB= RMSR
c** Zero out various arrays
   10     IF(NPARM.GT.0) THEN
              DO  I = 1,NPARM
c** PSS is the array of Saved Parameter Sensitivities from previous 
c   iteration to be carried into dyidpj subroutine - used in predicting
c   increment for derivatives by differences.
                  PSS(I)= PS(I)
c** PCS is the saved array of parameter changes from previous iteration
c   to be used (if necessary) to attempt to stablize fit
                  PCS(I)= PC(I)
                  PS(I) = 0.D0
                  PU(I) = 0.D0
                  PX(I) = 0.D0
                  PY(I) = 0.D0
                  DO  J = 1,NPARM
                      CM(I,J) = 0.D0
                      ENDDO
                  ENDDO
              ENDIF
c========Beginning of core linear least-squares step====================
c** Begin by forming the Jacobian Matrix from partial derivative matrix
          DO  I = 1,NDATA
c** User-supplied subroutine DYIDPJ uses current (trial) parameter 
c  values {PV} to generate predicted datum # I [y(calc;I)=YC] and its
c  partial derivatives w.r.t. each of the parameters, returning the 
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* NOTE 1: if more convenient, DYIDPJ could prepare the y(calc) values 
c     and derivatives for all data at the same time (when I=1), but only
c     returned the values here one datum at a time (for I > 1).]
c* NOTE 2: the partial derivative array PC returned by DYIDPJ must have
c     an entry for every parameter in the model, though for parameters 
c     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
              CALL DYIDPJ(I,NDATA,NPTOT,YC,PV,PC,PSS)
              IF(NPARM.LT.NPTOT) THEN
c** For constrained parameter or sequential rounding, collapse partial 
c   derivative array here
                  DO  J= NPTOT,1,-1
                      IF(JFXP(J).GT.0) THEN
c!! First ... move derivative for special constrained-parameter POTFIT case
cc                        IF(JFXP(J).GT.1) THEN
cc                            write(6,666) I,J,PC(J),JFXP(J),PC(JFXP(J))
cc                            PC(JFXP(J))= PC(JFXP(J))+ PC(J)
cc666 FORMAT(' For  IDAT=',I5,'  add PC(',I3,') =',1pD15.8,
cc   1  '  to PC(',0pI3,') =',1pD15.8)
cc                            ENDIF
c  ... now continue collapsing partial derivative array
                          IF(J.LT.NPTOT) THEN
                              DO  K= J,NPTOT-1
                                  PC(K)= PC(K+1)
                                  ENDDO
                              ENDIF
                          PC(NPTOT)= 0.d0
                          ENDIF
                      ENDDO
                  ENDIF
              YD(I)= YC - YO(I)
              S = 1.D0/YU(I)
cc *** For 'Robust' fitting, adjust uncertainties here
              IF(Zthrd.GT.0.d0) S= 1.d0/DSQRT(YU(I)**2 + Zthrd*YD(I)**2)
              YC= -YD(I)*S
              DSE= DSE+ YC*YC
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,YC,PX,PY)
                  ENDIF
              ENDDO
          RMSR= DSQRT(DSE/NDATA)
          IF(NPARM.LE.0) GO TO 60
c** Compute the inverse of  CM 
          CM(1,1) = 1.D0 / CM(1,1)
          DO  I = 2,NPARM
              L = I - 1
              DO  J = 1,L
                  S = 0.D0
                  DO  K = J,L
                      S = S + CM(K,I) * CM(J,K)
                      ENDDO
                  CM(J,I) = -S / CM(I,I)
                  ENDDO
              CM(I,I) = 1.D0 / CM(I,I)
              ENDDO
c** Solve for parameter changes  PC(j)
          DO  I = 1,NPARM
              J = NPARM - I + 1
              PC(J) = 0.D0
              DO  K = J,NPARM
                  PC(J) = PC(J) + CM(J,K) * PU(K)
                  ENDDO
              ENDDO
c** Get (upper triangular) "dispersion Matrix" [variance-covarience 
c  matrix  without the sigma^2 factor].
          DO  I = 1,NPARM
              DO  J = I,NPARM
                  YC = 0.D0
                  DO  K = J,NPARM
                      YC = YC + CM(I,K) * CM(J,K)
                      ENDDO
                  CM(I,J) = YC
                  ENDDO
              ENDDO
c** Generate core of Parameter Uncertainties  PU(j) and (symmetric)
c   correlation matrix  CM
          DO  J = 1,NPARM
              PU(J) = DSQRT(CM(J,J))
              DO  K= J,NPARM
                  CM(J,K)= CM(J,K)/PU(J)
                  ENDDO
              DO  K= 1,J
                  CM(K,J)= CM(K,J)/PU(J)
                  CM(J,K)= CM(K,J)
                  ENDDO
              ENDDO
c** Generate standard error  DSE = sigma^2,  and prepare to calculate 
c  Parameter Sensitivities PS
          IF(NDATA.GT.NPARM) THEN
              DSE= DSQRT(DSE/(NDATA-NPARM))
            ELSE
              DSE= 0.d0
            ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would 
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
          YC= DSE*0.1d0/DFLOAT(NPARM)
          S= DSE*TFACT
          DO  J = 1,NPARM
              PU(J)= S* PU(J)
              PS(J)= YC*DSQRT(NDATA/PS(J))
              ENDDO
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          IF((ITER.GT.1).AND.(RMSR.GT.2.0d0*RMSRB).AND.(ISCAL.LE. 3))
     1                                                            THEN
c** LeRoy's Marquardt-like scheme to damp changes if RMSR increases ...
              ISCAL= ISCAL+ 1
              IF(LPRINT.GE.0) THEN
                  WRITE(6,620) ITER,RMSR,RMSR/RMSRB,ISCAL
  620 FORMAT(' At Iteration',i3,'  RMSD=',1PD8.1,'  RMSD/RMSDB=',D8.1,
     1 "  Scale PC by  (1/4)**",i1)
ccc               WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
                  ENDIF
              DO  J= 1,NPTOT
                  PC(J)= 0.25d0*PCS(J)
                  PV(J)= PV(J)- 3.d0*PC(J)
                  ENDDO
              GOTO 10
              ENDIF
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c========End of core linear least-squares step==========================
c ... early exit if Rounding cycle finished ... 
          IF(QUIT.GT.0) GO TO 54
c
c** Next test for convergence 
          TSTPS= 0.D0
          TSTPU= 0.D0
          DO  J= 1, NPARM
              TSTPS= MAX(TSTPS,DABS(PC(J)/PS(J)))
              TSTPU= MAX(TSTPU,DABS(PC(J)/PU(J)))
              ENDDO
          IF(LPRINT.GE.4) WRITE(6,604) ITER,RMSR,TSTPS,TSTPU
c** Now ... update parameters (careful about rounding)
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN
cc               IF(JFXP(J).GT.1) THEN       !! a special PotFit option
c** If this parameter constrained to equal some earlier parameter ....
cc                   PV(J)= PV(JFXP(J))      
cc                   WRITE(6,668) J,JFXP(J),PV(J),ITER
cc                   ENDIF
cc668 FORMAT(' Constrain  PV('i3,') = PV(',I3,') =',1pd15.8,
cc   1 '  on cycle',i3)
c** If parameter held fixed (by input or rounding process), shift values
c   of change, sensitivity & uncertainty to correct label.
                  IF(J.LT.NPTOT) THEN
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      ENDIF
                  PC(J)= 0.d0
                  PS(J)= 0.d0
                  PU(J)= 0.d0
                ELSE
                  PV(J)= PV(J)+ PC(J)
                ENDIF
              ENDDO
          IF(LPRINT.GE.5) WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),
     1                                                      J=1,NPTOT)
          IF(ITER.GT.1) THEN
c** New Convergence test requires  RMSD to be constant to 1 part in 10^7
c     in adjacent cycles (unlikely to occur by accident!)
c** Replaces less severe requirement that  TSTPS < 1.0
              IF(ABS((RMSR/RMSRB)-1.d0).LT.1.d-07) THEN
                  IF(LPRINT.GE.3) WRITE(6,607) ITER,
     1                                      ABS(RMSR/RMSRB-1.d0),TSTPS
                  GO TO 54
                  ENDIF
              ENDIF
cc        CALL FLUSH(6) 
          IF(ROBUST.GT.0) Zthrd= 1.d0/3.d0
   50     CONTINUE
      WRITE(6,610) NPARM,NDATA,ITER,RMSR,TSTPS,TSTPU
c** End of iterative convergence loop for (in general) non-linear case.
c======================================================================
c
   54 IF(NPARM.LT.NPTOT) THEN
c** If necessary, redistribute correlation matrix elements to full 
c  NPTOT-element correlation matrix
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN  
c* If parameter J was held fixed
                  IF(J.LT.NPTOT) THEN
c ... then move every lower CM element down one row:
                      DO  I= NPTOT,J+1,-1
c ... For  K < J, just shift down or over to the right
                          IF(J.GT.1) THEN
                              DO  K= 1,J-1
                                  CM(I,K)= CM(I-1,K) 
                                  CM(K,I)= CM(I,K)
                                  ENDDO
                              ENDIF
c ... while for  K > J  also shift elements one column to the right
                          DO  K= NPTOT,J+1,-1
                              CM(I,K)= CM(I-1,K-1)
                              ENDDO
                          ENDDO
                      ENDIF
c ... and finally, insert appropriate row/column of zeros ....
                  DO  I= 1,NPTOT
                      CM(I,J)= 0.d0
                      CM(J,I)= 0.d0
                      ENDDO 
                  CM(J,J)= 1.d0
                  ENDIF
              ENDDO
          ENDIF
      IF(QUIT.GT.0) GOTO 60
      IF(NPARM.EQ.NPFIT) THEN
c** If desired, print unrounded parameters and fit properties
          IF(LPRINT.NE.0) THEN
              WRITE(6,616) NDATA,NPARM,RMSR,TSTPS
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      IF(IROUND.EQ.0) RETURN
c** Automated 'Sequential Rounding and Refitting' section:  round 
c  selected parameter, fix it, and return (above) to repeat fit.
      IF(IROUND.LT.0) THEN
c ... if IROUND < 0, sequentially round off 'last' remaining parameter
          DO  J= 1, NPTOT
              IF(JFXP(J).LE.0) THEN
                  JFIX= J
                  ENDIF
              ENDDO
        ELSE
c ... if IROUND > 0, sequentially round off remaining parameter with
c                    largest relative uncertainty.
c ... First, select parameter JFIX with the largest relative uncertainty
          JFIX= NPTOT
          K= 0
          TSTPS= 0.d0
          DO  J= 1,NPTOT
              IF(JFXP(J).LE.0) THEN
                  K= K+1
                  TSTPSB= DABS(PU(J)/PV(J))
                  IF(TSTPSB.GT.TSTPS) THEN
                      JFIX= J
                      TSTPS= TSTPSB
                      ENDIF
                  ENDIF
              ENDDO 
        ENDIF
      YC= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPTOT,NPTOT,JFIX,PV,PU,PS,CM)
      JFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1       WRITE(6,614) JFIX,YC,PU(JFIX),PS(JFIX),JFIX,PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all non-fixed 
c  parameters set free to get full correct final correlation matrix, 
c  uncertainties & sensitivities.  Don't update parameters on this pass!
          NPARM= NPFIT
          QUIT= 1
          DO  J= 1,NPTOT
              JFXP(J)= IFXP(J)
              ENDDO
c ... reinitialize for derivative-by-differences calculation
          RMSR= 0.d0
          ENDIF
      GO TO 6
c
c** If no parameters varied or sequential rounding completed - simply 
c   calculate DSE from RMS residuals and return.
   60 DSE= 0.d0
      IF(NDATA.GT.NPFIT) THEN
          DSE= RMSR*DSQRT(DFLOAT(NDATA)/DFLOAT(NDATA-NPFIT))
        ELSE
          DSE= 0.d0
        ENDIF
      IF(NPFIT.GT.0) THEN
          IF(LPRINT.GT.0) THEN
c** Print final rounded parameters with original Uncert. & Sensitivities
              IF(QUIT.LT.1) WRITE(6,616) NDATA, NPFIT, RMSR, TSTPS
              IF(QUIT.EQ.1) WRITE(6,616) NDATA, NPFIT, RMSR
              DO  J= 1, NPTOT
                  IF(JFXP(J).GT.0) THEN
c** If parameter held fixed (by rounding process), shift values of
c   change, sensitivity & uncertainty to correct absolute number label.
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      PC(J)= 0.d0
                      PS(J)= 0.d0
                      PU(J)= 0.d0
                      ENDIF
                  ENDDO
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      RETURN
c
  602 FORMAT(/' *** NLLSSRR problem:  [NPTOT=',i4,'] > min{MXPdim=',
     1  i4,' NPMAX=',i4,', NDATA=',i6,'}')
  604 FORMAT(' After Cycle #',i2,':  DRMSD=',1PD14.7,'    test(PS)=',
     1  1PD8.1,'   test(PU)=',D8.1)
  606 FORMAT(/' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',
     1  1PD8.1,' < 0.01   DRMSD=',D10.3)
  607 FORMAT(/' Full',i3,'-cycle convergence:  {ABS(RMSR/RMSRB)-1}=',
     1  1PD9.2,'  TSTPS=',D8.1)
  610 FORMAT(/ ' !! CAUTION !! fit of',i5,' parameters to',I6,' data not
     1 converged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((3x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d9.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/4x,'fix PV(',I4,') as ',D19.11,
     2 '  & refit:  DRMS(deviations)=',D12.5)
  616 FORMAT(/i6,' data fit to',i5,' param. yields  DRMS(devn)=',
     1 1PD14.7:'  tst(PS)=',D8.1)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.        
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF    
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.    
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED  
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS. 
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(DABS(Z(1)).LE.0.D0) GOTO 10
          IF(DABS(A(I,I)) .LT. DABS(Z(1))) THEN
              Z(2) = A(I,I) / Z(1)
              GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GC(I) = Z(2) * GS(I)
            ELSE
              Z(2) = Z(1) / A(I,I)
              GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GS(I) = Z(2) * GC(I)
            ENDIF
          A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
          Z(2) = GC(I) * F(I) + GS(I) * B
          B = GC(I) * B - GS(I) * F(I)
          F(I) = Z(2)
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ROUND(IROUND,NPMAX,NPARM,NPTOT,IPAR,PV,PU,PS,CM)
c** Subroutine to round off parameter # IPAR with value PV(IPAR) at the
c  |IROUND|'th significant digit of:  [its uncertainty  PU(IPAR)] . 
c** On return, the rounded value replaced the initial value  PV(IPAR).
c** Then ... use the correlation matrix CM and the uncertainties PU(I)
c  in the other (NPTOT-1) [or (NPARM-1) free] parameters to calculate 
c  the optimum compensating changes PV(I) in their values.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER    IROUND,NPMAX,NPARM,NPTOT,IPAR,I,IRND,KRND
      REAL*8  PU(NPMAX),PS(NPMAX),PV(NPMAX),CM(NPMAX,NPMAX),CNST,
     1        CRND,XRND,FCT,Z0
      DATA Z0/0.d0/
      CNST= PV(IPAR)
      XRND= DLOG10(PU(IPAR))
c** If appropriate, base last rounding step on sensitivity (not uncert.)
      IF((NPARM.EQ.1).AND.(PS(IPAR).LT.PU(IPAR))) XRND= DLOG10(PS(IPAR))
c** First ... fiddle with log's to perform the rounding
      IRND= INT(XRND)
      IF(XRND.GT.0) IRND=IRND+1
      IRND= IRND- IROUND
      FCT= 10.D0**IRND
      CRND= PV(IPAR)/FCT
      XRND= Z0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
      IF(DABS(CRND).GE.1.D+16) THEN
          WRITE(6,601) IROUND,IPAR
           RETURN
           ENDIF
      IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
          KRND= NINT(CRND/1.D+8)
          XRND= KRND*1.D+8
          CRND= CRND-XRND
          XRND= XRND*FCT
          END IF
      IRND= NINT(CRND)
      CNST= IRND*FCT+ XRND
c????????????????
c** Zero parameters more aggressively ... if unc. > 2* value
        if(dabs(PU(IPAR)/PV(IPAR)).GT.2.d0) then
            CNST= 0.d0
            endif
c????????????????
c** Now ... combine rounding change in parameter # IPAR, together with
c  correlation matrix CM and parameter uncertainties PU to predict
c  changes in other parameters to optimally compensate for rounding off
c  of parameter-IPAR.  Method pointed out by Mary Thompson (Dept. of
c  Statistics, UW),
      IF(IPAR.GT.1) THEN
          XRND= (CNST-PV(IPAR))/PU(IPAR)
          DO  I= 1,NPTOT
              IF(I.NE.IPAR) THEN
                  PV(I)= PV(I)+ CM(IPAR,I)*PU(I)*XRND
                  ENDIF
              ENDDO
          ENDIF
      PV(IPAR)= CNST
      RETURN
  601 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,IFXP,YC,PV,PD,PS)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i, 
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t. 
c  each of the free polynomial coefficients varied in the fit 
c  (for j=1 to NPTOT).  **  Elements of the integer array IFXP indicate
c  whether parameter j is being held fixed [IFXP(j) > 0] or varied in
c  the fit [IFXP(j).le.0].  If the former, the partial derivative 
c  for parameter j should be  PD(j)= 0.0. 
c=====================================================================
c** Use COMMON block(s) to bring in values of the independent variable 
c  [here XX(i)] and any other parameters or variables needeed to
c  calculate YC and the partial derivatives. 
c=====================================================================
c     INTEGER  I,J,NDATA,NPTOT,MXDATA,IFXP(NPTOT)
c     PARAMETER  (MXDATA= 501)
c     REAL*8  RMSR,YC,PV(NPTOT),PD(NPTOT),PS(NPTOT),POWER,XX(MXDATA)
c     COMMON /DATABLK/XX
c=====================================================================
c** NOTE BENE(!!) for non-linear fits, need to be sure that the
c  calculations of YC and PD(j) are based on the current UPDATED PV(j)
c  values.  If other (than PV) parameter labels are used internally
c  in the calculations, UPDATE them whenever (say)  I = 1 .
c=====================================================================
c     POWER= 1.D0
c     YC= PV(1)
c     PD(1)= POWER
c     DO 10 J= 2,NPTOT
c         POWER= POWER*XX(I)
c         YC= YC+ PV(J)*POWER
c         PD(J)= POWER
c  10     CONTINUE
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

