Include "mymodule.f03"
      program transitionDipoleMomentTest
!
!     This program is a pilot code for debugging other programs being written to
!     compute transition dipole moments between non-orthogonal Slater
!     determinants.
!
!
!     H. P. Hratchian
!     Department of Chemistry & Chemical Biology
!     Center for Chemical Computation and Theory
!     University of California, Merced
!     hhratchian@ucmerced.edu
!
!
!
!     USE Connections
!
      use myModule
!
!     Variable Declarations
!
      implicit none
      character(len=512)::matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1,GMatrixFile2
      integer::i,j,k,numCmdLineArgs,nBasis1,nBasis2,nAlpha1,nBeta1,nAlpha2,  &
        nBeta2,nBasisUse1,nBasisUse2,nOccAlpha1,nOccBeta1,nVirtAlpha1,  &
        nVirtBeta1,nOccAlpha2,nOccBeta2,nVirtAlpha2,nVirtBeta2,nDets
      integer,dimension(:,:),allocatable::determinantList
      real(kind=real64)::scfEnergy1,scfEnergy2,deltaSCFEnergy,  &
        transitionDipoleStrength,oscillatorStrength
      real(kind=real64),dimension(3)::transitionDipoleMoment
      real(kind=real64),dimension(:),allocatable::dipoleIntegrals0kX,  &
        dipoleIntegrals0kY,dipoleIntegrals0kZ,overlapIntegrals0k
      type(MQC_Variable)::SMatrixAO,dipoleAOx,dipoleAOy,dipoleAOz,  &
        MOCoefficientsAlpha1,MOCoefficientsBeta1,MOCoefficientsAlpha2,  &
        MOCoefficientsBeta2,PMatrixAlpha1,PMatrixBeta1,PMatrixTotal1,  &
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2
      type(MQC_Variable)::tmpMQC1,tmpMQC2
!
!     Format Statements
!
 1000 Format(1x,'Enter densityDifference.')
 1010 Format(1x,'Matrix File 1: ',A)
 1011 Format(1x,'Matrix File 2: ',A,/)
 1100 format(1x,'File ',i1,': nBasis=',i4,'  nElectrons=',i4,'  nAlphaEl=',i4,'  nBetaEl=',i4)
 1200 Format(/,1x,'SCF Energies',/,  &
        3x,'SCF 1:     ',F15.8,' a.u.',/,  &
        3x,'SCF 2:     ',F15.8,' a.u.',/,  &
        3x,'Delta-SCF: ',F15.8,' a.u.',2x,'=',2x,F15.8,' eV',/,  &
        36x,'=',2x,F15.8,' cm^-1',/,  &
        36x,'=',2x,F15.8,' nm')
 1500 format(/,1x,'Number of Singles Determinants: ',i10)
 5000 format(/,1x,80('-'),/,34x,'Integral Summary',/,  &
        6x,'i --> a',24x,'<0|r|i->a>',19x,'<k|Psi2>',/,1x,80('-'))
 5100 format(4x,i3,' --> ',i3,5x,3(5x,f8.5),8x,f8.5)
 5200 format(1x,80('-'),/)
 5500 format(/,1x,'Transition Dipole Stength: ',f8.5,' a.u.',  &
        3x,'Oscillator Strength: ',f8.5,' a.u.')
!
!
      write(IOut,1000)
      call setDEBUG(.true.)
!
!     Get the name of the two matrix files.
!
      numCmdLineArgs = command_argument_count()
      if(numCmdLineArgs.ne.2)  &
        call mqc_error('Wrong number of command line arguments found. Use 2 arguments.')
      call get_command_argument(1,matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      call GMatrixFile1%load(matrixFilename)
      call get_command_argument(2,matrixFilename)
      write(IOut,1011) TRIM(matrixFilename)
      call GMatrixFile2%load(matrixFilename)
!
!     Pick up the number of basis functions and electrons from each matrix file
!     and report those.
!
      nBasis1 = GMatrixFile1%getVal('nBasis')
      nBasisUse1 = GMatrixFile1%getVal('nBasisUse')
      nAlpha1 = GMatrixFile1%getVal('nAlpha')
      nBeta1  = GMatrixFile1%getVal('nBeta')
      write(iOut,1100) 1,nBasis1,nAlpha1+nBeta1,nAlpha1,nBeta1
      nBasis2 = GMatrixFile2%getVal('nBasis')
      nBasisUse2 = GMatrixFile2%getVal('nBasisUse')
      nAlpha2 = GMatrixFile2%getVal('nAlpha')
      nBeta2  = GMatrixFile2%getVal('nBeta')
      write(iOut,1100) 2,nBasis2,nAlpha2+nBeta2,nAlpha2,nBeta2
!
!     Get the SCF energies from the two jobs and report them.
!
      scfEnergy1 = GMatrixFile1%getValReal('scfEnergy')
      scfEnergy2 = GMatrixFile2%getValReal('scfEnergy')
      deltaSCFEnergy = scfEnergy2-scfEnergy1
      write(iOut,1200) scfEnergy1,scfEnergy2,deltaSCFEnergy,  &
        deltaSCFEnergy*evPHartree,deltaSCFEnergy*cmM1PHartree,  &
        deltaSCFEnergy*evPHartree*nmPev
!
!     Figure out the number of occupied and virtual MOs for each determinant.
!
      nOccAlpha1 = nAlpha1
      nVirtAlpha1 = nBasisUse1 - nOccAlpha1
      nOccBeta1 = nBeta1
      nVirtBeta1 = nBasisUse1 - nOccBeta1
      nOccAlpha2 = nAlpha2
      nVirtAlpha2 = nBasisUse2 - nOccAlpha2
      nOccBeta2 = nBeta2
      nVirtBeta2 = nBasisUse2 - nOccBeta2
!
!     Build the determinant data structure for the reference and all singles
!     substituted determinants of the bra determinant. The data structure used
!     here is a simple pair based swap list.
!
      nDets = nOccAlpha1*nVirtAlpha1+nOccBeta1*nVirtBeta1+1
      write(iOut,1500) nDets
      ALLOCATE(determinantList(2,nDets))
      determinantList = 0
      k = 1
      do i = 1,nOccAlpha1
        do j = nOccAlpha1+1,nBasisUse1
          k = k+1
          determinantList(1,k) = i
          determinantList(2,k) = j
        endDo
      endDo
      do i = 1,nOccBeta1
        do j = nOccBeta1+1,nBasisUse1
          k = k+1
          determinantList(1,k) = -i
          determinantList(2,k) = -j
        endDo
      endDo
      call mqc_print(TRANSPOSE(determinantList),iOut,header='determinant list:',blank_at_top=.true.)
!
!     Pick up the MO coefficient matrices, the overlap matrix, and the AO dipole
!     integral matrix.
!
      call GMatrixFile1%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile1%getArray('DIPOLE INTEGRALS',mqcVarOut=dipoleAOx,arraynum=1)
      call GMatrixFile1%getArray('DIPOLE INTEGRALS',mqcVarOut=dipoleAOy,arraynum=2)
      call GMatrixFile1%getArray('DIPOLE INTEGRALS',mqcVarOut=dipoleAOz,arraynum=3)
      call GMatrixFile1%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha1)
      call GMatrixFile1%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=MOCoefficientsAlpha1)
      if(GMatrixFile1%isUnrestricted()) then
        call GMatrixFile1%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta1)
        call GMatrixFile1%getArray('BETA MO COEFFICIENTS',mqcVarOut=MOCoefficientsBeta1)
      else
        PMatrixBeta1 = PMatrixAlpha1
        MOCoefficientsBeta1 = MOCoefficientsAlpha1
      endIf
      PMatrixTotal1 = PMatrixAlpha1+PMatrixBeta1
      if(DEBUG) then
        call MOCoefficientsAlpha1%print(iOut,header='MOCoefficientsAlpha1')
        call MOCoefficientsBeta1%print(iOut,header='MOCoefficientsBeta1')
      endIf
!
      call GMatrixFile2%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha2)
      call GMatrixFile2%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=MOCoefficientsAlpha2)
      if(GMatrixFile2%isUnrestricted()) then
        call GMatrixFile2%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta2)
        call GMatrixFile2%getArray('BETA MO COEFFICIENTS',mqcVarOut=MOCoefficientsBeta2)
      else
        PMatrixBeta2 = PMatrixAlpha2
        MOCoefficientsBeta2 = MOCoefficientsAlpha2
      endIf
      PMatrixTotal2 = PMatrixAlpha2+PMatrixBeta2
      if(DEBUG) then
        call MOCoefficientsAlpha2%print(iOut,header='MOCoefficientsAlpha2')
        call MOCoefficientsBeta2%print(iOut,header='MOCoefficientsBeta2')
      endIf
!
!     Form the dipole vector for each <0|r|k> where the kets |k> are the
!     substitued determinants in determinantList.
!
      ALLOCATE(dipoleIntegrals0kX(nDets),dipoleIntegrals0kY(nDets),  &
        dipoleIntegrals0kZ(nDets))
      do i = 1,nDets
        if(determinantList(1,i).eq.0.and.determinantList(2,i).eq.0) then
          dipoleIntegrals0kX(i) = contraction(PMatrixTotal1,dipoleAOx) 
          dipoleIntegrals0kY(i) = contraction(PMatrixTotal1,dipoleAOy) 
          dipoleIntegrals0kZ(i) = contraction(PMatrixTotal1,dipoleAOz) 
          call mqc_print(dipoleIntegrals0kX(i),iOut,header='HF dipole X:',blank_at_top=.true.)
          call mqc_print(dipoleIntegrals0kY(i),iOut,header='HF dipole Y:')
          call mqc_print(dipoleIntegrals0kZ(i),iOut,header='HF dipole Z:')
        elseIf(determinantList(1,i).gt.0.and.determinantList(2,i).gt.0) then
          tmpMQC1 = MOCoefficientsAlpha1%column(determinantList(1,i))
          tmpMQC2 = MOCoefficientsAlpha1%column(determinantList(2,i))
          dipoleIntegrals0kX(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOx,tmpMQC2))
          dipoleIntegrals0kY(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOy,tmpMQC2))
          dipoleIntegrals0kZ(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOz,tmpMQC2))
          if(DEBUG) then
            call mqc_print(dipoleIntegrals0kX(i),iOut,header='dipole X:',blank_at_top=.true.)
            call mqc_print(dipoleIntegrals0kY(i),iOut,header='dipole Y:')
            call mqc_print(dipoleIntegrals0kZ(i),iOut,header='dipole Z:')
          endIf
        elseIf(determinantList(1,i).lt.0.and.determinantList(2,i).lt.0) then
          tmpMQC1 = MOCoefficientsBeta1%column(-determinantList(1,i))
          tmpMQC2 = MOCoefficientsBeta1%column(-determinantList(2,i))
          dipoleIntegrals0kX(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOx,tmpMQC2))
          dipoleIntegrals0kY(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOy,tmpMQC2))
          dipoleIntegrals0kZ(i) = dot_product(tmpMQC1,  &
            MQC_Variable_MatrixVector(dipoleAOz,tmpMQC2))
          if(DEBUG) then
            call mqc_print(dipoleIntegrals0kX(i),iOut,header='dipole X:',blank_at_top=.true.)
            call mqc_print(dipoleIntegrals0kY(i),iOut,header='dipole Y:')
            call mqc_print(dipoleIntegrals0kZ(i),iOut,header='dipole Z:')
          endIf
        else
          call mqc_error('Found a detminant in determinantList that doesn''t make sense.')
        endIf
      endDo
!
!     Form the vector of overlaps <0|k>.
!
      ALLOCATE(overlapIntegrals0k(nDets))

!      overlapIntegrals0k(1) = singles_hf_overlap(8,  &
!        10,nBasis1,nOccAlpha1,nOccBeta1,  &
!        nOccAlpha2,nOccBeta2,MOCoefficientsAlpha1,MOCoefficientsBeta1,  &
!        MOCoefficientsAlpha2,MOCoefficientsBeta2,SMatrixAO)
!      write(iOut,*)' Hrant - STOP!'
!      goto 999

      do i = 1,nDets
        overlapIntegrals0k(i) = singles_hf_overlap(determinantList(1,i),  &
          determinantList(2,i),nBasis1,nOccAlpha1,nOccBeta1,  &
          nOccAlpha2,nOccBeta2,MOCoefficientsAlpha1,MOCoefficientsBeta1,  &
          MOCoefficientsAlpha2,MOCoefficientsBeta2,SMatrixAO)
      endDo
!
!     Report the integrals solved above.
!
      write(iOut,5000)
      do i = 1,nDets
        write(iOut,5100) determinantList(1,i),determinantList(2,i),  &
          dipoleIntegrals0kX(i),dipoleIntegrals0kY(i),dipoleIntegrals0kZ(i),overlapIntegrals0k(i)
      endDo
      write(iOut,5200)
!
!     Put things together to get the final transition dipole moment vector and
!     oscillator strength.
!
      transitionDipoleMoment(1) = dot_product(dipoleIntegrals0kX,overlapIntegrals0k)
      transitionDipoleMoment(2) = dot_product(dipoleIntegrals0kY,overlapIntegrals0k)
      transitionDipoleMoment(3) = dot_product(dipoleIntegrals0kZ,overlapIntegrals0k)
      call mqc_print(transitionDipoleMoment,iOut,header='Transition Dipole Moment',  &
        blank_at_top=.true.)
      transitionDipoleStrength = dot_product(transitionDipoleMoment,transitionDipoleMoment)
      oscillatorStrength = (float(2)*deltaSCFEnergy/float(3))*transitionDipoleStrength
      write(iOut,5500) transitionDipoleStrength,oscillatorStrength
!
 999  Continue
      end program transitionDipoleMomentTest
