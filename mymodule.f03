      module mymodule
!
!     USE Connections
!
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
      logical::DEBUG=.true.
!
!
!     Module Procedures...
!
      CONTAINS
!
!PROCEDURE Singles_HF_Overlap
      function singles_hf_overlap(i,a,nBasis,nOccAlpha1,nOccBeta1,  &
          nOccAlpha2,nOccBeta2,CMatrixAlpha1,CMatrixBeta1,CMatrixAlpha2,  &
          CMatrixBeta2,overlapMatrix) result(overlapValue)
!
!
      implicit none
      integer::i,a,nBasis,nOccAlpha1,nOccBeta1,nOccAlpha2,nOccBeta2
      type(MQC_Variable)::CMatrixAlpha1,CMatrixAlpha2,CMatrixBeta1,  &
        CMatrixBeta2,overlapMatrix
      real(kind=real64)::overlapValue
      real(kind=real64)::overlapAlpha,overlapBeta
      type(MQC_Variable)::tmpMatrix1,tmpMatrix2
!
 1000 format(1x,'Overlap terms: ',i4,' --> ',i4,'  OvA=',f10.7,' OvB=',f10.7,' Overlap=',f10.7)
!
!
!     Compute the alpha component.
!
      tmpMatrix1 = CMatrixAlpha1
      if(i.gt.0.and.a.gt.0) call MQC_Variable_MatrixPermuteColumns(tmpMatrix1,i,a) 
      tmpMatrix2 = MatMul(Transpose(tmpMatrix1%subMatrix([1,nBasis],[1,NOccAlpha1])),  &
        MatMul(overlapMatrix,CMatrixAlpha2%subMatrix([1,nBasis],[1,NOccAlpha2])))
      overlapAlpha = tmpMatrix2%det()
!
!     Compute the beta component.
!
      tmpMatrix1 = CMatrixBeta1
      if(i.lt.0.and.a.lt.0)  &
        call MQC_Variable_MatrixPermuteColumns(tmpMatrix1,-i,-a) 
      tmpMatrix2 = MatMul(Transpose(tmpMatrix1%subMatrix([1,nBasis],[1,NOccBeta1])),  &
        MatMul(overlapMatrix,CMatrixBeta2%subMatrix([1,nBasis],[1,NOccBeta2])))
      overlapBeta = tmpMatrix2%det()
!
      overlapValue = overlapAlpha*overlapBeta
      if(DEBUG) write(iOut,1000) i,a,overlapAlpha,overlapBeta,overlapValue
      return
      end function singles_hf_overlap


      end module mymodule
