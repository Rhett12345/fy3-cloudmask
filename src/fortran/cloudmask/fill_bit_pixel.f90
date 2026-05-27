subroutine fill_bit_pixel(nmtests,nbands,bad_value,bad_geo,           &
                          snglnt,desert,testbits,qa_bits,bitarray,    &
                          qa_bitarray)

!use cloudmask_data_arrays

      implicit none
      save

!----------------------------------------------------------------------
!!F77 
!
!!Description:
!     Routine for placing results of the cloud mask product and
!     qa array into a line of data values.
!
!!Input parameters:
! nc            Current processing element
! nmtests       Number of tests applied to this pixel
! nbands        Number of bands successfully read for this pixel
! bad_value     Logical value indicating band radiance or reflectance 
!               value
! bad_geo       Logical variable flagging bad lat/long data
! snglnt        Logical variables where true indicates sun glint
! desert        Logical varibable whre true indicates desert ecosystem
! testbits      four-byte integer containing bit results
! qa_bits       Byte array containing qa bit results
!
!!Output Parameters:
! bitarray      Array containing line of 48 bit test results
! qa_bitarray   Array containing line of 10 byte qa results
!
!!Revision History:
!
!!Team-unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!END
!----------------------------------------------------------------------

!      include 'global.inc'

!     scalar arguments
      integer nc,nmtests,nbands
      logical bad_value,snglnt,desert,bad_geo

!     scalar arrays
      byte testbits(6),bitarray(6),qa_bits(10),qa_bitarray(10) 

!     local scalars 
      integer i,debug,h_output

!     external routines
      external set_bit,set_qa_bit

! ... Common statement for debug purposes
!      common / bug / debug, h_output

! ... debug statement ............................................
!      if (debug .gt. 0) then
!        write(h_output,'(10x/,''Within fill_bit_line routine '',/)')
!        write(h_output,'(10x/,''Bad_value,Bad_geo, desert, sunglint= '',
!     +        4L5,/)') bad_value, bad_geo, desert, snglnt
!      endif
! .............................................................

! ... Fill in final pixel values before putting into line
!     array.  This is where the quality of the cloud mask
!     product is determined.  

! ... Now decide quality of pixel confidence based upon
! ...  number of tests used and processing path

! ... Next, how many tests and bands were used? If none
! ... set all output product bits to zero.
      if (nmtests .eq. 0 .or. nbands .eq. 0 .or. bad_geo) then
        do 100 i = 1 , 6
          testbits(i) = 0
  100   continue
! ...   and set the usable qa_bit to 0 (not useful)
        qa_bits(1) = 0

! ... If there were still some bands that were useful,
! ... then scale the qaulity accordingly
! ... (< 3 set quality bit to 4)
      else if (nmtests .lt. 3) then
        call set_bit(testbits,0)  
        call set_qa_bit(qa_bits,0)
        call set_qa_bit(qa_bits,3)

! ... (< 7 set quality bit to 6)
      else if (nmtests .lt. 7) then
        call set_bit(testbits,0) 
        call set_qa_bit(qa_bits,0)
        call set_qa_bit(qa_bits,2)
        call set_qa_bit(qa_bits,3)

!     Else set qaulity to highest value of 7
      else
         call set_bit(testbits,0)
         do i = 0 , 3
           call set_qa_bit(qa_bits,i)
         enddo
      endif

! ... Now if area is in difficult processing path region then
! ...  reduce quality to 6
      if (snglnt) then
         if (qa_bits(1) .eq. 15) qa_bits(1) = 13
      endif

! ... save bit flags for the current element in the line array
      do 200 i = 1 , 10
        if (i .le. 6) bitarray(i) = testbits(i)
        qa_bitarray(i) = qa_bits(i)
  200 continue

      return
      end
