      subroutine conf_test_2val(val,locut,hicut,power,midpt,nmval,
     *                          conflev)

      implicit none
      save

c----------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for determining the level of confidence of a
c     particular clear sky test. Input two threshold values 
c     ('midpt') with associated confidence limits ('locut', 'hicut')
c     which define 
c     a range of values for a test.  "Low" and "high" cutoffs
c     refer to low or high confidence ends of an interval and 
c     not necessarily to absolute value.  Routine calculates the
c     confidence based on an "S" function.  One may change the shape
c     of the function by changing 'power' and/or 'midpt'.
c     This version accepts 2 threshold values.
c
c!Input Parameters:
c val      current individual test value
c locut    low confidence cutoff value (less then this is 0 conf.)
c hicut    high confidence cutoff value (greather than this is
c               100% conf.)
c power    S function curve power
c midpt    midpoint of S curve (currently 50% - test straight threshold)
c nmval    1 or 2 threshold values for this test?
c
c!Output Parameters:
c conflev  calculated confidence that fov is unobstructed for this test
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!End
c----------------------------------------------------------------------

c     scalar arguments
      integer nmval
      real power,val,conflev
c     array arguments
      real locut(2),hicut(2),midpt(2)

c     local scalars
      real range,coeff,s1,c,alpha1,alpha2,
     *     beta1,beta2,gamma1,gamma2

      c = 0.0
      coeff = 2.0 ** (power - 1.0)

c     Check if testing a single threshold or a range of values.
      if(nmval .eq. 2) then

c        Range.
         gamma1 = hicut(1)
         gamma2 = hicut(2)
         alpha1 = locut(1)
         alpha2 = locut(2)
         beta1 = midpt(1)
         beta2 = midpt(2)

c        Find if interval between inner cutoffs passes test 
c        or fails.
         if( (alpha1-gamma1) .gt. 0.0) then

c           Inner region fails test.          
c           Check for value beyond function range.   
            if(val .gt. alpha1 .and. val .lt. alpha2) then

               c = 0.0

            else if(val .lt. gamma1 .or. val .gt. gamma2) then

               c = 1.0

            else if(val .le. alpha1) then

c              Value is within range of lower set of limits.
               if(val .ge. beta1) then
                  range = 2.0 * (beta1 - alpha1)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - alpha1) / range
                     c = coeff * s1**power
                  end if
               else
                  range = 2.0 * (beta1 - gamma1)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = abs(val - gamma1) / range
                     c = 1.0 - (coeff * s1**power)
                  end if
               end if

            else

c              Value is within range of upper set of limits.
               if(val .le. beta2) then
                  range = 2.0 * (beta2 - alpha2)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - alpha2) / range
                     c = coeff * s1**power
                  end if
               else
                  range = 2.0 * (beta2 - gamma2)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - gamma2) / range
                     c = 1.0 - (coeff * s1**power)
                  end if
               end if

            end if

         else

c           Inner region passes test.
c           Check for value beyond function range.   
            if(val .gt. gamma1 .and. val .lt. gamma2) then

               c = 1.0

            else if(val .lt. alpha1 .or. val .gt. alpha2) then

               c = 0.0

            else if(val .le. gamma1) then

c              Value is within range of lower set of limits.
                  if(val .le. beta1) then
                  range = 2.0 * (beta1 - alpha1)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - alpha1) / range
                     c = coeff * s1**power
                  end if
               else
                  range = abs(2.0 * (beta1 - gamma1))
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = abs((val - gamma1) / range)
                     c = 1.0 - (coeff * s1**power)
                  end if
               end if

            else

c              Value is within range of upper set of limits.
               if(val .ge. beta2) then
                  range = 2.0 * (beta2 - alpha2)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - alpha2) / range
                     c = coeff * s1**power
                  end if
               else
                  range = 2.0 * (beta2 - gamma2)
                  if(abs(range) .lt. 1.0e-12) then
                     c = 0.5
                  else
                     s1 = (val - gamma2) / range
                     c = 1.0 - (coeff * s1**power)
                  end if
               end if
    
            end if

         end if

      else 

      end if

c     Force confidence values to be between 0 and 1.
c     Also catch NaN (NaN .ne. NaN is true per IEEE 754).
      if(c .ne. c) c = 0.0
      if(c .gt. 1.0) c = 1.0
      if(c .lt. 0.0) c = 0.0

      conflev = c

      return
      end
