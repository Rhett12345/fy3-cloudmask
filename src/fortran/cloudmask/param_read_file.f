C--------------------------------------------------------------------
C  Copyright (C) 2002,  Space Science and Engineering Center, 
C  University C  of Wisconsin-Madison, Madison WI.
C      
C  This program is free software; you can redistribute it 
C  and/or modify it under the terms of the GNU General 
C  Public License as published by the Free Software Foundation; 
C  either version 2 of the License, or (at your option) any 
C  later version.
C
C  This program is distributed in the hope that it will be 
C  useful, but WITHOUT ANY WARRANTY; without even the implied 
C  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
C  See the  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public 
C  License along with this program; if not, write to the Free 
C  Software Foundation, Inc., 59 Temple Place, Suite 330, 
C  Boston, MA  02111-1307 USA
C--------------------------------------------------------------------
C
C
      INTEGER FUNCTION PARAM_READ_FILE( PCF_NUM, PARAM_MAX,
     &  PARAM_NUM, PARAM_LIST )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Read a parameter file. A parameter file is an ASCII text file
C     containing 1 or more name/value pairs of the form
C
C     NAME : VALUE
C
C     A valid name/value pair must contain
C     - a name containing at least one character,
C     - a colon,
C     - at least one value. More than one value
C     may be defined by using commas to separate values, e.g.
C
C     ANGLES : 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0
C
C     Comments are identified by the '!' character, which may occur
C     at the beginning of a line, or after a name/value pair, thus
C
C     ! This is a comment
C     PI : 3.1415    ! This is also a comment
C
C     are both valid comments. Blank lines are ignored. 
C
C !INPUT PARAMETERS:
C     PCF_NUM       PCF number for parameter file
C     PARAM_MAX     Maximum number of parameters
C                   (dimension of output array PARAM_LIST)
C
C !OUTPUT PARAMETERS:
C     PARAM_NUM     Number of parameters read from FILE
C     PARAM_LIST    Array of parameter strings read from FILE
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C     Original version by Liam.Gumley@ssec.wisc.edu
C
C !END
C--------------------------------------------------------------------
      
      IMPLICIT NONE
      SAVE

c --- parameters
      integer 		param_max
      character*(*) 	param_list( param_max )
      integer 		pcf_num
      integer 		param_num

c --- internal variables
      character*255 	string
      integer           lun
      integer           param_len
      integer 		count

c --- Set number of parameters found
      param_num = 0

c --- file already opened by file_open 
      lun = pcf_num

c ... Get string length of parameter list
      param_len = len( param_list( 1 ) )

c ... Check that string length of parameter list does not exceed
c ... internal string length

      if ( param_len .gt. len( string ) ) then
        param_read_file = -2
        return
      endif

c ... Read all lines from the input file, checking that maximum
c ... parameter element number is not exceeded

      count = 0
20    continue
        READ( lun, '(a)', end = 40 ) string
        count = count + 1
        if ( count .gt. param_max ) then
          param_read_file = -3
          return
        endif
        param_list( count ) = string( 1 : param_len )
      goto 20
40    continue

c --- rewind the file
      REWIND( lun )

c ... Set return values

      param_num = count
      param_read_file = 0
      
      END
