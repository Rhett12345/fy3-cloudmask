!F77-INC************************************************************
!
!!Description: Include file containing parameters which will be used
!              by a number of routines in the cloud mask processing
!              package.
!
!!Input Parameters: N/A
!
!!Output Parameters: N/A
!
!!END****************************************************************

!       Include file to set main cloud mask processing
!       parameters

!       Number of 1km pixels in a MODIS line of data
!        integer,parameter :: npixel = 1500

!       Number of lines in the processing box (context)
        integer, parameter :: nlcntx = 3

!       Number of elements in the processing box (context)
        integer, parameter :: necntx = 3

!       Number of lines per scan cube
!        integer scans_cube
!        parameter(scans_cube = 10)

!       Number of bands used in cloud mask processing
        integer, parameter :: inband = 25

!       Number of 250m visible band reflectances read
        integer, parameter :: vis_band = 4

!       Actual number of bands used in algorithm
        integer, parameter :: bands_used = 19

!       Number of bands used in sun-glint algorithm
        integer, parameter :: sg_bands_used = 5

!       Number of IR channels used
!        The number of ir channels is really 17, but one is an embedded
!       1.38 micron vis channel.
        integer, parameter :: ir_band = 6

!       Number of bands to do spatial variability tests on
        integer, parameter :: var_band = 2

!       Number of 250km fov's within a 1km footprint
        integer, parameter :: num250_per_1km = 16


!       Maximum number of elements in a MODIS line of data
!       (2 visible channels)
!        integer nx
!        parameter(nx = 4 * npixel)

!       Maximum number of lines in a scan cube of data
!       (2 visible channels)
!        integer ny
!        parameter(ny = 4 * scans_cube)
 
!       Bad data value
        real, parameter :: bad_data = -999.0

!       Bad data value
        real, parameter :: misg = -99.0

!       Metadata PSA missing value
        real, parameter :: Meta_misg = -99999.0