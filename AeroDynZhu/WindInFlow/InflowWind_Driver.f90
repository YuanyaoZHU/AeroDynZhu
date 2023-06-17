!****************************************************************************
!
!  PROGRAM: InflowWind_Driver  - This program tests the inflow wind module
!
!****************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************

Module InflowWind_Driver

!   USE NWTC_Library       !NOTE: Not sure why this doesn't need to be specified
   USE InflowWind
   USE InflowWind_Types
   USE InflowWind_Driver_Types    ! Contains types and routines for handling the input arguments
   USE InflowWind_Driver_Subs     ! Contains subroutines for the driver program
   
   IMPLICIT NONE
   !-----------------------------------modified by zhuyy-----------------------------------
   !====================type InflowWind_ZHU=========================================
   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("InflowWind_Driver","","")
    CHARACTER(9), PARAMETER  :: ValidParamAry(27) =  (/ &
                            "WIND1VELX","WIND1VELY","WIND1VELZ","WIND2VELX","WIND2VELY","WIND2VELZ","WIND3VELX", &
                            "WIND3VELY","WIND3VELZ","WIND4VELX","WIND4VELY","WIND4VELZ","WIND5VELX","WIND5VELY", &
                            "WIND5VELZ","WIND6VELX","WIND6VELY","WIND6VELZ","WIND7VELX","WIND7VELY","WIND7VELZ", &
                            "WIND8VELX","WIND8VELY","WIND8VELZ","WIND9VELX","WIND9VELY","WIND9VELZ"/)
     !--------------------------------------------------------------------
   TYPE, PUBLIC     :: InflowWind_TYPE_ZHU
            ! Info on this code
        
        INTEGER(IntKi)                                     :: IfWDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots

            ! Types needed here (from InflowWind module)
        TYPE(InflowWind_InitInputType)                     :: InflowWind_InitInp      ! Data for initialization -- this is where the input info goes
        TYPE(InflowWind_InputType)                         :: InflowWind_u1           ! input     -- contains xyz coords of interest -- set 1
        TYPE(InflowWind_InputType)                         :: InflowWind_u2           ! input     -- contains xyz coords of interest -- set 2
        TYPE(InflowWind_InputType)                         :: InflowWind_u3           ! input     -- contains xyz coords of interest -- set 3
        TYPE(InflowWind_ParameterType)                     :: InflowWind_p            ! Parameters
        TYPE(InflowWind_ContinuousStateType)               :: InflowWind_x            ! Continous State Data  (not used here)
        TYPE(InflowWind_DiscreteStateType)                 :: InflowWind_xd           ! Discrete State Data   (not used here)
        TYPE(InflowWind_ConstraintStateType)               :: InflowWind_z            ! Constraint State Data (not used here)
        TYPE(InflowWind_OtherStateType)                    :: InflowWind_OtherState   ! Other State Data      (Wind data is stored in here)
        TYPE(InflowWind_OutputType)                        :: InflowWind_y1           ! Output Data -- contains the velocities at xyz -- set 1
        TYPE(InflowWind_OutputType)                        :: InflowWind_y2           ! Output Data -- contains the velocities at xyz -- set 2
        TYPE(InflowWind_OutputType)                        :: InflowWind_y3           ! Output Data -- contains the velocities at xyz -- set 3
        TYPE(InflowWind_MiscVarType)                       :: InflowWind_MiscVars     ! misc/optimization data
        TYPE(InflowWind_InitOutputType)                    :: InflowWind_InitOut      ! Output Data -- contains the names and units


            ! Local variables for this code
        TYPE(IfWDriver_Flags)                              :: CLSettingsFlags         ! Flags indicating which command line arguments were specified
        TYPE(IfWDriver_Settings)                           :: CLSettings              ! Command line arguments passed in
        TYPE(IfWDriver_Flags)                              :: SettingsFlags           ! Flags indicating which settings were specified (includes CL and ipt file)
        TYPE(IfWDriver_Settings)                           :: Settings                ! Driver settings
        REAL(DbKi)                                         :: Timer(1:2)              ! Keep track of how long this takes to run
        REAL(DbKi)                                         :: TimeNow                 ! The current time
        INTEGER(IntKi)                                     :: NumTotalPoints          ! Number of points for this iteration
        LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence
        CHARACTER(11)                                      :: TmpNumString            ! Temporary string for holding a number
        INTEGER(IntKi)                                     :: ITime                   ! Generic counter for keeping track of the timestep index


    !FIXME: may want to borrow some of the type storage concepts from WAMIT2
    !FIXME: look at Waves.f90 for ideas on other things needed for this
            ! Local variables for the FFT calculations
        REAL(ReKi),    ALLOCATABLE                         :: FFTDataSetVel(:,:)      ! Velocity dataset for FFT calcs. Indices of (NumTimeSteps,3). Index 2 gives dimension U,V,W
        COMPLEX(ReKi), ALLOCATABLE                         :: FFTDataSetFrq(:,:)      ! Complex frequency information for the FFT (NumFreqs,3).  Index 2 gives dimension X,Y,Z
        REAL(ReKi),    ALLOCATABLE                         :: TimeArray(:)            ! Time array information.     (NumTimeSteps)
        REAL(ReKi),    ALLOCATABLE                         :: FreqArray(:)            ! Frequency array information (NumFreqs)




            ! Local variables for storing the arrays
        REAL(ReKi),ALLOCATABLE                             :: PointsXYZ(:,:)          !< (X,Y,Z) coordinates read from Points input file.  Velocity results stored in y2
        !INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices

            ! Temporary variables
        CHARACTER(1024)                                    :: TmpChar                 ! Temporary character variable

            ! Local Error Handling
        !INTEGER(IntKi)                                     :: ErrStat
        !CHARACTER(1024)                                    :: ErrMsg
        !INTEGER(IntKi)                                     :: ErrStatTmp
        !CHARACTER(2048)                                    :: ErrMsgTmp
        !INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp



            ! Temporary array for testing -- This is contained within a the InflowWind_Subs module. It appears here only for testing purposes

    END TYPE InflowWind_TYPE_ZHU
!-------------------------------------------------------------------------------------------------------------------------------------------
    contains
    
!===================subroutine InflowWind_int_zhu============================
!this subroutine is to initialize the inflowWind module
    
subroutine InflowWind_init_zhu( InflowWind_ZHU,startTime_zhu, position_zhu, pointNumber_zhu)

    type(InflowWind_TYPE_ZHU),      intent(inout)       :: InflowWind_ZHU    
    real(DbKi),                     intent(in   )       :: startTime_zhu
    real(ReKi),                     intent(in   )       :: position_zhu(:,:)
    integer(4),                     intent(in   )       :: pointNumber_zhu
	

    
    
        ! Local Error Handling
    INTEGER(IntKi)                                     :: ErrStat
    CHARACTER(1024)                                    :: ErrMsg
    INTEGER(IntKi)                                     :: ErrStatTmp
    CHARACTER(2048)                                    :: ErrMsgTmp
    INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp
    INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices
   !--------------------------------------------------------------------------
   !-=-=- Initialize the Library -=-=-
   !--------------------------------------------------------------------------

   CALL NWTC_Init
   CALL DispNVD(ProgInfo)

!   Beep = .FALSE.



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Setup the program -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

      ! Start the timer
   CALL CPU_TIME( InflowWind_ZHU%Timer(1) )


      ! Set some CLSettings to null/default values
   InflowWind_ZHU%CLSettings%DvrIptFileName           =  ""             ! No input name name until set
   InflowWind_ZHU%CLSettings%IfWIptFileName           =  ""             ! No IfW input file name until set
   InflowWind_ZHU%CLSettings%SummaryFileName          =  ""             ! No summary file name until set
   InflowWind_ZHU%CLSettings%NumTimeSteps             =  0_IntKi
   InflowWind_ZHU%CLSettings%DT                       =  0.0_DbKi
   InflowWind_ZHU%CLSettings%TStart                   =  0.0_ReKi
   InflowWind_ZHU%CLSettings%FFTcoord                 =  0.0_ReKi       ! Set to origin
   InflowWind_ZHU%CLSettings%GridDelta                =  0.0_ReKi       ! No stepsize
   InflowWind_ZHU%CLSettings%GridN                    =  1_IntKi        ! No grid points to calculate -- center of grid only
   InflowWind_ZHU%CLSettings%XRange                   =  0.0_ReKi       ! No xrange points
   InflowWind_ZHU%CLSettings%YRange                   =  0.0_ReKi       ! No Yrange points
   InflowWind_ZHU%CLSettings%ZRange                   =  0.0_ReKi       ! No Zrange points
   InflowWind_ZHU%CLSettings%PointsFileName           =  ""             ! No points file name until set
   InflowWind_ZHU%CLSettings%PointsOutputName         =  ""             ! No points file name until set
   InflowWind_ZHU%CLSettings%FFTOutputName            =  ""             ! No FFT output file name until set
   InflowWind_ZHU%CLSettings%WindGridOutputName       =  ""             ! No WindGrid output file name until set
   InflowWind_ZHU%CLSettings%WindGridOutputUnit       =  -1_IntKi       ! No WindGrid output unit set
   InflowWind_ZHU%CLSettings%FFTOutputUnit            =  -1_IntKi       ! No FFT output unit set
   InflowWind_ZHU%CLSettings%PointsOutputUnit         =  -1_IntKi       ! No Points file output unit set
   InflowWind_ZHU%CLSettings%ProgInfo                 =  ProgInfo       ! Driver info

      ! Set some CLSettingsFlags to null/default values
   InflowWind_ZHU%CLSettingsFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   InflowWind_ZHU%CLSettingsFlags%IfWIptFile          =  .FALSE.        ! InflowWind input filename given as command line argument
   InflowWind_ZHU%CLSettingsFlags%Summary             =  .FALSE.        ! create a summary at command line? (data extents in the wind file)
   InflowWind_ZHU%CLSettingsFlags%SummaryFile         =  .FALSE.        ! create a summary file of the output?
   InflowWind_ZHU%CLSettingsFlags%TStart              =  .FALSE.        ! specified time to start at
   InflowWind_ZHU%CLSettingsFlags%NumTimeSteps        =  .FALSE.        ! specified a number of timesteps
   InflowWind_ZHU%CLSettingsFlags%NumTimeStepsDefault =  .FALSE.        ! specified 'DEFAULT' for number of timesteps
   InflowWind_ZHU%CLSettingsFlags%DT                  =  .FALSE.        ! specified a resolution in time
   InflowWind_ZHU%CLSettingsFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   InflowWind_ZHU%CLSettingsFlags%FFTcalc             =  .FALSE.        ! do an FFT
   InflowWind_ZHU%CLSettingsFlags%WindGrid            =  .FALSE.        ! Requested output of wind data on a grid -- input file option only
   InflowWind_ZHU%CLSettingsFlags%XRange              =  .FALSE.        ! specified a range of x      -- command line option only -- stored as GridCtrCoord and GridDelta
   InflowWind_ZHU%CLSettingsFlags%YRange              =  .FALSE.        ! specified a range of y      -- command line option only -- stored as GridCtrCoord and GridDelta
   InflowWind_ZHU%CLSettingsFlags%ZRange              =  .FALSE.        ! specified a range of z      -- command line option only -- stored as GridCtrCoord and GridDelta
   InflowWind_ZHU%CLSettingsFlags%Dx                  =  .FALSE.        ! specified a resolution in x -- command line option only, 0.0 otherwise
   InflowWind_ZHU%CLSettingsFlags%Dy                  =  .FALSE.        ! speficied a resolution in y
   InflowWind_ZHU%CLSettingsFlags%Dz                  =  .FALSE.        ! specified a resolution in z
   InflowWind_ZHU%CLSettingsFlags%PointsFile          =  .FALSE.        ! points filename to read in  -- command line option only
   InflowWind_ZHU%CLSettingsFlags%WindGridOutputInit  =  .FALSE.        ! Wind Grid output file not started
   InflowWind_ZHU%CLSettingsFlags%FFTOutputInit       =  .FALSE.        ! FFT output file not started
   InflowWind_ZHU%CLSettingsFlags%PointsOutputInit    =  .FALSE.        ! Points output file not started
   InflowWind_ZHU%CLSettingsFlags%Verbose             =  .FALSE.        ! Turn on verbose error reporting?
   InflowWind_ZHU%CLSettingsFlags%VVerbose            =  .FALSE.        ! Turn on very verbose error reporting?
   InflowWind_ZHU%CLSettingsFlags%WrHAWC              =  .FALSE.        ! don't convert to HAWC format
   InflowWind_ZHU%CLSettingsFlags%WrBladed            =  .FALSE.        ! don't convert to Bladed format
   InflowWind_ZHU%CLSettingsFlags%WrVTK               =  .FALSE.        ! don't convert to VTK format

      ! Initialize the driver settings to their default values (same as the CL -- command line -- values)
   InflowWind_ZHU%Settings       =  InflowWind_ZHU%CLSettings
   InflowWind_ZHU%SettingsFlags  =  InflowWind_ZHU%CLSettingsFlags


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Parse the command line inputs -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   CALL RetrieveArgs( InflowWind_ZHU%CLSettings, InflowWind_ZHU%CLSettingsFlags, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL WrScr( NewLine//ErrMsg )
      ErrStat  =  ErrID_None
   ENDIF


      ! Check if we are doing verbose error reporting
   IF ( InflowWind_ZHU%CLSettingsFlags%VVerbose ) THEN
      InflowWind_ZHU%IfWDriver_Verbose =  10_IntKi
   ENDIF
   IF ( InflowWind_ZHU%CLSettingsFlags%Verbose ) THEN
      InflowWind_ZHU%IfWDriver_Verbose =  7_IntKi
   ENDIF



      ! Verbose error reporting
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Settings from the command line: ---')
      CALL printSettings( InflowWind_ZHU%CLSettingsFlags, InflowWind_ZHU%CLSettings )
      CALL WrSCr(NewLine)
   ENDIF


      ! Verbose error reporting
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Driver settings (before reading driver ipt file): ---')
      CALL printSettings( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings )
      CALL WrScr(NewLine)
   ENDIF




      ! Copy the input file information from the CLSettings to the Settings.
      ! At this point only one input file type can be set.
   IF ( InflowWind_ZHU%CLSettingsFlags%DvrIptFile ) THEN
      InflowWind_ZHU%SettingsFlags%DvrIptFile   =  InflowWind_ZHU%CLSettingsFlags%DvrIptFile
      InflowWind_ZHU%Settings%DvrIptFileName    =  InflowWind_ZHU%CLSettings%DvrIptFileName
   ELSE
      InflowWind_ZHU%SettingsFlags%IfWIptFile   =  InflowWind_ZHU%CLSettingsFlags%IfWIptFile
      InflowWind_ZHU%Settings%IfWIptFileName    =  InflowWind_ZHU%CLSettings%IfWIptFileName
   ENDIF


      ! If the filename given was not the IfW input file (-ifw option), then it is treated
      ! as the driver input file (flag should be set correctly by RetrieveArgs).  So, we must
      ! open this.
   IF ( InflowWind_ZHU%SettingsFlags%DvrIptFile ) THEN

         ! Read the driver input file
      CALL ReadDvrIptFile( InflowWind_ZHU%CLSettings%DvrIptFileName, InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings, ProgInfo, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF


         ! VVerbose error reporting
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after reading the driver ipt file: ---')
         CALL printSettings( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings )
         CALL WrScr(NewLine)
      ENDIF


         ! VVerbose error reporting
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('Updating driver settings with command line arguments')


         ! Now that we have read in the driver input settings, we need to override these with any
         ! values from the command line arguments.  The .TRUE. indicates that a driver input file
         ! was read.
      CALL UpdateSettingsWithCL( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings, InflowWind_ZHU%CLSettingsFlags, InflowWind_ZHU%CLSettings, .TRUE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

         ! Verbose error reporting
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrSCr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings )
         CALL WrScr(NewLine)
      ENDIF


   ELSE

      InflowWind_ZHU%SettingsFlags%WrHAWC = .FALSE.
      InflowWind_ZHU%SettingsFlags%WrBladed = .FALSE.
      InflowWind_ZHU%SettingsFlags%WrVTK = .FALSE.

         ! VVerbose error reporting
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('No driver input file used. Updating driver settings with command line arguments')


         ! Since there were no settings picked up from the driver input file, we need to copy over all
         ! the CLSettings into the regular Settings.  The .FALSE. is a flag indicating that the driver
         ! input file was not read.
      CALL UpdateSettingsWithCL( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings, InflowWind_ZHU%CLSettingsFlags, InflowWind_ZHU%CLSettings, .FALSE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

         ! Verbose error reporting
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings )
         CALL WrScr(NewLine)
      ENDIF

   ENDIF



      ! Sanity check: if an input points file is specified, make sure it actually exists. Open it if specified

   IF ( InflowWind_ZHU%SettingsFlags%PointsFile ) THEN
      INQUIRE( file=TRIM(InflowWind_ZHU%Settings%PointsFileName), exist=InflowWind_ZHU%TempFileExist )
      IF ( InflowWind_ZHU%TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find the points file "//TRIM(InflowWind_ZHU%Settings%PointsFileName))

         ! Now read the file in and save the points
      CALL ReadPointsFile( InflowWind_ZHU%Settings%PointsFileName, InflowWind_ZHU%PointsXYZ, ErrStat,ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

         ! Make name for output
      CALL GetRoot( InflowWind_ZHU%Settings%PointsFileName, InflowWind_ZHU%Settings%PointsOutputName )
      InflowWind_ZHU%Settings%PointsOutputName  =  TRIM(InflowWind_ZHU%Settings%PointsOutputName)//'.Velocity.dat'

      CALL WrScr(NewLine//"Read "//TRIM(Num2LStr(SIZE(InflowWind_ZHU%PointsXYZ,DIM=2)))//" points from '"//TRIM(InflowWind_ZHU%Settings%PointsFileName)//   &
         "'.  Results output to '"//TRIM(InflowWind_ZHU%Settings%PointsOutputName)//"'.")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(InflowWind_ZHU%Settings%PointsOutputName), exist=InflowWind_ZHU%TempFileExist )
      IF ( InflowWind_ZHU%TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(InflowWind_ZHU%Settings%PointsOutputName))

   ENDIF


      ! Setup the FFT output
   IF ( InflowWind_ZHU%SettingsFlags%FFTcalc ) THEN
         ! Make name for output
      IF ( InflowWind_ZHU%SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( InflowWind_ZHU%Settings%DvrIptFileName, InflowWind_ZHU%Settings%FFTOutputName )
      ELSE
         CALL GetRoot( InflowWind_ZHU%Settings%IfWIptFileName, InflowWind_ZHU%Settings%FFTOutputName )
      ENDIF

      InflowWind_ZHU%Settings%FFTOutputName  =  TRIM(InflowWind_ZHU%Settings%FFTOutputName)//'_'// &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(1)))//'x_'// &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(2)))//'y_'// &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(3)))//'z'//'.FFT'

      CALL WrScr(NewLine//"Writing FFT results to '"//TRIM(InflowWind_ZHU%Settings%FFTOutputName)//"' for coordinate ( "//    &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(1)))//", "// &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(2)))//", "// &
         TRIM(Num2LStr(InflowWind_ZHU%Settings%FFTcoord(3)))//" ).")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(InflowWind_ZHU%Settings%FFTOutputName), exist=InflowWind_ZHU%TempFileExist )
      IF ( InflowWind_ZHU%TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(InflowWind_ZHU%Settings%FFTOutputName))


   ENDIF


      ! Setup WindGrid output
   IF ( InflowWind_ZHU%SettingsFlags%WindGrid ) THEN

         ! Create WindGrid output name
      IF ( InflowWind_ZHU%SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( InflowWind_ZHU%Settings%DvrIptFileName, InflowWind_ZHU%Settings%WindGridOutputName )
      ELSE
         CALL GetRoot( InflowWind_ZHU%Settings%IfWIptFileName, InflowWind_ZHU%Settings%WindGridOutputName )
      ENDIF

      InflowWind_ZHU%Settings%WindGridOutputName   =  TRIM(InflowWind_ZHU%Settings%WindGridOutputName)//'.WindGrid.out'

         ! Output message if some verbosity.
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WindGridMessage(InflowWind_ZHU%Settings, .FALSE., ErrMsgTmp, LenErrMsgTmp )         ! .FALSE. for no comment characters.  ErrMsgTmp holds the message.

         CALL WrScr(NewLine//TRIM(ErrMsgTmp)//NewLine)
      ENDIF

   ENDIF


      ! Summary file output
   IF ( InflowWind_ZHU%SettingsFlags%SummaryFile ) THEN

         ! Create SummaryFile output name
      IF ( InflowWind_ZHU%SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( InflowWind_ZHU%Settings%DvrIptFileName, InflowWind_ZHU%Settings%SummaryFileName )
      ELSE
         CALL GetRoot( InflowWind_ZHU%Settings%IfWIptFileName, InflowWind_ZHU%Settings%SummaryFileName )
      ENDIF

      InflowWind_ZHU%Settings%SummaryFileName   =  TRIM(InflowWind_ZHU%Settings%SummaryFileName)//'.sum'

      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('Driver summary output file: '//TRIM(InflowWind_ZHU%Settings%SummaryFileName))

   ENDIF


      ! Give status update of the driver flags, if verbose
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after finalizing: ---')
      CALL printSettings( InflowWind_ZHU%SettingsFlags, InflowWind_ZHU%Settings )
      CALL WrScr(NewLine)
   ENDIF


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Initialize the Module -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Initialize the InflowWind module --> it will initialize all its pieces

   !DATA HANDLING:
   !
   !  During the initialization call to the IfW_Init, we specify how many points we expect to need.  For this,
   !  we will ask for the WindGrid data (if that was set), or for a single point at 50 m (random height, not reported).
   !  The arrays in the IfW_InputType and IfW_OutputType will be allocated by IfW_Init based on this number. The
   !  u1 and y1 are used for this.
   !
   !  After initialization, we will create a second set of data corresponding to the Points file data.  This set
   !  of arrays will be allocated by the driver code and copied over.  The number of points will usually be
   !  different than the number asked for in the WindGrid data. The u2 and y2 are used for this.
   !
   !  During the calculations in the driver code, we will use the MOVE_ALLOC statement to switch between the set
   !  of arrays (position and velocity) for the WindGrid data and the Points file data.  The reason for doing
   !  this is to demonstrate one method for handling AeroDynamic loading wind data (blades etc) and LIDAR requested
   !  wind information, which may be on a different timestep basis.
   !
   !  The one complication that occurs with this is that if the


      ! Set the number of points we are expecting to ask for initially
   InflowWind_ZHU%InflowWind_InitInp%NumWindPoints    =  1_IntKi                       ! Default for only asking at a single height
   !-----------------------------modified by zhuyy----------------------------------
   !IF ( SettingsFlags%WindGrid ) THEN
         ! Just in case the number of gridpoints in a direction is <= 0.
      !InflowWind_InitInp%NumWindPoints =  ( MAX(Settings%GridN(1),1_IntKi) ) * ( MAX(Settings%GridN(2),1_IntKi) ) *   &
      !                                    ( MAX(Settings%GridN(3),1_IntKi) )
   !ENDIF
    InflowWind_ZHU%InflowWind_InitInp%NumWindPoints = pointNumber_zhu !ZHU
   !---------------------------------------------------------------------------------
  

!      ! Should be able to either allocate InflowWind_u1%PositionXYZ here, or in IfW
!   CALL AllocAry( InflowWind_u1%PositionXYZ, 3, InflowWind_InitInp%NumWindPoints, &
!         "Array of positions at which to find wind velocities for the WindGrid", ErrStat, ErrMsg )



      ! Some other settings
   InflowWind_ZHU%InflowWind_InitInp%InputFileName    =  InflowWind_ZHU%Settings%IfWIptFileName       ! For now, IfW cannot work without an input file.
   !InflowWind_InitInp%DT               =  Settings%DT
   InflowWind_ZHU%InflowWind_InitInp%UseInputFile     =  .TRUE.
   IF ( InflowWind_ZHU%SettingsFlags%DvrIptFile )  THEN
      CALL GetRoot( InflowWind_ZHU%Settings%DvrIptFileName, InflowWind_ZHU%InflowWind_InitInp%RootName )
   ELSE
      CALL GetRoot( InflowWind_ZHU%InflowWind_InitInp%InputFileName, InflowWind_ZHU%InflowWind_InitInp%RootName )
      !InflowWind_InitInp%RootName = ""
   END IF
   InflowWind_ZHU%InflowWind_InitInp%RootName = trim(InflowWind_ZHU%InflowWind_InitInp%RootName)//'.IfW'
   
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) CALL WrScr('Calling InflowWind_Init...')


   CALL InflowWind_Init( InflowWind_ZHU%InflowWind_InitInp, InflowWind_ZHU%InflowWind_u1, InflowWind_ZHU%InflowWind_p, &
                  InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                  InflowWind_ZHU%InflowWind_y1, InflowWind_ZHU%InflowWind_MiscVars, InflowWind_ZHU%Settings%DT,  InflowWind_ZHU%InflowWind_InitOut, ErrStat, ErrMsg )


      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup(InflowWind_ZHU)
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                           ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF



      ! Let user know we returned from the InflowWind code if verbose
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'InflowWind_Init CALL returned without errors.'//NewLine)


      ! Convert InflowWind file to HAWC format
   IF (InflowWind_ZHU%SettingsFlags%WrHAWC) THEN
      
      CALL InflowWind_Convert2HAWC( InflowWind_ZHU%InflowWind_InitInp%RootName, InflowWind_ZHU%InflowWind_p, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg )
      
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup(InflowWind_ZHU)
            CALL ProgAbort( ErrMsg )
         ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' InflowWind_Convert2HAWC returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE
         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'InflowWind_Convert2HAWC CALL returned without errors.'//NewLine)
      END IF
      
   END IF ! Write HAWC2 files
   

      ! Convert InflowWind file to Native Bladed format
   IF (InflowWind_ZHU%SettingsFlags%WrBladed) THEN
      CALL InflowWind_Convert2Bladed( InflowWind_ZHU%InflowWind_InitInp%RootName, InflowWind_ZHU%InflowWind_p, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg )
      
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup(InflowWind_ZHU)
            CALL ProgAbort( ErrMsg )
         ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' InflowWind_Convert2Bladed returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE
         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'InflowWind_Convert2Bladed CALL returned without errors.'//NewLine)
      END IF
   END IF

   IF (InflowWind_ZHU%SettingsFlags%WrVTK) THEN
      CALL InflowWind_Convert2VTK( InflowWind_ZHU%InflowWind_InitInp%RootName, InflowWind_ZHU%InflowWind_p, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg )
      
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup(InflowWind_ZHU)
            CALL ProgAbort( ErrMsg )
         ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' InflowWind_Convert2VTK returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE
         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'InflowWind_Convert2VTK CALL returned without errors.'//NewLine)
      END IF
   
   END IF
   

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things

if (InflowWind_ZHU%SettingsFlags%WindGrid .or. InflowWind_ZHU%SettingsFlags%PointsFile .or. InflowWind_ZHU%SettingsFlags%FFTcalc) then ! we can skip all of this if we haven't asked for any output

      ! Timestep -- The timestep for the calling InflowWind_CalcOutput may need to be changed to what is in the file if the
      !  DT = DEFAULT option was used in the driver input file.  This does not need to be changed in the InflowWind_Parameters
      !  since IfW doesn't care what the timestep is.

   IF ( InflowWind_ZHU%SettingsFlags%DTDefault ) THEN
      IF ( InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%ConstantDT ) THEN
         InflowWind_ZHU%Settings%DT =  InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%DT

         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5 ) CALL WrScr(' DEFAULT requested for DT. Setting to '//TRIM(Num2LStr(InflowWind_ZHU%Settings%DT))//' as given in the wind file.')

      ELSE IF ( InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%NumTSteps <= 1 ) THEN

         InflowWind_ZHU%Settings%DT =  0.0_ReKi

         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5 ) CALL WrScr(' DEFAULT requested for DT. Setting to 0.0 since only one timestep given in the wind file.')

      ELSE
         CALL DriverCleanup(InflowWind_ZHU)
         CALL ProgAbort(NewLine//" Cannot use setting of 'DT = DEFAULT'.  The timesteps are not uniform within the data file,"// &
                                 " or there are less than 3 timesteps in the file.")
      ENDIF
   ENDIF


      ! If the number of timesteps was set to default, set it to the minimum of either number of timesteps in the file,
      ! or the largest number of complete timesteps for the given DT and time range of the file.
   IF ( InflowWind_ZHU%SettingsFlags%NumTimeStepsDefault ) THEN

         ! Set the value to the number of timesteps given in the file
      InflowWind_ZHU%Settings%NumTimeSteps   =  InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%NumTSteps

         ! Tell user what we set the number of timesteps to
      IF ( InflowWind_ZHU%IfWDriver_Verbose >= 7 ) CALL WrScr(' DEFAULT requested for NumTSteps. Setting to '//TRIM(Num2LStr(InflowWind_ZHU%Settings%NumTimeSteps))//' as given in the wind file.')

         ! Reduce to an integer number of timesteps for the given DT and time we are giving results for
      IF ( (InflowWind_ZHU%Settings%NumTimeSteps * InflowWind_ZHU%Settings%DT) > (InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%TRange(2) - InflowWind_ZHU%Settings%TStart) ) THEN

            ! Set the value.  Add 1 to get the first data point
         IF ( InflowWind_ZHU%Settings%DT > 0.0_DbKi ) THEN
            InflowWind_ZHU%Settings%NumTimeSteps = FLOOR( (InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%TRange(2) - InflowWind_ZHU%Settings%TStart)/InflowWind_ZHU%Settings%DT) + 1_IntKi
         ELSE
            InflowWind_ZHU%Settings%NumTimeSteps   =  1_IntKi
         ENDIF

            ! Tell the user that we reset things again
         IF ( InflowWind_ZHU%IfWDriver_Verbose >= 7 ) CALL WrScr(' Resetting NumTSteps to '//TRIM(Num2LStr(InflowWind_ZHU%Settings%NumTimeSteps))//   &
                                                  ' to stay within the time range given in the file (with rounding errors).')
      ENDIF

            ! Simple message if not extra verbose.
         IF ( InflowWind_ZHU%IfWDriver_Verbose == 5 ) CALL WrScr(' DEFAULT requested for NumTSteps. Setting to '//TRIM(Num2LStr(InflowWind_ZHU%Settings%NumTimeSteps))//'.')

   ENDIF



      ! WindGrid output file and points array
   IF ( InflowWind_ZHU%SettingsFlags%WindGrid ) THEN

         ! Write the header for the WindGrid output file
      CALL WindGridVel_OutputWrite( InflowWind_ZHU%Settings%WindGridOutputUnit, InflowWind_ZHU%Settings%WindGridOutputName, InflowWind_ZHU%SettingsFlags%WindGridOutputInit, &
               InflowWind_ZHU%Settings, InflowWind_ZHU%InflowWind_u1%PositionXYZ, InflowWind_ZHU%InflowWind_y1%VelocityUVW,    &
               InflowWind_ZHU%TimeNow, ErrStat, ErrMsg )

!-----------------------------------------------------MODIFIED BY ZHU--------------------------------------------
         ! Setup the actual grid points -- scan order, Y,Z,X
     ! Counter  =  1_IntKi
     ! DO I = 1,Settings%GridN(1)          ! Slowest scan over X
     !    DO K = 1,Settings%GridN(3)       ! Single Z row
     !       DO J = 1,Settings%GridN(2)    ! Step through Y fastest
     !          InflowWind_u1%PositionXYZ(1,Counter)   =  Settings%XRange(1) + Settings%GridDelta(1)*( I - 1 )
     !          InflowWind_u1%PositionXYZ(2,Counter)   =  Settings%YRange(1) + Settings%GridDelta(2)*( J - 1 )
     !          InflowWind_u1%PositionXYZ(3,Counter)   =  Settings%ZRange(1) + Settings%GridDelta(3)*( K - 1 )
     !          Counter = Counter+1
     !       ENDDO
     !    ENDDO
     ! ENDDO
     
      InflowWind_ZHU%InflowWind_u1%PositionXYZ = position_zhu
     
!---------------------------------------------------------------------------------------------------------------------

   ELSE

         ! We are going to set the WindGrid with only a single point at the RefHt so that we can calculate something
      InflowWind_ZHU%InflowWind_u1%PositionXYZ(1,1)  =  0.0_ReKi                                ! X
      InflowWind_ZHU%InflowWind_u1%PositionXYZ(2,1)  =  0.0_ReKi                                ! Y
      InflowWind_ZHU%InflowWind_u1%PositionXYZ(3,1)  =  InflowWind_ZHU%InflowWind_InitOut%WindFileInfo%RefHt   ! Z

   ENDIF



      ! If we read in a list of points from the Points input file, setup the arrays for it.
   IF ( InflowWind_ZHU%SettingsFlags%PointsFile ) THEN

         ! Move the points list
      CALL MOVE_ALLOC( InflowWind_ZHU%PointsXYZ, InflowWind_ZHU%InflowWind_u2%PositionXYZ )

         ! Allocate the array for the velocity results -- 3 x Npoints
      CALL AllocAry( InflowWind_ZHU%InflowWind_y2%VelocityUVW, 3, SIZE(InflowWind_ZHU%InflowWind_u2%PositionXYZ, DIM=2 ),    &
         "Array of velocities corresponding to Points file", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup(InflowWind_ZHU)
         CALL ProgAbort( ErrMsg )
      ENDIF

         ! WriteOutput info
      IF ( ALLOCATED(InflowWind_ZHU%InflowWind_y1%WriteOutput) ) THEN
         ALLOCATE( InflowWind_ZHU%InflowWind_y2%WriteOutput(SIZE(InflowWind_ZHU%InflowWind_y1%WriteOutput)), STAT=ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL DriverCleanup(InflowWind_ZHU)
            CALL ProgAbort( 'Error allocating InflowWind_y2%WriteOutput array.' )
         ENDIF
         InflowWind_ZHU%InflowWind_y2%WriteOutput  =  InflowWind_ZHU%InflowWind_y1%WriteOutput
      ENDIF

         ! Now create the output file.  Write header information
      CALL PointsVel_OutputWrite( InflowWind_ZHU%Settings%PointsOutputUnit, InflowWind_ZHU%Settings%PointsOutputName, InflowWind_ZHU%SettingsFlags%PointsOutputInit, &
               InflowWind_ZHU%Settings, InflowWind_ZHU%InflowWind_u2%PositionXYZ, InflowWind_ZHU%InflowWind_y2%VelocityUVW,    &
               InflowWind_ZHU%TimeNow, ErrStat, ErrMsg )

   ENDIF







      ! FFT setup
   IF ( InflowWind_ZHU%SettingsFlags%FFTcalc ) THEN

         ! Allocate arrays for passing data into and out of IfW -- 3 x 1 size
      CALL AllocAry( InflowWind_ZHU%InflowWind_u3%PositionXYZ, 3, 1,    &
         "Array for FFT position information", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup(InflowWind_ZHU)
         CALL ProgAbort( ErrMsg )
      ENDIF

      CALL AllocAry( InflowWind_ZHU%InflowWind_y3%VelocityUVW, 3, 1,    &
         "Array for FFT position information", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup(InflowWind_ZHU)
         CALL ProgAbort( ErrMsg )
      ENDIF

         ! WriteOutput info
      IF ( ALLOCATED(InflowWind_ZHU%InflowWind_y1%WriteOutput) ) THEN
         ALLOCATE( InflowWind_ZHU%InflowWind_y3%WriteOutput(SIZE(InflowWind_ZHU%InflowWind_y1%WriteOutput)), STAT=ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL DriverCleanup(InflowWind_ZHU)
            CALL ProgAbort( 'Error allocating InflowWind_y2%WriteOutput array.' )
         ENDIF
         InflowWind_ZHU%InflowWind_y3%WriteOutput  =  InflowWind_ZHU%InflowWind_y1%WriteOutput
      ENDIF



         ! Allocate storage for the FFT information
!TODO: allocate space for the FFT data storage at driver level
!        FFTDataSetVel(TSteps,3)             (TimeIdx,DimIdx)     Where DimIdx:= 1=U (velocity in X), 2=V (velocity in Y), 3=W (velocity in Z)
!        FFTDataSetFrq(Freqs,3)              (FreqIdx,DimIdx)     Where DimIdx:= 1=X, 2=Y, 3=Z
!        TimeArray(TSteps)                   (TimeIdx)            Array of values for the time -- IS THIS NECESSARY? NOT FOR FFT.
!        FreqArray(NumFreqs)                 (FreqIdx)            Array of values for the frequency
!TODO: Use the WindFileTRange and WindFileNumTSteps for the FFT
   ENDIF


      ! Report the rotation of the coordinates.
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi .AND. InflowWind_ZHU%InflowWind_p%NWindVel > 0_IntKi )   THEN
      CALL WrScr(NewLine//NewLine//'  Rotation of coordinates to prime (wind file) coordinates by rotating '//   &
                  TRIM(Num2LStr(R2D*InflowWind_ZHU%InflowWind_p%PropagationDir))// &
                  ' degrees (meteorological wind direction change) ...'//NewLine)
      CALL WrScr('          ------ WindViXYZ ---------    ----- WindViXYZprime -----')

      DO I = 1,InflowWind_ZHU%InflowWind_p%NWindVel
         ErrMsgTmp   =  ''
         ErrMsgTmp   =  '   '//TRIM(Num2LStr(I))//'  '
         ErrMsgTmp   =  TRIM(ErrMsgTmp)//'      ('//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(1,I)))//      &
                        ', '//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(2,I)))//', '//                      &
                        TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(3,I)))//')'
         ErrMsgTmp   =  ErrMsgTmp(1:40)//'('//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZprime(1,I)))//   &
                        ', '//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZprime(2,I)))//', '//             &
                        TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZprime(3,I)))//')'
         CALL WrScr(TRIM(ErrMsgTmp))
      ENDDO
      CALL WrScr(NewLine)
   ENDIF

end if

end subroutine InflowWind_init_zhu

!*****************************************************************************************************************************
!============================================InflowWind_cal_zhu===============================================================
subroutine InflowWind_cal_zhu(time,InflowWind_ZHU)

    real(DbKi),        intent(in   )   :: time !time
	type(InflowWind_TYPE_ZHU),intent(inout) :: InflowWind_ZHU
    
    ! Local Error Handling
    INTEGER(IntKi)                                     :: ErrStat
    CHARACTER(1024)                                    :: ErrMsg
    INTEGER(IntKi)                                     :: ErrStatTmp
    CHARACTER(2048)                                    :: ErrMsgTmp
    INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp
    INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices

       !--------------------------------------------------------------------------------------------------------------------------------
       !-=-=- Time stepping loop -=-=-
       !--------------------------------------------------------------------------------------------------------------------------------

    !TODO: save FFT info
    !TODO: output file with data from IfW (on time steps)
    !TODO: Summary info

       IF ( InflowWind_ZHU%IfWDriver_Verbose >= 5_IntKi )    CALL WrScr(NewLine//'Calling InflowWind_CalcOutput...'//NewLine)


       !DO ITime =  0, MAX( Settings%NumTimeSteps, 1_IntKi )

            InflowWind_ZHU%TimeNow = time   
       !TimeNow  =  Settings%TStart + Settings%DT*(ITime)
            ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'flag1::TimeNow =',InflowWind_ZHU%TimeNow
             ! Get results for WindGrid data from IfW -- WindGrid may contain only a single point at the hub if the WindGrid flag isn't set.
          CALL InflowWind_CalcOutput( InflowWind_ZHU%TimeNow,  InflowWind_ZHU%InflowWind_u1, InflowWind_ZHU%InflowWind_p, &
                      InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                      InflowWind_ZHU%InflowWind_y1, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg)



             ! Make sure no errors occured that give us reason to terminate now.
          IF ( ErrStat >= AbortErrLev ) THEN
             CALL DriverCleanup(InflowWind_ZHU)
             CALL ProgAbort( ErrMsg )
          ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) ) THEN
             CALL WrScr(NewLine//' Timestep '//TRIM(Num2LStr(InflowWind_ZHU%ITime))//   &
                        ' InflowWind_Calc returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//        &
                        NewLine//'                           ErrMsg:  '//TRIM(ErrMsg)//NewLine)
          ENDIF


             ! Write the WindGrid results to a file for this timestep
          IF ( InflowWind_ZHU%SettingsFlags%WindGrid ) THEN

             CALL WindGridVel_OutputWrite( InflowWind_ZHU%Settings%WindGridOutputUnit, InflowWind_ZHU%Settings%WindGridOutputName, InflowWind_ZHU%SettingsFlags%WindGridOutputInit, &
                      InflowWind_ZHU%Settings, InflowWind_ZHU%InflowWind_u1%PositionXYZ, InflowWind_ZHU%InflowWind_y1%VelocityUVW,    &
                      InflowWind_ZHU%TimeNow, ErrStat, ErrMsg )

          ENDIF



             ! Calculate results for the Points and export them for this timestep
          IF ( InflowWind_ZHU%SettingsFlags%PointsFile ) THEN

                ! Get results for Points data from IfW
             CALL InflowWind_CalcOutput( InflowWind_ZHU%TimeNow,  InflowWind_ZHU%InflowWind_u2, InflowWind_ZHU%InflowWind_p, &
                      InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                      InflowWind_ZHU%InflowWind_y2, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg)

                ! Output the Points results for this timestep
             CALL PointsVel_OutputWrite( InflowWind_ZHU%Settings%PointsOutputUnit, InflowWind_ZHU%Settings%PointsOutputName, InflowWind_ZHU%SettingsFlags%PointsOutputInit, &
                     InflowWind_ZHU%Settings, InflowWind_ZHU%InflowWind_u2%PositionXYZ, InflowWind_ZHU%InflowWind_y2%VelocityUVW,    &
                      InflowWind_ZHU%TimeNow, ErrStat, ErrMsg )

          ENDIF



             ! Calculate results for FFT if we are performing one
          IF ( InflowWind_ZHU%SettingsFlags%FFTcalc ) THEN

                ! Get the results from IfW
    !        CALL InflowWind_CalcOutput()

                ! Copy results over to the array for storage
    !        FFTdata(ITime,:)  =


          ENDIF


      !  output table of results for the outlist comparison and check if very verbose -- print statements are
      !  used because we don't want linewrapping.
   IF ( InflowWind_ZHU%IfWDriver_Verbose >= 10_IntKi ) THEN

      write(*,'(A)') NewLine//NewLine//'   DiskVel:  ( '//TRIM(num2lstr(InflowWind_ZHU%inflowwind_y1%diskvel(1)))//', '//             &
                  TRIM(num2lstr(InflowWind_ZHU%inflowwind_y1%diskvel(2)))//', '//TRIM(num2lstr(InflowWind_ZHU%inflowwind_y1%diskvel(3)))//' )'
      write(*,'(A)') NewLine//NewLine//'   Requested wind points and writeoutput results at last timestep (t='//       &
                  TRIM(Num2LStr(InflowWind_ZHU%TimeNow))//'):'//NewLine
      write(*,'(A)') '          ------ WindViXYZ ---------    ----- WindViUVW ---------           -- AllOuts --     '//      &
                  '------------- WriteOutput -------------'
      write(*,'(A)') ' Index,      coord,        name             Vector value                        Value         '//      &
                  ' Name        Unit   OutIndex     Value'
      DO I = 1,27
         ErrMsgTmp   =  ''
         ErrMsgTmp   =  '   '//TRIM(Num2LStr(I))//'  '
         IF ( InflowWind_ZHU%InflowWind_p%NWindVel >= (I-1)/3+1 ) THEN
            ErrMsgTmp   =  TRIM(ErrMsgTmp)//'      ('//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(1,(I-1)/3+1)))//      &
                           ', '//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(2,(I-1)/3+1)))//', '//                      &
                           TRIM(Num2LStr(InflowWind_ZHU%InflowWind_p%WindViXYZ(3,(I-1)/3+1)))//')'
            ErrMsgTmp   =  ErrMsgTmp(1:25)//ValidParamAry(I)
            ErrMsgTmp   =  ErrMsgTmp(1:40)//'('//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_MiscVars%WindViUVW(1,(I-1)/3+1)))//   &
                           ', '//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_MiscVars%WindViUVW(2,(I-1)/3+1)))//', '//             &
                           TRIM(Num2LStr(InflowWind_ZHU%InflowWind_MiscVars%WindViUVW(3,(I-1)/3+1)))//')'

         ENDIF
         ErrMsgTmp   =  ErrMsgTmp(1:73)//' '//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_MiscVars%AllOuts(I)))
         DO J  = 1, InflowWind_ZHU%InflowWind_p%NumOuts
            IF ( InflowWind_ZHU%InflowWind_p%OutParam(J)%Indx == I ) THEN
               ErrMsgTmp   =  ErrMsgTmp(1:94)//InflowWind_ZHU%InflowWind_InitOut%WriteOutputHdr(J)//'  '//     &
                              InflowWind_ZHU%InflowWind_InitOut%WriteOutputUnt(J)//TRIM(Num2LStr(J))
               ErrMsgTmp   =  ErrMsgTmp(1:126)//TRIM(Num2LStr(InflowWind_ZHU%InflowWind_y1%WriteOutput(J)))
            ENDIF 
         ENDDO
         write(*,'(A)') TRIM(ErrMsgTmp)
      ENDDO
   ENDIF
   

!end if


end subroutine InflowWind_cal_zhu

   !ENDDO    ! ITime loop
!*************************************************************************************************************!
!==========================================InflowWind_end_zhu=================================================!

subroutine InflowWind_end_zhu(InflowWind_ZHU)

    TYPE(InflowWind_TYPE_ZHU), intent(inout)  :: InflowWind_ZHU

        ! Local Error Handling
    INTEGER(IntKi)                                     :: ErrStat
    CHARACTER(1024)                                    :: ErrMsg
    INTEGER(IntKi)                                     :: ErrStatTmp
    CHARACTER(2048)                                    :: ErrMsgTmp
    INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp
    INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices
	
	
   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Calculate OtherStates -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Iff we add in some averaging / TI / mean etc, it would be in OtherStates. Right now it doesn't look like we need to do that.
   !     -- Test that here.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Output results -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------


!FFT calculations occur here.  Output to file.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- We are done, so close everything down -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'flag::1 end_inflow'
   CALL InflowWind_DestroyInitOutput( InflowWind_ZHU%InflowWind_InitOut,    ErrStat, ErrMsg )
   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'flag::2 end_inflow'
   CALL InflowWind_End( InflowWind_ZHU%InflowWind_u1, InflowWind_ZHU%InflowWind_p, &
                  InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                  InflowWind_ZHU%InflowWind_y1, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg )
   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'flag::3 end_inflow'
      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup(InflowWind_ZHU)
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (1/3) returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//' InflowWind_End call 1 of 3:    ok')
   ENDIF
   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'flag::4 end_inflow'


   CALL InflowWind_End( InflowWind_ZHU%InflowWind_u2, InflowWind_ZHU%InflowWind_p, &
                  InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                  InflowWind_ZHU%InflowWind_y2, InflowWind_ZHU%InflowWind_MiscVars,  ErrStat, ErrMsg )
   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率     write(*,*)'flag::5 end_inflow'
      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup(InflowWind_ZHU)
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (2/3) returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(' InflowWind_End call 2 of 3:    ok')
   ENDIF


   ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率     write(*,*)'flag::6 end_inflow'
   CALL InflowWind_End( InflowWind_ZHU%InflowWind_u3, InflowWind_ZHU%InflowWind_p, &
                  InflowWind_ZHU%InflowWind_x, InflowWind_ZHU%InflowWind_xd, InflowWind_ZHU%InflowWind_z, InflowWind_ZHU%InflowWind_OtherState, &
                  InflowWind_ZHU%InflowWind_y3, InflowWind_ZHU%InflowWind_MiscVars, ErrStat, ErrMsg )

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup(InflowWind_ZHU)
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (3/3)  returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                 ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( InflowWind_ZHU%IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( InflowWind_ZHU%IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(' InflowWind_End call 3 of 3:    ok')
   ENDIF



   CALL DriverCleanup(InflowWind_ZHU)
   
end subroutine InflowWind_end_zhu

!CONTAINS

   SUBROUTINE DriverCleanup(InflowWind_ZHU)
    type(InflowWind_TYPE_ZHU), intent(inout)  ::InflowWind_ZHU

      if (InflowWind_ZHU%Settings%WindGridOutputUnit  > -1_IntKi ) CLOSE( InflowWind_ZHU%Settings%WindGridOutputUnit )
      if (InflowWind_ZHU%Settings%PointsOutputUnit    > -1_IntKi ) CLOSE( InflowWind_ZHU%Settings%PointsOutputUnit )
      if (InflowWind_ZHU%Settings%FFTOutputUnit       > -1_IntKi ) CLOSE( InflowWind_ZHU%Settings%FFTOutputUnit )


         ! Find out how long this actually took
      CALL CPU_TIME( InflowWind_ZHU%Timer(2) )
      CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(InflowWind_ZHU%Timer(2)-InflowWind_ZHU%Timer(1)))//' seconds')


   END SUBROUTINE DriverCleanup


END MODULE InflowWind_Driver




