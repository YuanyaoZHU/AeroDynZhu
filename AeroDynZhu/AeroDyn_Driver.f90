!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
program AeroDyn_Driver

   use AeroDyn_Driver_Subs
   use InflowWind_Driver
 
   implicit none
   
   !--------------------------------   添加Dll中函数
   interface
      function gettheta(theta_dll,towerRootMotion,towerRootAngle,towerRootVel,towerRootAngleVel) bind(c,name='GETTHETA')
          use,intrinsic::iso_c_binding
          real(8)::theta_dll
          real(8)::towerRootMotion(3)
          real(8)::towerRootAngle(3)
          real(8)::towerRootVel(3)
          real(8)::towerRootAngleVel(3)          
      end function
   end interface
  !-------------------------------------------------- 
   interface
       function writeforce(time1,force_a,force_b,force_c, force_t, mom) bind(c,name='WRFO1')
           use,intrinsic::iso_c_binding
           !real(8):: force_dll
            real(8)                                      :: time1
            real(8), dimension(3,19)                     :: force_a
            real(8), dimension(3,19)                     :: force_b
            real(8), dimension(3,19)                     :: force_c
            real(8), dimension(3,10)                     :: force_t
            real(8)                                      :: mom
       end function
   end interface
   
   !--------------------------------------------------
   INTERFACE
   
      SUBROUTINE DISCON ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C, NAME='DISCON')

         USE, INTRINSIC :: ISO_C_Binding
         
         REAL(4),          INTENT(INOUT) :: avrSWAP   (*)  ! DATA 

         INTEGER(2),         INTENT(INOUT) :: aviFAIL        ! FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  ! INFILE
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(*)  ! OUTNAME (Simulation RootName)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  ! MESSAGE (Message from DLL to simulation code [ErrMsg])         
      END SUBROUTINE DISCON
   END INTERFACE   
   !--------------------------------------------------
      ! Program variables

   real(DbKi)                                     :: time                 !< Variable for storing time, in seconds 
   real(DbKi)                                     :: dT_Dvr               !< copy of DT, to make sure AD didn't change it
                                                    
   type(Dvr_SimData)                              :: DvrData              ! The data required for running the AD driver
   type(AeroDyn_Data)                             :: AD                   ! AeroDyn data 
   type(Type_Zhu)                                 :: ZHU                  ! ZHU:this varible for input and output data from modelica
   
   integer(IntKi)                                 :: iCase                ! loop counter (for driver case)
   integer(IntKi)                                 :: nt                   ! loop counter (for time step)
   integer(IntKi)                                 :: j                    ! loop counter (for array of inputs)
   integer(IntKi)                                 :: numSteps             ! number of time steps in the simulation
   integer(IntKi)                                 :: errStat              ! Status of error message
   character(ErrMsgLen)                           :: errMsg               ! Error message if ErrStat /= ErrID_None
   
   integer(IntKi)                                 :: handle               !ZHU:用于接收返回值
   real(DbKi)                                     :: InputTheta 
   real(DbKi)                                     :: Force1               !ZHU:该变量用于把force从单精度转化为双精度
   real(DbKi), dimension(3,19)                     :: bla
   real(DbKi), dimension(3,19)                     :: blb
   real(DbKi), dimension(3,19)                     :: blc
   real(DbKi), dimension(3,10)                     :: force_Tower
   integer(IntKi)                                 :: i_z
   integer(IntKi)                                 :: ii_z
   real(DbKi), dimension(3)                       :: towerRootMotion       !towerRootMotion
   real(DbKi), dimension(3)                       :: towerRootAngle        !towerRootAngle
   real(DbKi), dimension(3)                       :: towerRootVel          !towerRoot 速度
   real(DbKi), dimension(3)                       :: towerRootAngleVel     !towerRoot 角速度
   integer(IntKi)                                 :: zhu_i
   type(InflowWind_TYPE_ZHU)					  :: InflowWind_ZHU		   !DEFINE FOR INFLOWWIND
   real(DbKi)									  :: planeForce            !zhu: define to sum the force act on the normal plane
   !integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   !integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   !real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   !real                                           :: UsrTime1                                ! User CPU time for simulation initialization
   !real                                           :: UsrTime2                                ! User CPU time for simulation (without intialization)
   !real                                           :: UsrTimeDiff                             ! Difference in CPU time from start to finish of program execution
   !real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   !real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   !integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)
   logical                                        :: AD_Initialized
   !****************ZHU:DEFINE BY ZHU**************************************
   !---------------ZHU: use for discon-----------------------------------
   INTEGER(2)                                     :: aviFAIL
   CHARACTER(KIND=C_CHAR)                         :: accINFILE(3)
   CHARACTER(KIND=C_CHAR)                         :: avcOUTNAME(3)
   CHARACTER(KIND=C_CHAR)                         :: avcMSG(3)
   real(8)                                        :: moment
   !-----------------------------------------------------------------------
   real(8)                                        :: lengthBladeNode(19)                      !ZHU: Length of the aerodynamic node
   real(8)                                        :: lengthTowerNode(11)                      !ZHU: Length of tower node
   real(8)                                        :: dis2rot(19)                              !ZHU:distance of note to roter
   real(8)                                        :: momentOfNode(3,19)                         !ZHU:Every moment of node
   real(8)                                        :: momentOfBlade(3)                         !ZHU: moment of blade
   real(8)                                        :: momentOfPlate                            !ZHU: moment of plate
   real(8)                                        :: momentOfPlateFromGenerator               !ZHU: moment of plate calculate from gerator
   real(8)                                        :: errMoment                                !ZHU: error of moment between generator and plate
   real(8)                                        :: verifyForce(3)                           !ZHU: to verify the force
   !---------------------------------------------------------------------
   !*******************************zhu:WindInFlow*****************************
   !--------------------------------------------------------------------------
   
   
   
   
   !----------------------------------------------------------------------
   errStat     = ErrID_None
   errMsg      = ''
   AD_Initialized = .false.
   
   time        = 0.0 ! seconds
   open(unit=200,file='test.txt',status = 'REPLACE')
   write(200,'(7A)') 'Time',char(9),'Moment of Plate from Generate',char(9),'Moment of Plate',char(9),'ERR of Mment'
   data dis2rot /0,1.36665,4.09996,6.83325,10.24995,14.34995,18.44995,22.54995,26.64995,30.74995,34.84995,38.94995,43.04995,47.14995,51.24995,54.66665,57.39995,60.13325,61.4999/          
   data lengthBladeNode /0,2.7333,2.7333,2.7333,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,2.7333,2.7333,2.7333,2.7333/
   data lengthTowerNode /3.88,7.76,7.76,7.76,7.76,7.76,7.76,7.76,7.76,7.76,3.88/
   ! Get the current time
   !call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   !call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   
   
      ! initialize this driver:
   call Dvr_Init( DvrData, ErrStat, ErrMsg)
      call CheckError()
   
   
   
   do iCase = 1, DvrData%NumCases
      call WrScr( NewLine//'Running case '//trim(num2lstr(iCase))//' of '//trim(num2lstr(DvrData%NumCases))//'.' )
   
      
      !dT = TwoPi/DvrData%Cases(iCase)%RotSpeed / DvrData%NumSect ! sec
      
      numSteps = ceiling( DvrData%Cases(iCase)%TMax / DvrData%Cases(iCase)%dT)      
      dT_Dvr   = DvrData%Cases(iCase)%dT
      
      call WrScr ('   WndSpeed='//trim(num2lstr(DvrData%Cases(iCase)%WndSpeed))//&
               ' m/s; ShearExp='//trim(num2lstr(DvrData%Cases(iCase)%ShearExp))//&
                   '; RotSpeed='//trim(num2lstr(DvrData%Cases(iCase)%RotSpeed*RPS2RPM))//&
                  ' rpm; Pitch='//trim(num2lstr(DvrData%Cases(iCase)%Pitch*R2D))//&
                    ' deg; Yaw='//trim(num2lstr(DvrData%Cases(iCase)%Yaw*R2D))//&
                     ' deg; dT='//trim(num2lstr(DvrData%Cases(iCase)%dT))//&
                     ' s; Tmax='//trim(num2lstr(DvrData%Cases(iCase)%Tmax))//&
                 ' s; numSteps='//trim(num2lstr(numSteps)) )
      
      
         ! Set the Initialization input data for AeroDyn based on the Driver input file data, and initialize AD
         ! (this also initializes inputs to AD for first time step)
      call Init_AeroDyn(iCase, DvrData, AD, dT_Dvr, errStat, errMsg,ZHU,InflowWind_ZHU)
         call CheckError()
         AD_Initialized = .true.
         
         if (.not. EqualRealNos( dT_Dvr, DvrData%Cases(iCase)%dT ) ) then
            ErrStat = ErrID_Fatal
            ErrMsg = 'AeroDyn changed the time step for case '//trim(num2lstr(iCase))//'. Change DTAero to "default".'
            call CheckError()
         end if
                                    
      
      call Dvr_InitializeOutputFile( iCase, DvrData%Cases(iCase), DvrData%OutFileData, errStat, errMsg)
         call CheckError()
         
     
     !-------------ZHU:Init discon----------------------
      ZHU%InputTheta(2)=0.0_ReKi   !ZHU: init InputTheta(2)
      
      ZHU%pitch_In=DvrData%Cases(iCase)%pitch
      call avrSWAP_Init(ZHU%avrSWAP,dT_Dvr,ZHU%pitch_In) !ZHU:init sevo parameter 'avrSWAP'
      
      ZHU%avrSWAP(1)=0  !ZHU:first time call discon
      
      accINFILE='AB'
      avcOUTNAME='AB'
      
      call DISCON ( ZHU%avrSWAP,aviFAIL, accINFILE, avcOUTNAME, avcMSG)
      open(unit=101,FILE='bladeMotion.txt',STATUS='REPLACE') !ZHU:本通道用于输出叶片节点的位置，将在Set_AD_Inputs子函数调用
     !--------------------------------------------------    
      do nt = 1, numSteps
         
         !...............................
         ! set AD inputs for nt (and keep values at nt-1 as well)
         !...............................
         
         handle=gettheta(InputTheta,towerRootMotion,towerRootAngle, towerRootVel, towerRootAngleVel)
         !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度 write(*,*)'================================================'
         !write(*,*)'InputTheta=  ',InputTheta
         !write(*,*)'TowerRootMotion[1] = ', towerRootMotion(1)
         !write(*,*)'TowerRootMotion[2] = ', towerRootMotion(2)
         !write(*,*)'TowerRootMotion[3] = ', towerRootMotion(3)
         ZHU%towerRootMotion(1) = towerRootMotion(1)
         ZHU%towerRootMotion(2) = towerRootMotion(2)
         ZHU%towerRootMotion(3) = towerRootMotion(3)
         ZHU%towerRootAngle(1) = towerRootAngle(1)
         ZHU%towerRootAngle(2) = towerRootAngle(2)
         ZHU%towerRootAngle(3) = towerRootAngle(3)
         do zhu_i = 1, 3
             ZHU%towerRootVel(zhu_i)=towerRootVel(zhu_i)
             ZHU%towerRootAngleVel(zhu_i) = towerRootAngleVel(zhu_i)
         end do
         
         
         ZHU%InputTheta(1) = -InputTheta
         ZHU%PitchAngle = 0.0_ReKi  !ZHU:get from out
         ZHU%RotorMoment = 0.0_ReKi !ZHU:get from out
         ZHU%RotorSpeed = (ZHU%InputTheta(1)-ZHU%InputTheta(2))/dT_Dvr
         ZHU%InputTheta(2) = ZHU%InputTheta(1)
         ZHU%GenSpeed = ZHU%RotorSpeed*97
         !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*) 'GenSpeed',ZHU%GenSpeed
         !==========================================
         ZHU%avrSWAP(1)=1                !ZHU: tell discon begin to simulate
         ZHU%InTime = nt*dT_Dvr          !ZHUU: tell discon the simulation time
         call avrSWAP_Set(ZHU%avrSWAP,ZHU%InTime,ZHU%GenSpeed)    !ZHU:set the avrSWAP parameters for simulation
         
         call DISCON ( ZHU%avrSWAP,aviFAIL, accINFILE, avcOUTNAME, avcMSG)  ! call discon for control
         ZHU%GenFlag=int(ZHU%avrSWAP(7))
         ZHU%avrSWAP(4)=ZHU%avrSWAP(45)   ! set the next step pitch angle
         ZHU%avrSWAP(33)=ZHU%avrSWAP(45)  ! set the next step pitch angle
         ZHU%avrSWAP(34)=ZHU%avrSWAP(45)  ! set the next step pitch angle
         !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*) 'pitch_angle=',ZHU%avrSWAP(4)
         call Set_AD_Inputs(iCase,nt,DvrData,AD,errStat,errMsg,ZHU,InflowWind_ZHU)
            call CheckError()
   
         time = AD%InputTime(2)
         
         !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*)'time=',time
            
            ! Calculate outputs at nt - 1

         call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat, errMsg )
            call CheckError()
            
        Force1 = real(AD%y%BladeLoad(1)%Force(3,19),kind=8)
            
        !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*) 'force' , size(AD%y%BladeLoad(1)%Force)
        !write(*,*) 'force dim1' , size(AD%m%X,dim=1)
        !write(*,*) 'force dim2' , size(AD%m%X,dim=2)
            
            !do i_z = 1, 3
               do ii_z = 1, 19
                  bla(1,ii_z) = AD%m%X(ii_z,1)*lengthBladeNode(ii_z)
                  bla(2,ii_z) = AD%m%Y(ii_z,1)*lengthBladeNode(ii_z)
                  bla(3,ii_z) = 0
                  blb(1,ii_z) = AD%m%X(ii_z,2)*lengthBladeNode(ii_z)
                  blb(2,ii_z) = AD%m%Y(ii_z,2)*lengthBladeNode(ii_z)
                  blb(3,ii_z) = 0
                  blc(1,ii_z) = AD%m%X(ii_z,3)*lengthBladeNode(ii_z)
                  blc(2,ii_z) = AD%m%Y(ii_z,3)*lengthBladeNode(ii_z)
                  blc(3,ii_z) = 0 
                  !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*)'bla(1, ',ii_z,')=',bla(1,ii_z)
                  
                   !bla(i_z,ii_z) = abs(real(AD%y%BladeLoad(1)%Force(i_z,ii_z),kind=8)*lengthBladeNode(ii_z))
                   !blb(i_z,ii_z) = abs(real(AD%y%BladeLoad(2)%Force(i_z,ii_z),kind=8)*lengthBladeNode(ii_z))
                   !blc(i_z,ii_z) = abs(real(AD%y%BladeLoad(3)%Force(i_z,ii_z),kind=8)*lengthBladeNode(ii_z))
               end do
               
               do ii_z = 1,11
                   force_tower(1,ii_z) =  AD%m%X_Twr(ii_z)*lengthTowerNode(ii_z)
                   force_tower(2,ii_z) =  AD%m%Y_Twr(ii_z)*lengthTowerNode(ii_z)
                   force_tower(3,ii_z) =  0
                   !write(*,*) 'force_tower(1,',ii_z,')=',force_tower(1,ii_z)
                   !write(*,*) 'force_tower(2,',ii_z,')=',force_tower(2,ii_z)
                   !write(*,*) 'force_tower(3,',ii_z,')=',force_tower(3,ii_z)
               end do    
            !end do
			planeForce = 0
			do ii_z = 1,19
				planeForce =  planeForce+bla(1,ii_z)+blb(1,ii_z)+blc(1,ii_z)
			end do
            !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*)'planeForce =	',planeForce
            !verifyforce(1)= bla(2,9)
            !verifyforce(2)= real(AD%m%X(2,9),kind=8)
            !verifyforce(3)= real(AD%m%X(3,9),kind=8)
        !------------------ZHU: Calculate total moment--------------------    
            do ii_z = 1, 19
                momentOfNode(1,ii_z) = bla(1,ii_z)*(dis2rot(ii_z)+0.15)
                momentOfNode(2,ii_z) = blb(1,ii_z)*(dis2rot(ii_z)+0.15)
                momentOfNode(3,ii_z) = blc(1,ii_z)*(dis2rot(ii_z)+0.15)
            end do
               momentOfBlade(1)=0
               momentOfBlade(2)=0
               momentOfBlade(3)=0
            do i_z = 1,3
               do ii_z = 1,19
                momentOfBlade(i_z) = momentOfNode(i_z,ii_z) + momentOfBlade(i_z)
               end do
            end do
                momentOfPlate=0
            do i_z = 1,3
                momentOfPlate = momentOfBlade(i_z) + momentOfPlate
            end do

           
        !-----------------------------------------------------------------
         moment=real(ZHU%avrSWAP(47),kind=8)
         !if (ZHU%GenFlag==1) then
         !    ZHU%MoBeGenToSh=MissAngle*867637000.00
         !end if
         
         
         
         
         handle = writeforce(time,bla,blb,blc,force_tower,moment)!输出外力
         
           !zhu:2022年2月18日注释掉原代码，原因是降低显示频率，加快计算速度    write(*,*) 'MOMENT=', moment, '          swap(47)=', ZHU%avrSWAP(47)
           
           !write(200,*) 'Time = ',time
           momentOfPlateFromGenerator = moment*97
           
           !write(200,*) 'Moment of Plate from Generator = ',momentOfPlateFromGenerator
           !write(200,*) 'Moment of Plate= ', momentOfPlate  
           errMoment = momentOfPlateFromGenerator - momentOfPlate
           !write(200,*) 'ERR of moment = ',errMoment,'       +:阻力大于风力；-:阻力小于风力' 
           
           IF(.NOT.ALLOCATED(ZHU%output)) ALLOCATE(ZHU%output(5))
                      
           ZHU%output(1)=time
           ZHU%output(2)=momentOfPlateFromGenerator
           ZHU%output(3)=momentOfPlate
           ZHU%output(4)=errMoment
           ZHU%output(5)=verifyforce(1)
           !ZHU%output(6)=verifyforce(2)
           !ZHU%output(7)=verifyforce(3)
          
           write(200,'(F7.3,a,\)') ZHU%output(1),char(9)
           write(200,'(E20.14,a,\)')((ZHU%output(i_z),char(9)),i_z=2,5)
           !write(200,'(E20.14,a,\)') ((ZHU%output(i_z),char(9)),i_z=1,7)
           write(200,*)
           !write(200,'(E30.12)') time,momentOfPlateFromGenerator,momentOfPlate,errMoment, verifyforce(1), verifyforce(2),verifyforce(3)
           !write(200,'(F7.3,a)',ADVANCE='NO') time,char(9)
           !write(200,'(11E20.14)',ADVANCE='NO') momentOfPlateFromGenerator, char(9),momentOfPlate,char(9),errMoment,char(9),verifyforce(1),char(9),verifyforce(2),char(9),verifyforce(3)
           
           
         call Dvr_WriteOutputLine(DvrData%OutFileData, time, AD%y%WriteOutput, errStat, errMsg)
            call CheckError()
            
            
            ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt

         call AD_UpdateStates( time, nt-1, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat, errMsg )
            call CheckError()
      
                  
      end do !nt=1,numSteps
      
      close(200)
      
      call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat, errMsg )
         AD_Initialized = .false.         
         call CheckError()
         close( DvrData%OutFileData%unOutFile )
         
                     
      do j = 2, numInp
         call AD_DestroyInput (AD%u(j),  errStat, errMsg)
            call CheckError()
      end do
      

   end do !iCase = 1, DvrData%NumCases
   
   
   call Dvr_End()
   
contains
!................................   
   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError
!................................   
   subroutine Dvr_End()
   
         ! Local variables
      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      
      character(*), parameter                       :: RoutineName = 'Dvr_End'
         ! Close the output file
      if (DvrData%OutFileData%unOutFile > 0) close(DvrData%OutFileData%unOutFile)
            
      if ( AD_Initialized ) then
         call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
           
      call AD_Dvr_DestroyDvr_SimData( DvrData, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      call AD_Dvr_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
      
      
   end subroutine Dvr_End
!................................   
end program AeroDyn_Driver