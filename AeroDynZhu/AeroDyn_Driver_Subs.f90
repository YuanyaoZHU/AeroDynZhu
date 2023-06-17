!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD
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
module AeroDyn_Driver_Subs
   
   use AeroDyn_Driver_Types   
   use AeroDyn
   use VersionInfo
   use InflowWind_Driver
   
   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_driver', '', '' )  ! The version number of this program.
                                                    
    contains

!----------------------------------------------------------------------------------------------------------------------------------
subroutine avrSWAP_Init(avrSWAP,dT_Dvr,pitch_In)  !ZHU:Write by Zhu Yuanyao
   real(4),dimension(108)      , intent(inout) :: avrSWAP       !
   real(8)                     , intent(in   ) :: dT_Dvr        !interval time
   real(4)                     , intent(in   ) :: pitch_In
  ! real(4)                     , intent(in   ) :: pitch_Out
   !----init----
   avrSWAP(1)=-2        !avrSWAP(1)=0,First call at time zero;avrSWAP(1)=1, All subsequent timesteps; avrSWAP(1)=-1, Final call at the end of the simulation; avrSWAP(1)=2, Real Time uodate step; avrSWAP(1)=-2, init only.
   avrSWAP(2)=0         !first time at 0 time point. [s]
   avrSWAP(3)=dT_Dvr    !communication interval time [s]
   avrSWAP(4)=pitch_In  !pitch angle from last
   
   avrSWAP(5:108)=0.0
   avrSWAP(33)=pitch_In
   avrSWAP(34)=pitch_In
   
END SUBROUTINE avrSWAP_Init  
!-------------------------------------------------------------------
subroutine avrSWAP_Set(avrSWAP,InTime,GenSpeed)
    real(4),dimension(108)    , intent(inout) :: avrSWAP
    real(4)                   , intent(in   ) :: InTime
    real(4)                   , intent(in   ) :: GenSpeed
 
    
    avrSWAP(2)=InTime
    
    avrSWAP(20)=GenSpeed
    
end subroutine avrSWAP_Set
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_Init(DvrData,errStat,errMsg )

   type(Dvr_SimData),            intent(  out) :: DvrData       ! driver data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Dvr_Init'

   CHARACTER(1000)                             :: inputFile     ! String to hold the file name.
   CHARACTER(200)                              :: git_commit    ! String containing the current git commit hash

   TYPE(ProgDesc), PARAMETER                   :: version   = ProgDesc( 'AeroDyn Driver', '', '' )  ! The version number of this program.

   ErrStat = ErrID_None
   ErrMsg  = ""


   DvrData%OutFileData%unOutFile   = -1
   
   CALL NWTC_Init()
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )   
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, ErrStat2 )
   IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
      call SetErrStat(ErrID_Fatal, 'The required input file was not specified on the command line.', ErrStat, ErrMsg, RoutineName) 
      
         !bjj:  if people have compiled themselves, they should be able to figure out the file name, right?         
      IF (BITS_IN_ADDR==32) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_Win32.exe' )
      ELSEIF( BITS_IN_ADDR == 64) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_x64.exe' )
      ELSE
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver.exe' )
      END IF
         
      return
   END IF        
         
      ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      if (errStat >= AbortErrLev) return
      
      ! validate the inputs
   call ValidateInputs(DvrData, errStat2, errMsg2)      
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      
end subroutine Dvr_Init 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AeroDyn(iCase, DvrData, AD, dt, errStat, errMsg, ZHU, InflowWind_ZHU)

   integer(IntKi),               intent(in   ) :: iCase         ! driver case
   type(Dvr_SimData),            intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
         
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   type(Type_Zhu)              , intent(  out) :: ZHU            !ZHU: use for Set_AD_Inputs 
   type(InflowWind_TYPE_ZHU)   , intent(inout) :: InflowWind_ZHU  !ZHU:

      ! locals
   real(reKi)                                  :: theta(3)
   integer(IntKi)                              :: j, k   
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Init_AeroDyn'
   integer(IntKi)                              :: NumNode_zhu   ! zhu:total number of turbine node
   integer(IntKi)                              :: counter_zhu = 1 ! zhu: count the number of turbine node
                                                  
   ! local data                                
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
   real(DbKi)                                  :: time = 0.0
      
  
      
   errStat = ErrID_None
   errMsg  = ''
   
   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBlades
   InitInData%RootName       = DvrData%outFileData%Root
   InitInData%Gravity        = 9.80665_ReKi                
   
      ! set initialization data:
   call AllocAry( InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
      
   InitInData%HubPosition = (/ DvrData%Overhang * cos(DvrData%shftTilt), 0.0_ReKi, DvrData%HubHt /)
   theta(1) = 0.0_ReKi
   theta(2) = -DvrData%shftTilt
   theta(3) = 0.0_ReKi
   InitInData%HubOrientation = EulerConstruct( theta )
     
   
   do k=1,InitInData%numBlades
                     
      theta(1) = (k-1)*TwoPi/real(InitInData%numBlades,ReKi)
      theta(2) = DvrData%precone
      theta(3) = 0.0_ReKi
      InitInData%BladeRootOrientation(:,:,k) = matmul( EulerConstruct( theta ), InitInData%HubOrientation )
                  
      InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + DvrData%hubRad * InitInData%BladeRootOrientation(3,:,k)      
      
   end do
      
      
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
 !-----------------------------------------------------------------------------------------     
   call AllocAry( ZHU%towerNodePosition, 3, AD%p%NumTwrNds, 'TowerNodeOriginPosition', errStat2, ErrMsg2 ) !ZHU:save the tower origion node position
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )                                  !ZHU:
      
   do k=1,AD%p%NumTwrNds
        ZHU%towerNodePosition(:,k)=AD%u(1)%TowerMotion%Position(:,k)
        ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*)'Zhu%position[',k,'] =',ZHU%towerNodePosition(1,k),',',ZHU%towerNodePosition(2,k),',',ZHU%towerNodePosition(3,k)
   end do
  !-----------------------------------------------------------------------------------------
   ZHU%NumNode = AD%p%NumTwrNds + AD%p%NumBlades * AD%p%NumBlNds
   call AllocAry(ZHU%NodePosition,3, ZHU%NumNode, 'TotalNodePosition',errStat2, ErrMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName)
      
   do k = 1, AD%p%NumTwrNds
       ZHU%NodePosition(:,counter_zhu) = AD%u(1)%TowerMotion%Position(:,k)
       counter_zhu = counter_zhu + 1
   end do
   do k = 1,AD%p%NumBlades
        do j = 1, AD%p%NumBLNds
            ZHU%NodePosition(:,counter_zhu) = AD%u(1)%BladeMotion(k)%Position(:,j)
            counter_zhu = counter_zhu + 1
        end do
   end do
   
   call InflowWind_init_zhu( InflowWind_ZHU,time, ZHU%NodePosition, ZHU%NumNode)
  !-----------------------------------------------------------------------------------------
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if   
         
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end do
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   
      ! we know exact values, so we're going to initialize inputs this way (instead of using the input guesses from AD_Init)
   AD%InputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(iCase,j,DvrData,AD,errStat2,errMsg2,ZHU,InflowWind_ZHU)   !ZHU:Add ZHU variable
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO              
   
      
      ! move AD initOut data to AD Driver
   call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
   call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )   
     
   DvrData%OutFileData%AD_ver = InitOutData%ver
   
contains
   subroutine cleanup()
      call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup
   
end subroutine Init_AeroDyn
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine returns time=(nt-1) * DvrData%Cases(iCase)%dT, and cycles values in the input array AD%InputTime and AD%u.
!! it then sets the inputs for nt * DvrData%Cases(iCase)%dT, which are index values 1 in the arrays.
subroutine Set_AD_Inputs(iCase,nt,DvrData,AD,errStat,errMsg,ZHU,InflowWind_ZHU)

    integer(IntKi)              , intent(in   ) :: iCase         ! case number 
    integer(IntKi)              , intent(in   ) :: nt            ! time step number
   
    type(Dvr_SimData),            intent(inout) :: DvrData       ! Driver data 
    type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
    integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
    character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   
    type(Type_Zhu)              , intent(inout) :: ZHU           ! ZHU:get variable from modelica
	type(InflowWind_TYPE_ZHU)   , intent(inout) :: InflowWind_ZHU
   
        ! local variables
    integer(IntKi)                              :: errStat2      ! local status of error message
    character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
    character(*), parameter                     :: RoutineName = 'Set_AD_Inputs'

    integer(intKi)                              :: j             ! loop counter for nodes
    integer(intKi)                              :: k             ! loop counter for blades
    integer(intKi)                              :: zhu_i         ! loop index for zhu

    real(ReKi)                                  :: z             ! height (m)
    !real(ReKi)                                  :: angle
    real(R8Ki)                                  :: theta(3)
    real(R8Ki)                                  :: position(3)
    real(R8Ki)                                  :: position_zhu  ! zhu
    real(R8Ki)                                  :: orientation(3,3)
    real(R8Ki)                                  :: orientation2(3,3)
    real(R8Ki)                                  :: rotateMat(3,3)
    real(R8Ki)                                  :: cross_production(3) !zhu
    integer(4)                                  :: counter_zhu=1
	
   
   
    errStat = ErrID_None
    errMsg  = ""
   
      ! note that this initialization is a little different than the general algorithm in FAST because here
      ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
      
    !................
    ! shift previous calculations:
    !................
    do j = numInp-1,1,-1
        call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
        AD%InputTime(j+1) = AD%InputTime(j)
    end do
    AD%inputTime(1) = nt * DvrData%Cases(iCase)%dT
         
    !................
    ! calculate new values
    !................
   
    ! Tower motions:
    do j=1,AD%u(1)%TowerMotion%nnodes
        !---------------------modified by zhu---------------------------------------!
        theta(1) = ZHU%towerRootAngle(1)
        theta(2) = ZHU%towerRootAngle(2)
        theta(3) = ZHU%towerRootAngle(3)

        orientation = EulerConstruct(theta)!zhu: 此orientation为塔基处的方位矩阵，2023年6月2日。
        AD%u(1)%TowerMotion%Orientation(:,:,j)=matmul( AD%u(1)%TowerMotion%RefOrientation(:,:,j) , orientation   )!zhu 被旋转的向量或原方位矩阵在前，在原基础上增加的角度构成的旋转矩阵（方位矩阵）在后。
        !-----------------------------------------------------------------------------!
        !AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j) ! identity
          
        ZHU%towerNodePosition(:,j) = ZHU%towerRootMotion(:) + matmul(AD%u(1)%TowerMotion%Position(:,j) ,AD%u(1)%TowerMotion%Orientation(:,:,j) )!zhu:这一步实现了塔筒各节点的实时位置
         
        AD%u(1)%TowerMotion%TranslationDisp(:,j) =ZHU%towerNodePosition(:,j) - AD%u(1)%TowerMotion%Position(:,j)                    !zhu:计算各点相对初始状态的位移
        !AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
        cross_production(:) = cross_product( ZHU%towerRootAngleVel(:),ZHU%towerNodePosition(:,j) )
        AD%u(1)%TowerMotion%TranslationVel( :,j) = ZHU%towerRootVel(:) +   cross_production(:)  !zhu:计算塔筒上各点在全局坐标系中的速度
         
        !do zhu_i = 1,3
        !write(*,*)'towerAngleVel[',zhu_i,'] = ', AD%u(1)%TowerMotion%TranslationVel( zhu_i,j)
        !write(*,*)'cross_production[',zhu_i,'] = ', cross_production(zhu_i)
        !end do
        !AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
    end do !j=nnodes
      
    ! Hub motions:
    theta(1) = 0.0_ReKi
    theta(2) = 0.0_ReKi
    theta(3) = DvrData%Cases(iCase)%Yaw
    orientation2 = EulerConstruct(theta)    !zhu:先局部，后整体
    orientation = matmul( orientation2,orientation) !zhu：此orientation为经过偏航后的orientation,2023年6月2日。
       
    AD%u(1)%HubMotion%TranslationDisp(:,1) = ZHU%towerRootMotion(:) + matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1)
    !AD%u(1)%HubMotion%TranslationDisp(:,1) = matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )
      
    !theta(1) = AD%inputTime(1) * DvrData%Cases(iCase)%RotSpeed                                                       !ZHU:instead by ZHU%InputTheta
    theta(1) = ZHU%InputTheta(1)
    theta(2) = 0.0_ReKi
    theta(3) = 0.0_ReKi
    AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%HubMotion%RefOrientation(:,:,1), orientation )
    orientation2 = EulerConstruct( theta )
    write(*,*)'theta = ',theta(1)
    write(*,*)'RotorSpeed = ', ZHU%RotorSpeed
    AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation2, AD%u(1)%HubMotion%Orientation(  :,:,1) )!zhu:这里之前写错了，把orientation2写成了orinetation导致实际风机叶片并没有转动
      
    !AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * DvrData%Cases(iCase)%RotSpeed    !ZHU:instead by ZHU%RotorSpeed
    AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * ZHU%RotorSpeed !这一行的计算非常具有迷惑性，实际上应该是RT*w，其中RT为Orientation的转置，w为包含RotorSpeed的速度的向量，但是按照它这种写法看似不对，实际上是和RT*w等效的。
                                                                                                    ! 这个AD%u(1)%HubMotion%RotationVel(    :,1)还应加上ZHU%towerRootVel(:)，原因可看《FAST及多体系统中的坐标转换问题（方位矩阵）》，2023年6月3日
    ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*) 'DvrRotSpeed=  ',DvrData%Cases(iCase)%RotSpeed
    ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率    write(*,*) 'ZHU%RotorSpeed=  ',Zhu%RotorSpeed
      
    ! Blade root motions:
    do k=1,DvrData%numBlades !此步将首先构建方位矩阵orientation2, 再通过与hub的方位矩阵相乘，获得叶根处的方位矩阵，叶片的局部坐标乘以         
        theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
        theta(2) =  DvrData%precone
        theta(3) = -ZHU%avrSWAP(45)          !ZHU:Instead DvrData%Cases(iCase)%pitch
       ! zhu:2022年2月18日修改，原因是降低输出频率提高计算效率     write(*,*) 'theta(3)=',theta(3)
        !theta(3) = -DvrData%Cases(iCase)%pitch
        orientation2 = EulerConstruct(theta)  !zhu:这里吧orientation改成了orientation2
         
        AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation2, AD%u(1)%HubMotion%Orientation(  :,:,1) )!zhu: 这里吧orientation改成了orientation2
         
    end do !k=numBlades
            
    ! Blade motions:
    do k=1,DvrData%numBlades
        rotateMat = transpose( AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) )!方位矩阵的转置等于其逆矩阵
        rotateMat = matmul( rotateMat, AD%u(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) ) 
        orientation = transpose(rotateMat)
         
        rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi !为什么要减去1，请查看语雀文档《FAST及多体系统中的坐标转换问题（方位矩阵）》，注：2023.5.21
        rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
        rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi
                  
        do j=1,AD%u(1)%BladeMotion(k)%nnodes        
            position = AD%u(1)%BladeMotion(k)%Position(:,j) - AD%u(1)%HubMotion%Position(:,1)
            
            AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = AD%u(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )   !这里rotateMat在它定义的时候就已经相对BladeRootMotion的Orientation修正过了，所以这里叶片的位置不用修改
            ! zhu:这一步先计算节点
            AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )
            
            
            position =  AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) &  !AD%u(1)%BladeMotion(k)%Position(:,j) 叶片节点原始位置，AD%u(1)%BladeMotion(k)%TranslationDisp(:,j)为叶片节点位移
                        - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)            !AD%u(1)%HubMotion%Position(:,1)为轮毂原始位置，AD%u(1)%HubMotion%TranslationDisp(:,1)为轮毂位移
                                                                                                              !position是叶片节点到轮毂在全局坐标系下的当前时刻的相对位置。
            !position_zhu = sqrt(position(1)**2+position(2)**2+position(3)**2) !zhu 
            !write(*,*)'position[',j,'] = ',position_zhu
            AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position ) + &                                                                   !zhu:局部运动；2023年6月3日，这一行中HubMotion%RotationVel仅考虑了轮毂旋转角速度，并未考虑塔基处的角速度影响，所以并不全面，具体修改方式可看上面计算HubMotion%RotationVel代码的注释，也可看语雀《FAST及多体系统中的坐标转换问题（方位矩阵）》
                                                            cross_product( ZHU%towerRootAngleVel(:), AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j)) +  & !zhu:局部坐标系相对全局坐标系的转动引起的速度 ；2023年6月3日，这一项应该是：cross_product( ZHU%towerRootAngleVel(:), matmul( AD%u(1)%HubMotion%Position(:,1), AD%u(1)%HubMotion%Orientation(  :,:,1)))，具体可看语雀《FAST及多体系统中的坐标转换问题（方位矩阵）》中的具体介绍
                                                            ZHU%towerRootVel(:)                                                                                                               !zhu:局部坐标系相对全局坐标系的运动
            !AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )
        end do !j=nnodes
                                    
    end do !k=numBlades
      
    !---------------------zhu: reset the position of the node-------------------------
    counter_zhu = 1                                   !初始化counter_zhu, 虽然前面已经初始化了，但是不知道为什么这里还是会延续上一步的counter_zhu,导致ZHU%NodePosition传输坐标位置发生错误。 2022年6月1日
    do k = 1, AD%p%NumTwrNds
        ZHU%NodePosition(:,counter_zhu) = AD%u(1)%TowerMotion%TranslationDisp(:,k)+AD%u(1)%TowerMotion%position(:,k)
        
        !if (k==1) then                                                          !输出到bladeMotion.txt来显示counter_zhu的数值
        !    write(101,*)'counter_zhu = ',counter_zhu    
        !end if
        
        counter_zhu = counter_zhu + 1
    end do
    do k = 1,AD%p%NumBlades
        do j = 1, AD%p%NumBLNds
            ZHU%NodePosition(:,counter_zhu) = AD%u(1)%BladeMotion(k)%TranslationDisp(:,j)+AD%u(1)%BladeMotion(k)%position(:,j)
            
            !if ( j ==  AD%p%NumBLNds) then                                      !输出到bladeMotion.txt来显示ZHU%NodePosition的位置
            !    write(101,*)ZHU%NodePosition(1,counter_zhu),'    ',ZHU%NodePosition(2,counter_zhu),'    ',ZHU%NodePosition(3,counter_zhu),'    ',counter_zhu
            !end if
            
            counter_zhu = counter_zhu + 1

            
        end do
    end do
    !write(101,*)AD%u(1)%BladeMotion(1)%TranslationDisp(1,19),'    ',AD%u(1)%BladeMotion(1)%TranslationDisp(2,19),'    ',AD%u(1)%BladeMotion(1)%TranslationDisp(3,19),'    ',AD%u(1)%BladeMotion(2)%TranslationDisp(1,19),'    ',AD%u(1)%BladeMotion(2)%TranslationDisp(2,19),'    ',AD%u(1)%BladeMotion(2)%TranslationDisp(3,19),'    ',AD%u(1)%BladeMotion(3)%TranslationDisp(1,19),'    ',AD%u(1)%BladeMotion(3)%TranslationDisp(2,19),'    ',AD%u(1)%BladeMotion(3)%TranslationDisp(3,19)
    
      
    InflowWind_ZHU%InflowWind_u1%PositionXYZ = ZHU%NodePosition
    !write(101,*)ZHU%NodePosition(1,98),'    ',ZHU%NodePosition(2,98),'    ',ZHU%NodePosition(3,98),'    ',ZHU%NodePosition(1,117),'    ',ZHU%NodePosition(2,117),'    ',ZHU%NodePosition(3,117),'    ',ZHU%NodePosition(1,136),'    ',ZHU%NodePosition(2,136),'    ',ZHU%NodePosition(3,136),'    ',AD%p%NumTwrNds
    !write(101,*)ZHU%NodePosition(1,AD%p%NumTwrNds+AD%p%NumBLNds*1),'    ',ZHU%NodePosition(2,AD%p%NumTwrNds+AD%p%NumBLNds*1),'    ',ZHU%NodePosition(3,AD%p%NumTwrNds+AD%p%NumBLNds*1),'    ',ZHU%NodePosition(1,AD%p%NumTwrNds+AD%p%NumBLNds*2),'    ',ZHU%NodePosition(2,AD%p%NumTwrNds+AD%p%NumBLNds*2),'    ',ZHU%NodePosition(3,AD%p%NumTwrNds+AD%p%NumBLNds*2),'    ',ZHU%NodePosition(1,AD%p%NumTwrNds+AD%p%NumBLNds*3),'    ',ZHU%NodePosition(2,AD%p%NumTwrNds+AD%p%NumBLNds*3),'    ',ZHU%NodePosition(3,AD%p%NumTwrNds+AD%p%NumBLNds*3)
  
    !---------------------zhu: set the wind velocity----------------------------------
    call InflowWind_cal_zhu(AD%inputTime(1),InflowWind_ZHU)
    !InflowOnTower_zhu
    counter_zhu = 1
    do j=1,AD%u(1)%TowerMotion%nnodes
        !AD%u(1)%InflowOnTower(1,j)=0
        !AD%u(1)%InflowOnTower(2,j)=0
        !AD%u(1)%InflowOnTower(3,j)=0
        AD%u(1)%InflowOnTower(:,j) = InflowWind_ZHU%InflowWind_y1%VelocityUVW(:,counter_zhu)
        counter_zhu = counter_zhu + 1
    end do
    !InflowOnBlade
    do k=1,DvrData%numBlades
        do j=1,AD%u(1)%BladeMotion(k)%nnodes
            !AD%u(1)%InflowOnBlade(1,j,k)=0
            !AD%u(1)%InflowOnBlade(2,j,k)=0
            !AD%u(1)%InflowOnBlade(3,j,k)=0
            AD%u(1)%InflowOnBlade(:,j,k) = InflowWind_ZHU%InflowWind_y1%VelocityUVW(:,counter_zhu)
            counter_zhu = counter_zhu + 1
        end do
    end do  
      
    ! Inflow wind velocities:
    ! InflowOnBlade
    !do k=1,DvrData%numBlades
    !   do j=1,AD%u(1)%BladeMotion(k)%nnodes
    !      z = AD%u(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
    !      AD%u(1)%InflowOnBlade(1,j,k) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
    !      AD%u(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
    !      AD%u(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W      
    !   end do !j=nnodes
    !end do !k=numBlades
      
    !InflowOnTower
    !do j=1,AD%u(1)%TowerMotion%nnodes
    !   z = AD%u(1)%TowerMotion%Position(3,j) + AD%u(1)%TowerMotion%TranslationDisp(3,j)
    !   AD%u(1)%InflowOnTower(1,j) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
    !   AD%u(1)%InflowOnTower(2,j) = 0.0_ReKi !V
    !   AD%u(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
    !end do !j=nnodes
                     
end subroutine Set_AD_Inputs
!----------------------------------------------------------------------------------------------------------------------------------
function GetU( URef, ZRef, PLExp, z ) result (U)
   real(ReKi), intent(in) :: URef
   real(ReKi), intent(in) :: ZRef
   real(ReKi), intent(in) :: PLExp
   real(ReKi), intent(in) :: z
   real(ReKi)             :: U
   
   U = URef*(z/ZRef)**PLExp

end function GetU
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_ReadInputFile(fileName, DvrData, errStat, errMsg )
   ! This routine opens the gets the data from the input files.

   character(*),                  intent( in    )   :: fileName
   type(Dvr_SimData),             intent(   out )   :: DvrData
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   

      ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: inpVersion                               ! String containing the input-version information.
   character(1024)              :: line                                     ! String containing a line of input.
   integer                      :: unIn, unEc
   integer                      :: ICase
   integer                      :: Sttus
   character( 11)               :: DateNow                                  ! Date shortly after the start of execution.
   character(  8)               :: TimeNow                                  ! Time of day shortly after the start of execution.
   
   integer, parameter           :: NumCols = 7                              ! number of columns to be read from the input file
   real(DbKi)                   :: InpCase(NumCols)                         ! Temporary array to hold combined-case input parameters. (note that we store in double precision so the time is read correctly)
   logical                      :: TabDel      
   logical                      :: echo   

   INTEGER(IntKi)               :: ErrStat2                                 ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                  ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'Dvr_ReadInputFile'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1
   
   ! Open the input file
   CALL GetPath( fileName, PriPath )     ! Input files will be relative to the path where the primary input file is located.

   call GetNewUnit( unIn )   
   call OpenFInpFile( unIn, fileName, errStat2, ErrMsg2 )
   call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if ( errStat >= AbortErrLev ) then
      call cleanup()
      return
   end if

   
   call WrScr( 'Opening input file:  '//fileName )

      ! Skip a line, read the run title information.

   CALL ReadStr( UnIn, fileName, inpVersion, 'inpVersion', 'File Header: (line 1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL ReadStr( UnIn, fileName, DvrData%OutFileData%runTitle, 'runTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   call WrScr1 ( ' '//DvrData%OutFileData%runTitle )
   
      ! Read in the title line for the input-configuration subsection.
   CALL ReadStr( UnIn, fileName, line, 'line', 'File Header: (line 3)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! See if we should echo the output.     
   call ReadVar ( unIn, fileName, echo, 'Echo', 'Echo Input', errStat2, errMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( echo )  then
         ! Get date and time.
      dateNow = CurDate()
      timeNow = CurTime()
      call GetNewUnit( unEc ) 
      call getroot(fileName,DvrData%OutFileData%Root)      
      call  OpenFOutFile ( unEc, trim( DvrData%OutFileData%Root )//'.ech', errStat2, errMsg2 )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
         if ( errStat >= AbortErrLev ) then
            call cleanup()
            return
         end if
      
      write (unEc,'(A)')      'Echo of Input File:'
      write (unEc,'(A)')      ' "'//fileName//'"'
      write (unEc,'(A)')      'Generated on: '//trim( dateNow )//' at '//trim( timeNow )//'.'
      write (unEc,'(A)')      inpVersion
      write (unEc,'(A)')      DvrData%OutFileData%runTitle
      write (unEc,'(A)')      line
      write (unEc,Ec_LgFrmt)  echo, 'Echo', 'Echo input parameters to "rootname.ech"?'
   end if


      ! Read the rest of input-configuration section.
      
   call ReadVar ( unIn, fileName, DvrData%AD_InputFile,   'AD_InputFile',   'Name of the AeroDyn input file', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%AD_InputFile ) ) DvrData%AD_InputFile = TRIM(PriPath)//TRIM(DvrData%AD_InputFile)

   
      ! Read the turbine-data section.

   call ReadCom ( unIn, fileName, 'the turbine-data subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%numBlades,'NumBlades','Number of blades', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubRad,   'HubRad',   'Hub radius (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubHt,    'HubHt',    'Hub height (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%Overhang, 'Overhang',  'Overhang (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%ShftTilt, 'ShftTilt',  'Shaft tilt (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%ShftTilt = DvrData%ShftTilt*D2R
   call ReadVar ( unIn, fileName, DvrData%precone, 'Precone',  'Precone (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%precone = DvrData%precone*D2R
            
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if           


      ! Read the I/O-configuration section.

   call ReadCom ( unIn, fileName, 'the I/O-configuration subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%OutFileData%Root, 'OutFileRoot', 'Root name for any output files', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if (len_trim(DvrData%OutFileData%Root) == 0) then
      call getroot(fileName,DvrData%OutFileData%Root)
   end if
   
   call ReadVar ( unIn, fileName, TabDel,   'TabDel',   'Make output tab-delimited (fixed-width otherwise)?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if (TabDel) then
         DvrData%OutFileData%delim = TAB
      else
         DvrData%OutFileData%delim = " "
      end if
               
      ! OutFmt - Format used for text tabular output (except time).  Resulting field should be 10 characters. (-):
   call ReadVar( UnIn, fileName, DvrData%OutFileData%OutFmt, "OutFmt", "Format used for text tabular output (except time).  Resulting field should be 10 characters. (-)", ErrStat2, ErrMsg2, UnEc)
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library            
   call ReadVar ( unIn, fileName, Beep,  'Beep',     'Beep on exit?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if


      ! Read the combined-case section.

   call ReadCom  ( unIn, fileName, 'the combined-case subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar  ( unIn, fileName, DvrData%NumCases, 'NumCases', 'Number of cases to run', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (names)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (units)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
      
   if ( DvrData%NumCases < 1 )  then
      call setErrStat( ErrID_Fatal,'Variable "NumCases" must be > 0.' ,errstat,errmsg,routinename)
      call cleanup()
      return
   end if
   
   allocate ( DvrData%Cases(DvrData%NumCases) , STAT=Sttus )
   if ( Sttus /= 0 )  then
      call setErrStat( ErrID_Fatal,'Error allocating memory for the Cases array.',errstat,errmsg,routinename)
      call cleanup()
      return
   end if

   do ICase=1,DvrData%NumCases

      call ReadAry ( unIn, fileName, InpCase,  NumCols, 'InpCase',  'parameters for Case #' &
                     //trim( Int2LStr( ICase ) )//'.', errStat2, errMsg2, UnEc )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
            
      DvrData%Cases(iCase)%WndSpeed        = InpCase( 1)
      DvrData%Cases(ICase)%ShearExp        = InpCase( 2)
      DvrData%Cases(ICase)%RotSpeed        = InpCase( 3)*RPM2RPS
      DvrData%Cases(ICase)%Pitch           = InpCase( 4)*D2R
      DvrData%Cases(ICase)%Yaw             = InpCase( 5)*D2R
      DvrData%Cases(iCase)%dT              = InpCase( 6)
      DvrData%Cases(iCase)%Tmax            = InpCase( 7)
               
   end do ! ICase
   
   call cleanup ( )


   RETURN
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ValidateInputs(DvrData, errStat, errMsg)

   type(Dvr_SimData),             intent(in)    :: DvrData
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None

   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Turbine Data:
   if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%numBlades > 3 ) call SetErrStat( ErrID_Fatal, "There can be no more than 3 blades (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubRad < 0.0_ReKi .or. EqualRealNos(DvrData%HubRad, 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, "HubRad must be a positive number.", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubHt < DvrData%HubRad ) call SetErrStat( ErrID_Fatal, "HubHt must be at least HubRad.", ErrStat, ErrMsg, RoutineName)
   
      
      ! I-O Settings:
      ! Check that DvrData%OutFileData%OutFmt is a valid format specifier and will fit over the column headings
   call ChkRealFmtStr( DvrData%OutFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( FmtWidth /= ChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
      TRIM(Num2LStr(FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.', ErrStat, ErrMsg, RoutineName )

      ! Combined-Case Analysis:
   do i=1,DvrData%NumCases
   
      if (DvrData%Cases(i)%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0 in case '//trim(num2lstr(i))//'.',ErrStat, ErrMsg,RoutineName)
      if (DvrData%Cases(i)%TMax < DvrData%Cases(i)%DT ) call SetErrStat(ErrID_Fatal,'TMax must be larger than dT in case '//trim(num2lstr(i))//'.',ErrStat, ErrMsg,RoutineName)
      
   end do
   
   
   
end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputLine(OutFileData, t, output, errStat, errMsg)

   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_OutputFile)   ,  intent(in   )   :: OutFileData
   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      
   ! Local variables.

   character(200)                   :: frmt                                      ! A string to hold a format specifier
   character(15)                    :: tmpStr                                    ! temporary string to print the time output as text
   integer :: numOuts
   
   errStat = ErrID_None
   errMsg  = ''
   numOuts = size(output,1)
   frmt = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.4)' ) t
   call WrFileNR( OutFileData%unOutFile, tmpStr )
   call WrNumAryFileNR ( OutFileData%unOutFile, output,  frmt, errStat, errMsg )
   if ( errStat >= AbortErrLev ) return
   
     ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')
      
end subroutine Dvr_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_InitializeOutputFile( iCase, CaseData, OutFileData, errStat, errMsg)
      type(Dvr_OutputFile),     intent(inout)   :: OutFileData 
      
      integer(IntKi)         ,  intent(in   )   :: iCase                ! case number (to write in file description line and use for file name)
      type(Dvr_Case),           intent(in   )   :: CaseData
      
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None

         ! locals
      integer(IntKi)                            ::  i      
      integer(IntKi)                            :: numOuts
      
      
      
      call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) then
            OutFileData%unOutFile = -1
            return
         end if
         

      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.'//trim(num2lstr(iCase))//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         
      write (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(version))
      write (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
      write (OutFileData%unOutFile,'()' )    !print a blank line
     ! write (OutFileData%unOutFile,'(A,11(1x,A,"=",ES11.4e2,1x,A))'   ) 'Case '//trim(num2lstr(iCase))//':' &
      write (OutFileData%unOutFile,'(A,11(1x,A,"=",A,1x,A))'   ) 'Case '//trim(num2lstr(iCase))//':' &
         , 'WndSpeed', trim(num2lstr(CaseData%WndSpeed)), 'm/s;' &
         , 'ShearExp', trim(num2lstr(CaseData%ShearExp)), ';' &
         , 'RotSpeed', trim(num2lstr(CaseData%RotSpeed*RPS2RPM)),'rpm;' &
         , 'Pitch',    trim(num2lstr(CaseData%Pitch*R2D)), 'deg;' &
         , 'Yaw',      trim(num2lstr(CaseData%Yaw*R2D)), 'deg;' &
         , 'dT',       trim(num2lstr(CaseData%dT)), 's;' &
         , 'Tmax',     trim(num2lstr(CaseData%Tmax)),'s'
      
      write (OutFileData%unOutFile,'()' )    !print a blank line
              

      numOuts = size(OutFileData%WriteOutputHdr)
         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      call WrFileNR ( OutFileData%unOutFile, '     Time           ' )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputHdr(i) )
      end do ! i

      write (OutFileData%unOutFile,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      call WrFileNR ( OutFileData%unOutFile, '      (s)           ' )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputUnt(i) )
      end do ! i

      write (OutFileData%unOutFile,'()')      
      

      
end subroutine Dvr_InitializeOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
end module AeroDyn_Driver_Subs