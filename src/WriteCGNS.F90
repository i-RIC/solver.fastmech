    MODULE WriteCGNS
    USE RivVarMod
    USE RivVarWMod
    USE RivVarVertMod
    USE CalcCond
    USE READCGNS
    USE iric
    IMPLICIT NONE

    INTEGER :: tmpnx, tmpny

#if 5835
    integer :: gridid_2d = 1
    integer :: gridid_3d
#endif

    CONTAINS
    !	SUBROUTINE Write_CGNS3D_SolGrid()
    !IMPLICIT NONE
    !DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    !INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    !INTEGER :: ier
    !
    !      ALLOCATE(tmpx3d(ns2, nn, nz), STAT = ier)
    !      ALLOCATE(tmpy3d(ns2, nn, nz), STAT = ier)
    !      ALLOCATE(tmpz3d(ns2, nn, nz), STAT = ier)
    !
    !      CALL GETXYZ3DOUT(tmpx3d, tmpy3d, tmpz3d)
    !      call cg_iRIC_Write_Sol_GridCoord3d_f(tmpx3d, tmpy3d, tmpz3d, ier)
    !
    !      DEALLOCATE(tmpx3d, STAT = ier)
    !      DEALLOCATE(tmpy3d, STAT = ier)
    !      DEALLOCATE(tmpz3d, STAT = ier)
    !
    !  END SUBROUTINE Write_CGNS3D_SolGrid

#if 0
    SUBROUTINE Write_CGNS3D_Grid()
    IMPLICIT NONE
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER :: ier
    CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME
    INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
    INTEGER :: NBASES, BASES_ITER
    INTEGER :: BASE3DID, ZONE3DID
    INTEGER :: NSOLs, SOL_ITER, sindex
    INTEGER :: NGRIDS, GRIDS_ITER
    INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER, FSOL_ITER
    INTEGER :: CELLDIM, PHYSDIM
    INTEGER :: GIndex, CIndex
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax
    CHARACTER(LEN = 250) :: name
    INTEGER, DIMENSION(2) :: idata
    
    ALLOCATE(tmpx3d(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpy3d(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpz3d(ns2, nn, nz), STAT = ier)
    
    CALL GETXYZ3DOUT(tmpx3d, tmpy3d, tmpz3d)
    
    CALL cg_nbases_f(FID, NBASES, IER)
    IF(NBASES.EQ.1) THEN
        celldim = 3
        physdim = 3
        isize(1,1) = ns2
        isize(2,1) = nn
        isize(3,1) = nz
    
        isize(1,2) =  isize(1,1)-1
        isize(2,2) = isize(2,1)-1
        isize(3,2) = isize(3,1)-1
    
        isize(1,3) = 0
        isize(2,3) = 0
        isize(3,3) = 0
    
        CALL cg_base_write_f(FID, "iRIC3D", celldim, physdim, BASE3DID, ier )
        CALL cg_zone_write_f(FID, BASE3DID, "iRICZone", isize, Structured, ZONE3DID , ier )
    
        CALL cg_grid_write_f(FID, BASE3DID, ZONE3DID, "GridCoordinates", GIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateX", tmpx3d, CIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateY", tmpy3d, CIndex, ier);
        CALL cg_coord_write_f(FID, BASE3DID, ZONE3DID, RealDouble, "CoordinateZ", tmpz3d, CIndex, ier);
    
        CAll cg_biter_write_f(FID, BASE3DID, 'BaseIterativeData', 1, ier);
    
    ELSE !iRIC3D already exists so delete BaseIterative, ZoneIterative and Solutions
        IF(NBASES.ne.2) THEN
            STOP 'Error in 3D CGNS FILE - Bases'
            !PAUSE
        ENDIF
    
        BASES_ITER = 2 !Hardwired should be iRIC3D
        ZONES_ITER = 1 !Hardwired should be iRICZone
    
        CALL cg_goto_f(FID, 2, ier, 'END')
        CALL cg_delete_node_f('BaseIterativeData', ier)
        CALL cg_nzones_f(FID, BASES_ITER, NZONES, IER)
        IF(NZONES.ne.1) THEN
            STOP 'Error in 3D CGNS File - Zones'
            !PAUSE
        ENDIF
        CALL cg_goto_f(FID, BASES_ITER, ier, 'Zone_t', ZONES_ITER, 'end')
        !Delete existing solutions
        CALL cg_nsols_F(FID, BASES_ITER, ZONES_ITER, nsols, ier);
        DO SOL_ITER = 1,nsols
            CALL cg_sol_info_f(FID, BASES_ITER, ZONES_ITER, 1, solname, sindex, ier)
            CALL cg_delete_node_f(TRIM(solname), ier);
        ENDDO
        !Delete existing solution grids
        CALL cg_ngrids_F(FID, BASES_ITER, ZONES_ITER, ngrids, ier);
        DO GRIDS_ITER = 1,ngrids-1
            CALL cg_grid_read_f(FID, BASES_ITER, ZONES_ITER, 2, solname, ier)
            CALL cg_delete_node_F(TRIM(solname), ier);
        ENDDO
    
        !Delete ZoneIterative Node
        CALL cg_goto_f(FID, BASES_ITER, ier, "Zone_t", ZONES_ITER, "end")
        CALL cg_delete_node_f('ZoneIterativeData', ier)
    
    ENDIF
    
    DEALLOCATE(tmpx3d, tmpy3d, tmpz3d, STAT = ier)
    
    END SUBROUTINE Write_CGNS3D_Grid
#endif

#if 5835
    subroutine write_cgns3d_grid()
    implicit none

    double precision, allocatable, dimension(:,:,:) :: tmpx3d
    double precision, allocatable, dimension(:,:,:) :: tmpy3d
    double precision, allocatable, dimension(:,:,:) :: tmpz3d
    integer :: ier

    allocate(tmpx3d(ns2, nn, nz), stat = ier)
    allocate(tmpy3d(ns2, nn, nz), stat = ier)
    allocate(tmpz3d(ns2, nn, nz), stat = ier)

    call getxyz3dout(tmpx3d, tmpy3d, tmpz3d)
    call cg_iric_write_grid3d_coords_withgridid(fid, ns2, nn, nz, tmpx3d, tmpy3d, tmpz3d, gridid_3d, ier)

    deallocate(tmpx3d, stat = ier)
    deallocate(tmpy3d, stat = ier)
    deallocate(tmpz3d, stat = ier)
    end subroutine write_cgns3d_grid
#endif

#if 0
    SUBROUTINE Write_CGNS3D_FixedBed(solIndex, time, disch)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: solIndex
    REAL(KIND=mp), INTENT(IN) :: time, disch
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tx,  ty, tz, tmpreal4, tmpreal5
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER :: ier, tmpibc, i,j,k
    CHARACTER(LEN = 14) :: flowsol3D
    CHARACTER(LEN = 26) :: gridsol3D
    CHARACTER(LEN = 5) :: tsindex
    CHARACTER(LEN = 32) :: tmp_flowsol3D, tmp_gridsol3D
    CHARACTER(LEN=250) :: ZONENAME, BASENAME, USERNAME, SOLNAME
    INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
    INTEGER :: NBASES, BASES_ITER
    INTEGER :: NSOL, SOL_ITER
    INTEGER :: F2DSOL_ITER, F3DSOL_ITER, F1DSOL_ITER, F2DSOL_EXT_ITER, FSOL_ITER
    INTEGER :: CELLDIM, PHYSDIM
    INTEGER :: GIndex, CIndex
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax
    CHARACTER(LEN = 250) :: name
    INTEGER, DIMENSION(2) :: idata
    INTEGER, DIMENSION(3) :: dimvec


    flowsol3D = '3DFlowSolution'
    gridsol3D = 'GridCoordinatesForSolution'
    WRITE(tsindex, '(I5)') solIndex
    tsindex = ADJUSTL(tsindex)
    tmp_flowsol3D = flowsol3D//TRIM(tsindex)
    tmp_flowsol3D = TRIM(tmp_flowsol3D)
    tmp_gridsol3D = gridsol3D//TRIM(tsindex)
    tmp_gridsol3D = TRIM(tmp_gridsol3D)

    SolNames3D(solIndex) = tmp_flowsol3D
    GridNames3D(solIndex) = tmp_gridsol3D
    TimeIncrements(solIndex) = time
    DischIncrements(solIndex) = disch

    ALLOCATE(tx(ns2, nn, nz), STAT = ier)
    ALLOCATE(ty(ns2, nn, nz), STAT = ier)
    ALLOCATE(tz(ns2, nn, nz), STAT = ier)
    ALLOCATE(tmpint(ns2, nn, nz), STAT = ier)

    !        // write base iterative data, with new solid to the current base.
    !	    ier = cg_biter_write(fileid, baseid, BINAME, solid);
    !	    if (ier != 0){return ier;}
    !	    ier = local_gotobaseiterative();
    !	    if (ier != 0){return ier;}
    !	    // write time.
    !	    dimVec[0] = solid;
    !	    ier = cg_array_write("TimeValues", RealDouble, 1, dimVec, timearray);



    CALL cg_nbases(FID, NBASES, IER)
    IF(NBASES.ne.2) THEN
        STOP 'ERROR in 3D base node'
        !PAUSE
    ENDIF
    ! Building iRIC3D node
    BASES_ITER = 2
    ZONES_ITER = 1

    !write base iter time data
    CALL cg_biter_write(FID, BASES_ITER, 'BaseIterativeData', solIndex, ier)
    CALL cg_goto(FID, BASES_ITER, ier, 'BaseIterativeData_t', 1, "end")
    CALL cg_array_write('TimeValues', RealDouble, 1, solindex, timeincrements, ier)

    CALL GETXYZ3DOUT(tx, ty, tz)
    CALL cg_grid_write(FID, BASES_ITER, ZONES_ITER, tmp_gridsol3D, GIndex, ier)
    CALL cg_goto(FID, BASES_ITER, ier, 'Zone_t', ZONES_ITER, 'GridCoordinates_t', solIndex+1, 'end')
    dimvec(1) = ns2
    dimvec(2) = nn
    dimvec(3) = nz
    CALL cg_array_write('CoordinateX', RealDouble, 3, dimVec, tx, ier)
    CALL cg_array_write('CoordinateY', RealDouble, 3, dimVec, ty, ier)
    CALL cg_array_write('CoordinateZ', RealDouble, 3, dimVec, tz, ier)

    Call GET3DVELOCITYOUT(tx, ty, tz)
    CALL cg_sol_write(FID, BASES_ITER, ZONES_ITER, tmp_flowsol3D, Vertex, F3DSOL_ITER, IER)
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityX', tx, FSOL_ITER, IER)
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityY', ty, FSOL_ITER, IER)
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'VelocityZ', tz, FSOL_ITER, IER)

    CALL GET3DVELOCITYSMAGOUT(tx)
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'SVelocity', tx, FSOL_ITER, IER)

    CALL GET3DVELOCITYNMAGOUT(tx)
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, RealDouble, 'NVelocity', tx, FSOL_ITER, IER)

    DO i = 1,ns2
        DO j = 1,nn
            DO k = 1,nz
                tmpibc = ibc(i,j)
                if(tmpibc > 0) then
                    tmpibc = -1
                endif
                tmpint(i,j,k) = tmpibc
            ENDDO
        ENDDO
    ENDDO
    CALL cg_field_write(FID, BASES_ITER, ZONES_ITER, F3DSOL_ITER, INTEGER, 'IBC', tmpint, FSOL_ITER, IER)

    CALL cg_ziter_write(FID,BASES_ITER,ZONES_ITER,'ZoneIterativeData',ier)
    call cg_goto(FID,BASES_ITER,ier,'Zone_t',ZONES_ITER,'ZoneIterativeData_t',1,'end')
    idata(1)=32
    idata(2)=solIndex
    call cg_array_write('FlowSolutionPointers',Character,2,idata,solnames3D,ier)
    call cg_array_write('GridCoordinatesPointers', Character, 2, idata, gridnames3D, ier)

    DEALLOCATE(tx, ty, tz, tmpint, STAT = ier)
    END SUBROUTINE Write_CGNS3D_FixedBed
#endif

#if 5835
    subroutine write_cgns3d_fixedbed(solindex, time, disch)
    implicit none

    integer, intent(in) :: solindex
    real(kind=mp), intent(in) :: time
    real(kind=mp), intent(in) :: disch

    real(kind=mp), allocatable, dimension(:,:,:) :: tx
    real(kind=mp), allocatable, dimension(:,:,:) :: ty
    real(kind=mp), allocatable, dimension(:,:,:) :: tz

    integer, allocatable, dimension(:,:,:) :: tmpint
    integer :: ier
    integer :: tmpibc
    integer :: i
    integer :: j
    integer :: k

    allocate(tx(ns2, nn, nz), stat = ier)
    allocate(ty(ns2, nn, nz), stat = ier)
    allocate(tz(ns2, nn, nz), stat = ier)
    allocate(tmpint(ns2, nn, nz), stat = ier)

    call getxyz3dout(tx, ty, tz)
    call cg_iric_write_sol_grid3d_coords(fid, tx, ty, tz, ier)

    call get3dvelocityout(tx, ty, tz)
    call cg_iric_write_sol_node_real_withgridid(fid, gridid_3d, 'VelocityX', tx, ier)
    call cg_iric_write_sol_node_real_withgridid(fid, gridid_3d, 'VelocityY', ty, ier)
    call cg_iric_write_sol_node_real_withgridid(fid, gridid_3d, 'VelocityZ', tz, ier)

    call get3dvelocitysmagout(tx)
    call cg_iric_write_sol_node_real_withgridid(fid, gridid_3d, 'SVelocity', tx, ier)

    call get3dvelocitynmagout(tx)
    call cg_iric_write_sol_node_real_withgridid(fid, gridid_3d, 'NVelocity', tx, ier)

    do i = 1, ns2
        do j = 1, nn
            do k = 1, nz
                tmpibc = ibc(i, j)
                if(tmpibc > 0) then
                    tmpibc = -1
                endif
                tmpint(i, j, k) = tmpibc
            enddo
        enddo
    enddo
    call cg_iric_write_sol_node_integer_withgridid(fid, gridid_3d, 'IBC', tmpint, ier)

    deallocate(tx, ty, tz, tmpint, stat = ier)
    end subroutine write_cgns3d_fixedbed
#endif


#if 0
    SUBROUTINE Write_CGNS3D_MoveableBed(time, disch)
    IMPLICIT NONE
    REAL(KIND=mp), INTENT(IN) :: time, disch
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tx,  ty, tz, tmpreal4, tmpreal5
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER :: ier

    !		CALL CG_IRIC_WRITE_SOL_TIME_F(time, ier)
    ALLOCATE(tx(ns2, nn, nz), STAT = ier)
    ALLOCATE(ty(ns2, nn, nz), STAT = ier)
    ALLOCATE(tz(ns2, nn, nz), STAT = ier)

    call cg_iric_write_sol_time(fid, time, ier);

    CALL GETXYZ3DOUT(tx, ty, tz)
    call cg_iRIC_Write_Grid3d_Coords(tx, ty, tz, ier)

    Call GET3DVELOCITYOUT(tx, ty, tz)
    call cg_iric_write_grid_real_node(fid, 'VelocityX', tx, ier)
    call cg_iric_write_grid_real_node(fid, 'VelocityY', ty, ier)
    call cg_iric_write_grid_real_node(fid, 'VelocityZ', tz, ier)

    CALL GET3DVELOCITYSMAGOUT(tx)
    call cg_iric_write_grid_real_node(fid, 'SVelocity', tx, ier)

    CALL GET3DVELOCITYNMAGOUT(tx)
    call cg_iric_write_grid_real_node(fid, 'NVelocity', tx, ier)

    DEALLOCATE(tx, ty, tz, STAT = ier)

    END SUBROUTINE Write_CGNS3D_MoveableBed
#endif

#if 5835
    ! write_cgns3d_moveablebed is never called
#endif

    SUBROUTINE Write_CGNS2(time, disch)
    IMPLICIT NONE
    REAL(KIND=mp), INTENT(IN) :: time, disch
    INTEGER :: IER, iret
    INTEGER :: I,J
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: tmpreal1,  tmpreal2
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpreal1a,  tmpreal2a, tmpreal3, tmpreal4, tmpreal5
    REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpx3d, tmpy3d, tmpz3d
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmpint
    INTEGER, ALLOCATABLE, DIMEnSION(:,:) :: tmp2dint
    !		CALL CG_IRIC_WRITE_SOL_TIME_F(time, ier)

    ALLOCATE(tmpreal1(ns2, nn), STAT = ier)
    ALLOCATE(tmpreal2(ns2, nn), STAT = ier)
    ALLOCATE(tmp2dint(ns2, nn), STAT = ier)

    IF(CALCQUASI3D.and.IO_3DOUTPUT) THEN
        CALL getXYOut(tmpreal1, tmpreal2)
        CALL cg_iRIC_Write_Sol_Grid2d_Coords(fid, tmpreal1,tmpreal2,IER)
    ENDIF

    CALL GETIBCOUT(tmp2dint)
    CALL cg_iric_write_sol_node_integer_withgridid(fid, gridid_2d, 'IBC', tmp2dint, IER)

    CALL GETFMIBCOUT(tmp2dint)
    CALL cg_iric_write_sol_node_integer_withgridid(fid, gridid_2d, 'FMIBC', tmp2dint, IER)

    CALL getVelocityOut(tmpreal1, tmpreal2)
    CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityX', tmpreal1, IER)
    CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityY', tmpreal2, IER)

    CALL getDepthOut(tmpreal1)
    CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Depth', tmpreal1, IER)

    CALL getWSEOut(tmpreal1)
    CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'WaterSurfaceElevation', tmpreal1, IER)

    CALL getElevationOut(tmpreal1)
    CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Elevation', tmpreal1, IER)

    IF(IO_VelSN) THEN
        CALL GETVELOCITYSNOUT(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityS', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityN', tmpreal2, IER)
    ENDIF

    IF(IO_UnitDisch) THEN
        CALL getUnitDischOut(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'UnitDischargeX', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'UnitDischargeY', tmpreal2, IER)
    ENDIF

    IF(IO_HArea) THEN
        CALL getHAreaOut(tmpreal1)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'HabitatArea', tmpreal1, IER)
    ENDIF

    IF(IO_InitVel) THEN
        CALL getInitVelOut(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityInitS', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'VelocityInitN', tmpreal2, IER)
    ENDIF

    IF(IO_ShearXY) THEN
        CALL GETSHEARSTRESSOUT(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'ShearStressX', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'ShearStressY', tmpreal2, IER)
    ENDIF

    IF(IO_ShearSN) THEN
        CALL GETSHEARSTRESSSNOUT(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'ShearStressS', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'ShearStressN', tmpreal2, IER)
    ENDIF

    IF(IO_CD) THEN
        CALL getCDOut(tmpreal1)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Drag_Coefficient', tmpreal1, IER)
    ENDIF

    IF(CALCCSED) THEN
        IF(TRANSEQTYPE == 2) THEN
            CALL GETSANDDEPTHOUT(tmpreal1)
            CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Sand_Depth', tmpreal1, IER)
            CALL GETSANDFRACOUT(tmpreal1)
            CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Sand_Fraction', tmpreal1, IER)
            CALL GETLSUBHOUT(tmpreal1)
            CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'LSub', tmpreal1, IER)

        ENDIF
        CALL GetTransportRateOut(tmpreal1, tmpreal2)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'SedFluxX', tmpreal1, IER)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'SedFluxY', tmpreal2, IER)
    ENDIF

    IF(IO_StressDiv) THEN
        CALL GetStressDivOut(tmpreal1)
        IF(CALCCSED) THEN
            CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Erosion Rate', tmpreal1, IER)
        ELSE
            CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Shear Stress Divergence', tmpreal1, IER)
        ENDIF
    ENDIF

    IF(CALCQUASI3D.and.IO_HELIX) THEN
        CALL getHelixOut(tmpreal1)
        CALL cg_iric_write_sol_node_real_withgridid(fid, gridid_2d, 'Helix Strength', tmpreal1, IER)
        !            CALL getRSOut(tmpreal1)
        !            CALL cg_iRIC_Write_Sol_Real_f("RS", tmpreal1, IER)
        !            CALL getThetaOut(tmpreal1)
        !            CALL cg_iRIC_Write_Sol_Real_f("Theta", tmpreal1, IER)
    ENDIF


    DEALLOCATE(tmpreal1, STAT=ier)
    DEALLOCATE(tmpreal2, STAT=ier)
    DEALLOCATE(tmp2dint, STAT=ier)
    !        IF(NZ.GT.0) THEN
    !!            CALL CG_IRIC_CREATE_3DBASEZONE_F(FID, NS, NN, NZ, X, Y, ZZ, IER)
    !!            CALL CG_IRIC_CREATE_3DFLOWSOL_F(FID, INDEX, NS, NN, NZ, IER)
    !!            CALL cg_iRIC_Write_3DBaseIter(FID, INDEX, TimeIncrements, DischIncrements)
    !!            CALL cg_iRIC_Write_3DZoneIter(FID, index, SolNames3D)
    !            ALLOCATE(tmpx3d(ns, nn, nz), STAT = ier)
    !            ALLOCATE(tmpy3d(ns, nn, nz), STAT = ier)
    !            ALLOCATE(tmpz3d(ns, nn, nz), STAT = ier)
    !            CALL GETXY3DOUT(tmpx3d, tmpy3d)
    !            CALL GET3DGRIDZOUT(tmpz3d)
    !
    !            CALL cg_iric_write_sol_time_f(time, ier)
    !            Call cg_iric_write_sol_gridcoord3d_f(tmpx3d, tmpy3d, tmpz3d, ier)
    !!
    !!            ALLOCATE(tmpreal2a(ns, nn, nz), STAT = ier)
    !!            ALLOCATE(tmpreal3(ns, nn, nz), STAT = ier)
    !!            ALLOCATE(tmpreal4(ns, nn, nz), STAT = ier)
    !!            ALLOCATE(tmpreal5(ns, nn, nz), STAT = ier)
    !!            ALLOCATE(tmpint(ns, nn, nz), STAT = ier)
    !!
    !!            CALL  GET3DVELOCITYOUT(tmpreal1a, tmpreal2a, tmpreal3, tmpreal4, tmpreal5)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "VelocityX", tmpreal1a, IER)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "VelocityY", tmpreal2a, IER)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "VelocityZ", tmpreal3, IER)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "VelocityS", tmpreal4, IER)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "VelocityN", tmpreal5, IER)
    !!
    !!            CALL GET3DIBCOUT(tmpint)
    !!            CALL cg_iRIC_Write_3DSolIntegerNode_f(Index, NS, NN, NZ, "IBC", tmpint, IER)
    !!
    !!            CALL GET3DGRIDZOUT(tmpreal1a)
    !!            CALL cg_iRIC_Write_3DSolRealNode_f(Index, NS, NN, NZ, "GridZ", tmpreal1a, IER)
    !!
    !!            DEALLOCATE(tmpreal1a, STAT=ier)
    !!            DEALLOCATE(tmpreal2a, STAT=ier)
    !!            DEALLOCATE(tmpreal3, STAT=ier)
    !!            DEALLOCATE(tmpreal4, STAT=ier)
    !!            DEALLOCATE(tmpreal5, STAT=ier)
    !!            DEALLOCATE(tmpint, STAT=ier)
    !!
    !        ENDIF
    !    call cg_close_f(FID, IER)
    END SUBROUTINE




    SUBROUTINE GET3DGRIDZOUT(gridz)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: gridz
    INTEGER :: I,J,K
    DO k = 1,nz
        DO j= 1,nn
            DO i=1,ns2
                gridz(i,j,k) = (zz(i,j,k)/100.) + elevoffset
            END DO
        END DO
    END DO


    END SUBROUTINE

    SUBROUTINE GETIBCOUT(iibc)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO j= 1,nn
        DO i=1,ns2
            if(ibc(i,j).gt.0) then
                iibc(i,j) = -1
            else
                iibc(i,j) = ibc(i,j)
            endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETFMIBCOUT(iibc)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO j= 1,nn
        DO i=1,ns2
            !		        if(ibc(i,j).gt.0) then
            !		            iibc(i,j) = -1
            !		        else
            iibc(i,j) = ibc(i,j)
            !		        endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DIBCOUT(iibc)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: iibc
    INTEGER :: I,J,K
    DO k = 1,nz
        DO j= 1,nn
            DO i=1,ns2
                iibc(i,j,k) = ibc(i,j)
            END DO
        END DO
    END DO

    END SUBROUTINE


    SUBROUTINE GET3DVELOCITYOUT(ivelx, ively, ivelz)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: I,J,K
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy

    DO i=1,ns2
        rcos = cos(phirotation(i))
        rsin = sin(phirotation(i))
        DO j= 1,nn
            DO k = 1,nz

                ux = (uz(i,j,k)*rcos - vz(i,j,k)*rsin)/100.
                uy = (uz(i,j,k)*rsin + vz(i,j,k)*rcos)/100.
                xx = ux*fcos - uy*fsin
                yy = ux*fsin + uy*fcos
                ivelx(i,j,k) = xx
                ively(i,j,k) = yy
                ivelz(i,j,k) = 0.0

            END DO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYSOUT(ivelx, ively, ivelz)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: i, j, k

    DO j= 1,nn
        DO i=1,ns2
            DO k = 1,nz
                ivelx(i,j,k) = uz(i,j,k)/100.
                ively(i,j,k) = 0.0
                ivelz(i,j,k) = 0.0
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYNOUT(ivelx, ively, ivelz)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx, ively, ivelz
    INTEGER :: i, j, k

    DO j= 1,nn
        DO i=1,ns2
            DO k=1,nz
                ivelx(i,j,k) = 0.0
                ively(i,j,k) = vz(i,j,k)/100.
                ivelz(i,j,k) = 0.0
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYNMAGOUT(ively)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ively
    INTEGER :: i, j, k

    DO j= 1,nn
        DO i=1,ns2
            DO k=1,nz
                ively(i,j,k) = vz(i,j,k)/100.
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GET3DVELOCITYSMAGOUT(ivelx)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:,:), INTENT(OUT) :: ivelx
    INTEGER :: i, j, k

    DO j= 1,nn
        DO i=1,ns2
            DO k=1,nz
                ivelx(i,j,k) = uz(i,j,k)/100.
            ENDDO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETVELOCITYSNOUT(ivelx, ively)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            ivelx(i,j) = u(i,j)/100.
            ively(i,j) = vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSHEARSTRESSSNOUT(issx, issy)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: issx, issy
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            issx(i,j) = taus(i,j)/10.
            issy(i,j) = taun(i,j)/10.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSHEARSTRESSOUT(SSXOut, SSYOut)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: SSXOut, SSYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,nn
        DO i=1,ns2
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*ns

            ux = (taus(i,j)*rcos - taun(i,j)*rsin)/10.
            uy = (taus(i,j)*rsin + taun(i,j)*rcos)/10.
            xx = ux*fcos - uy*fsin
            yy = ux*fsin + uy*fcos
            SSXOut(i,j) = xx
            SSYOut(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETVELOCITYOUT(velXOut, velYOut)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: velXOut, velYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,nn
        DO i=1,ns2
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*ns

            !				    ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
            !				    uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
            !if(uvinterp == 0) then
            ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
            uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
            !else
            !  ux = (uout(i,j)*rcos - vout(i,j)*rsin)/100.
            !  uy = (uout(i,j)*rsin + vout(i,j)*rcos)/100.
            !endif
            xx = ux*fcos - uy*fsin
            yy = ux*fsin + uy*fcos
            velXOut(i,j) = xx
            velYOut(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETTransportRateOUT(TRXOut, TRYOut)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: TRXOut, TRYOut
    INTEGER :: i,j,count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy


    DO j= 1,nn
        DO i=1,ns2
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*ns

            ux = (qs(i,j)*rcos - qn(i,j)*rsin)
            uy = (qs(i,j)*rsin + qn(i,j)*rcos)
            xx = ux*fcos - uy*fsin
            yy = ux*fsin + uy*fcos
            TRXOut(i,j) = xx/10000.
            TRYOut(i,j) = yy/10000.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETINITVELOUT(ivelx, ively)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            ivelx(i,j) = u(i,j)/100.
            ively(i,j) = v(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETUNITDISCHOUT(ivelx, ively)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: ivelx, ively
    INTEGER :: i, j, count
    REAL(KIND=mp) :: rcos, rsin, xx, yy, ux, uy
    DO j= 1,nn
        DO i=1,ns2
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            !					count = ((i-1)*nn)+j
            count = i + (j-1)*ns

            ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
            uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
            xx = ux*fcos - uy*fsin
            yy = ux*fsin + uy*fcos
            ivelx(i,j) = xx
            ively(i,j) = yy
            !				tmpvar3(count) = u(i,j)/100.
            !				tmpvar4(count) = vout(i,j)/100.
        END DO
    END DO
    DO j= 1,nn
        DO i=1,ns2
            ivelx(i,j) = ivelx(i,j) * hl(i,j)/100.
            ively(i,j) = ively(i,j) * hl(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSANDDEPTHOUT(sd)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: sd
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            sd(i,j) = hfine(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETLSUBHOUT(tlsub)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tlsub
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            tlsub(i,j) = (vsandsub(i,j)/(0.3*harea(i,j)))/100.
        END DO
    END DO

    END SUBROUTINE


    SUBROUTINE getRSOut(irs)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: irs
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            irs(i,j) = rs(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE getThetaOut(itheta)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: itheta
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            itheta(i,j) = theta(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETSANDFRACOUT(sf)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: sf
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            sf(i,j) = Fracs(i,j)
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETCDOUT(tcd)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tcd
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            tcd(i,j) = totcd(i,j)
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETDEPTHOUT(depth)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: depth
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            depth(i,j) = hl(i,j)/100.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETHAREAOUT(htArea)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: htArea
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            htArea(i,j) = harea(i,j)/10000.
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETWSEOUT(wse)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: wse
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            wse(i,j) = (e(i,j)/100.) + elevoffset
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETELEVATIONOUT(elev)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: elev
    INTEGER :: i, j

    DO j= 1,nn
        DO i=1,ns2
            elev(i,j) = (eta(i,j)/100.) + elevoffset
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETXYOUT(tx, ty)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: tx, ty
    REAL(KIND=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    INTEGER :: i, j
    !        DO j= 1,nn
    !		    DO i=1,ns
    ! 			    rcos = cos(phirotation(i))
    !			    rsin = sin(phirotation(i))
    !    !					count = ((i-1)*nn)+j
    !			    count = i + (j-1)*ns
    !
    !				    ux = (u(i,j)*rcos - vout(i,j)*rsin)/100.
    !				    uy = (u(i,j)*rsin + vout(i,j)*rcos)/100.
    !				    xx = ux*fcos - uy*fsin
    !				    yy = ux*fsin + uy*fcos
    !				    velXOut(i,j) = xx
    !				    velYOut(i,j) = yy
    !    !				tmpvar3(count) = u(i,j)/100.
    !    !				tmpvar4(count) = vout(i,j)/100.
    !		    END DO
    !	    END DO

    DO j= 1,nn
        DO i=1,ns2
            rcos = cos(phirotation(i))
            rsin = sin(phirotation(i))
            ux = x(i,j)*rcos - y(i,j)*rsin
            uy = x(i,j)*rsin + y(i,j)*rcos

            xx = x(i,j)*fcos - y(i,j)*fsin
            yy = x(i,j)*fsin + y(i,j)*fcos
            tx(i,j) = (xx + xshift)/100.
            ty(i,j) = (yy + yshift)/100.

        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GETXYZ3DOUT(tx, ty, tz)
    IMPLICIT NONE
    Double Precision, DIMENSION(:,:,:), INTENT(OUT) :: tx, ty, tz
    INTEGER :: i, j, k
    Double Precision :: xx, yy, ux, uy
    double precision :: rcos, rsin

    DO k = 1,nz
        DO j=1,nn
            DO i= 1,ns2
                rcos = cos(phirotation(i))
                rsin = sin(phirotation(i))

                ux = x(i,j)*rcos - y(i,j)*rsin
                uy = x(i,j)*rsin + y(i,j)*rcos

                xx = x(i,j)*fcos - y(i,j)*fsin
                yy = x(i,j)*fsin + y(i,j)*fcos

                tx(i,j,k) = (xx + xshift)/100.
                ty(i,j,k) = (yy + yshift)/100.
                tz(i,j,k) = (zz(i,j,k)/100.) + elevoffset
            END DO
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE getHelixOut(helix)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: helix
    REAL(KIND=mp) :: xx, yy
    INTEGER :: i, j
    DO j= 1,nn
        DO i=1,ns2
            xx = atan2(taun(i,j), taus(i,j))*(360./(2.*3.14159))
            yy = atan2(vz(i,j,nz), uz(i,j, nz))*(360./(2.*3.14159))
            !			    if(vz(i,j,2).eq.0.0.or.uz(i,j,2).eq.0) then
            !			        helix(i,j) = 0
            !			    else
            helix(i,j) = yy-xx
            if(helix(i,j) > 180) then
                helix(i,j) = helix(i,j)-360.
            elseif(helix(i,j) < -180) then
                helix(i,j) = helix(i,j)+360.
            endif
            !		        Endif
        END DO
    END DO

    END SUBROUTINE

    SUBROUTINE GetStressDivOut(stressdiv)
    IMPLICIT NONE
    REAL(KIND=mp), DIMENSION(:,:), INTENT(OUT) :: stressdiv
    INTEGER :: i, j
    DO j= 1,nn
        DO i=1,ns2
            stressdiv(i,j) = con(i,j)
        END DO
    END DO

    END SUBROUTINE
#if 0
    SUBROUTINE write_error_CGNS(OutputFile)
    !		USE USER32
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: OutputFile
    CHARACTER(250) :: ZONENAME, BASENAME, USERNAME

    INTEGER :: FID, BID, ZID, IER, iret
    INTEGER :: NUSER_DATA, USER_ITER, NZONES, ZONES_ITER
    INTEGER :: NBASES, BASES_ITER
    INTEGER :: CELLDIM, PHYSDIM
    INTEGER, DIMENSION(3,3) :: isize
    INTEGER, DIMENSION(3) :: irmin, irmax
    INTEGER, DIMENSION(6) :: irinddata
    CHARACTER(LEN = 250) :: name, errorstring
    INTEGER ::namelen
    INTEGER, DIMENSION(1) :: dimvals
    !		namelen = len(InputFile)
    !		name = TRIM(InputFile(1:namelen-3)//'cgns')
    FID = CGNSFILEID
    !		CALL cg_open_f(OutputFile, MODE_MODIFY, FID, IER)
    !			IF(IER .NE. 0) THEN
    !!				iret = MESSAGEBOX(0, "CGNS ERROR"C, "Error"C, MB_OK)
    !			ENDIF
    CALL cg_open(OutputFile, MODE_MODIFY, FID, IER)
    IF(IER .NE. 0) THEN
        call cg_error_print()
    ENDIF
    BID = 1
    CALL cg_nbases(FID, NBASES, IER)
    DO BASES_ITER = 1, NBASES
        CALL cg_base_read(FID, BASES_ITER, BASENAME, CELLDIM, PHYSDIM, IER)

        SELECT CASE(TRIM(BASENAME))
        CASE('iRIC')
            CALL cg_nzones(FID, BASES_ITER, NZONES, IER)
            DO ZONES_ITER = 1, NZONES
                CALL cg_zone_read(FID, BASES_ITER, ZONES_ITER, ZONENAME, isize, IER)

                SELECT CASE(TRIM(ZONENAME))
                CASE('iRICZone')
                    CALL cg_goto(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER, 'end')
                    CALL cg_nuser_data(NUSER_DATA, IER)
                    DO USER_ITER = 1, NUSER_DATA
                        CALL cg_user_data_read(USER_ITER, USERNAME, IER)

                        SELECT CASE(TRIM(USERNAME))
                        CASE('Error')
                            call cg_goto(FID, BASES_ITER, IER, 'Zone_t', ZONES_ITER,'UserDefinedData_t', USER_ITER, 'end')
                            dimvals(1) = 1
                            CALL cg_array_write('ErrorCode', Integer, 1, 1, errorcode, ier)
                            SELECT CASE (errorcode)
                            CASE (-1)
                                Write(errorstring, '(A70,I5)') 'Error: Parameters chosen do not result in convergence at iteration:', iter
                                dimvals(1) = len(errorstring)
                                CALL cg_array_write('ErrorString', Character, 1, dimvals, errorstring, ier)
                            CASE (-10)
                                !								status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, "r2d_ErrorText",85,
                                !						 +"At least one row of nodes are completely dry - check water surface
                                !						 + boundary condition")
                                Write(errorstring, '(A100)') 'Error: at least on row of nodes is dry'
                                dimvals(1) = len(errorstring)
                                CALL cg_array_write('ErrorString', Character, 1, dimvals, errorstring, ier)

                                CASE DEFAULT
                                !								status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, "r2d_ErrorText",
                                !						 +		23, "An unkown error occured")
                                Write(errorstring, '(A70)') 'An unkown error occured'
                                dimvals(1) = len(errorstring)
                                CALL cg_array_write('ErrorString', Character, 1, dimvals, errorstring, ier)
                            END SELECT
                            !                                CALL cg_close_f(FID, IER)
                            return
                        END SELECT

                    ENDDO
                END SELECT
            ENDDO
        END SELECT
    ENDDO

    CALL cg_iric_close(FID, IER)

    END SUBROUTINE
#endif

    END MODULE
