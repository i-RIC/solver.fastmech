    module fastmech
    use einitmod2
    use rivvarmod2
    use calccond2
    use rivvarwmod2
    use rivvartimemod2
    use gridcond2
    use csedmod2
    use rivvarvertmod2
    use writecgns2
    use rivvertmod2
    use stressdivmod2
    use uinitmod2
    use tridiag2
    use rivcalcinitcond2

    use rivconnectivitymod2
    use fm_helpers2
    implicit none

    type :: fastmech_model
        integer, pointer :: itermax
        real(kind=mp) :: dt
        real(kind=mp) :: t
        real(kind=mp) :: t_end
        integer, pointer :: n_x, n_y, n_z
        !real(kind=mp), pointer, dimension(:,:) :: fm_elevation, fm_depth, fm_wse
        !real(kind=mp), pointer, dimension(:,:) :: fm_velocity_x, fm_velocity_y, rvo%taus, fm_shearstress_y
        !real(kind=mp), pointer, dimension(:,:) :: fm_dragcoefficient, fm_roughness, fm_vegroughness
        type(rivvar) :: t_rivvar
        type(calccond) :: t_calccond
        type(riv_w_var) :: t_rivwvar
        type(rivvartime) :: t_rivvartime
        type(csed) :: t_csed
        type(rivvarvert) :: t_rivvarvert
    end type fastmech_model


    private :: initialize, solve_fm

    contains

    ! initializes the model with values read from a file.
    subroutine initialize_from_file(model, config_file)
    implicit none
    include "cgnslib_f.h"
    include "iriclib_f.h"
    type(rivvar), pointer :: rvo
    type(calccond), pointer :: cco
    type(riv_w_var), pointer :: rvwo
    type(rivvartime), pointer :: rvto
    character (len=*), intent (in) :: config_file
    type(fastmech_model), target, intent(out) :: model
    integer :: tmp, fid, i,j
    real(kind=mp) :: newdisch, newstage
    character(250) :: zonename, basename, username

    integer :: ier, iret, tmpns, tmpnn
    real(kind = mp), allocatable, dimension(:,:) :: tmpx1, tmpy1
    rvo => model%t_rivvar
    cco => model%t_calccond
    rvwo => model%t_rivwvar
    rvto => model%t_rivvartime
    !model%t_csed => m_csed
    !model%t_rivvarvert => m_rvvo

    rvo%g = 980.0d0
    rvo%rho = 1.0d0
    rvo%vkc = 0.40d0
    rvo%noldwetnodes = 0
    rvo%nwetnodes = 0
    rvo%errorcode = 0
    call welcome

    call cg_open_f(config_file, cg_mode_modify, fid, ier)
    if(ier .ne. 0) then
        call cg_error_print_f()
        !pause
        return
    endif
    rvo%cgnsfileid = fid
    ! uncomment lines below for split solution
    !call iric_initoption_f(iric_option_dividesolutions, ier)
    !    if (ier /=0) stop "*** initialize option error***"

    call cg_iric_init_f(fid, ier)
    if(ier .ne. 0) then
        call cg_error_print_f()
        !pause
        return
    endif
    call iric_initoption_f(iric_option_cancel, ier)
    if (ier /=0) stop "*** initialize option error***"

    call cgns2_read_cc_foralloc(cco, rvo, ier)

    call cg_iric_gotogridcoord2d_f(tmpns, tmpnn, ier)
    rvo%ns2 = tmpns
    rvo%ns = tmpns+rvo%nsext
    rvo%nn = tmpnn

    call alloc_common2d(rvo, tmpns, rvo%nsext, tmpnn)
    call alloc_init2d(rvo, rvo%ns2, rvo%nn)
    call alloc_working(rvwo, rvo%ns, rvo%nn, cco%itm)
    !block below represents old gridcordmod
    allocate(tmpx1(tmpns,tmpnn), stat=ier)
    allocate(tmpy1(tmpns,tmpnn), stat=ier)
    call cg_iric_getgridcoord2d_f(tmpx1,tmpy1,ier)
    do i=1,rvo%ns2
        do j=1,rvo%nn
            rvo%x(i,j) = tmpx1(i,j)*100.
            rvo%y(i,j) = tmpy1(i,j)*100.
        enddo
    enddo
    call calcgridmetrics(rvo)
    deallocate(tmpx1, tmpy1, stat=ier)

    call cgns2_read_gridcondition(rvo, ier)
    call cgns2_read_calccondition(cco, rvo, rvto, ier)
    if(cco%transeqtype == 2) then
        !       call alloc_csed_dt()
    else
        call alloc_csed(model%t_csed, rvo%ns, rvo%nn)
    endif
    if(cco%calcquasi3d) then
        call alloc_common3d(rvo)
        if(cco%transeqtype == 2) then
            !call alloc_csed3d_dt() !array for daniele tonina wilcock-kenworthy
        endif
    endif

    rvo%mo=rvo%stot*(rvo%ns-1)/(rvo%ns2-1)
    rvo%wmax=rvo%w2(1)
    rvo%dn=rvo%wmax/(rvo%nn-1)
    rvo%ds=rvo%mo/(rvo%ns-1)
    rvo%nm=(rvo%nn+1)/2.

    call initratingcurves(rvto)
    call inittimeseries(rvto)
    call initarrays(rvo, rvwo, cco)
    call calcwsinitcond(rvo, rvwo, cco, rvto)

    call initialize(model)

    end subroutine initialize_from_file

    ! initializes the model with default hardcoded values.
    subroutine initialize_from_defaults(model)
    implicit none
    type (fastmech_model), intent (out) :: model
    integer :: tmp
    model%t_rivvar%g = 980.0d0
    model%t_rivvar%rho = 1.0d0
    model%t_rivvar%vkc = 0.40d0
    tmp = 1
    !model%alpha = 0.75
    !model%t_end = 20.
    !model%n_x = 10
    !model%n_y = 20
    call initialize(model)
    end subroutine initialize_from_defaults

    !allocates memory and sets values for either initialization technique.
    subroutine initialize(model)
    implicit none
    type(rivvar), pointer :: rvo
    type(rivvartime), pointer :: rvto
    type(calccond), pointer :: cco
    type(rivvarvert), pointer :: rvvo
    type(riv_w_var), pointer :: rvwo
    type (fastmech_model), target, intent (inout) :: model
    integer :: tmp
    real(kind=mp) :: newdisch, newstage
    rvo => model%t_rivvar
    rvto => model%t_rivvartime
    cco => model%t_calccond
    rvvo => model%t_rivvarvert
    rvwo => model%t_rivwvar
    tmp = 1
    call calc_area(rvo%ns, rvo%nn,rvo%phirotation, &
        rvo%x, rvo%y, rvo%xo, rvo%yo, &
        rvo%nm, rvo%dn, rvo%harea)
    if(rvto%vardischtype == 1) then !discharge time series
        call getinterptimeseriesvalue(rvto, 1, rvto%vardischstarttime, newdisch)
        cco%q = newdisch*1e6
    else
        newdisch = cco%q/1e6
    endif
    if(rvto%varstagetype == 1) then !stage time series
        !in this case the times series is stored in the rating curve
        call getinterpratingcurvevalue(rvto, 1, rvto%vardischstarttime, newstage)
    else if (rvto%varstagetype == 2) then !stage rating curve
        call getinterpratingcurvevalue(rvto, 1, newdisch, newstage)
    endif
    if(.not.cco%hotstart) then
        call einit(rvo, cco, rvo%e, rvo%hl, rvo%eta, rvo%ibc, rvo%w, rvo%hav)
        call uinit(rvo, cco, rvo%u,rvo%v,rvo%hav,rvo%w,rvo%eta,cco%q,rvo%ibc)
        if(rvo%errorcode < 0) then
            !call write_error_cgns(str_in)
            call dealloc_common2d(rvo)
            if(cco%vbc == 1) then
                call dealloc_velbc(cco)
            endif
            if(cco%calcquasi3d) then
                call dealloc_common3d(rvo)
                !    if(cco%transeqtype == 2) then
                !        call dealloc_csed3d_dt()
                !    endif
            endif
            call dealloc_working(rvwo)
            !pause
            return
        endif
    endif
    if(cco%roughnesstype == 1) then !roughness as z0
        !		call zotocdtwo() !calculate cd from log profile
        call z0tocd(rvo, rvvo, cco) !calculate cd from 2-part profile
    endif
    if(rvto%vardischendtime == 0 .and. rvto%vardischstarttime == 0)then
        rvo%nsteps = 0
    else
        rvo%nsteps = (rvto%vardischendtime-rvto%vardischstarttime)/cco%fmdt
    endif

    if(cco%iterationout) then
        if(cco%itm > rvo%nsteps) then
            call alloc_tsnames(rvto, cco%itm)
        else
            call alloc_tsnames(rvto, rvo%nsteps+10)
        endif
    else
        call alloc_tsnames(rvto, rvo%nsteps+10)
    endif
    rvo%oldstage = rvo%hav(rvo%ns)
    rvo%dsstage = rvo%hav(rvo%ns)
    !    endif
    rvo%newstage = rvo%oldstage
    rvo%olddisch = cco%q

    rvo%nct = -1
    model%t = rvto%vardischstarttime
    model%dt = cco%fmdt
    model%t_end = rvto%vardischendtime
    rvo%ptime = model%t+cco%iplinc*model%dt

    if(cco%soltype == 0) then
        ! set these to run a single time step
        model%t = 0
        model%t_end = 1
        model%dt = 1
    endif
    model%itermax => cco%itm
    model%n_x => rvo%ns2
    model%n_y => rvo%nn
    model%n_z => rvo%nz
    !model%fm_elevation => rvo%eta
    !model%fm_depth => rvo%hl
    !model%fm_wse => rvo%e
    !model%fm_velocity_x => rvo%u
    !model%fm_velocity_y => rvo%v
    !model%fm_shearstress_x => rvo%taus
    !model%fm_shearstress_y => rvo%taun
    !model%rvo%totcd => rvo%totcd
    !model%fm_roughness => rvo%cd
    !model%fm_vegroughness => rvo%cdv
    !model%t = 0.
    !model%dt = 1.
    !model%dx = 1.
    !model%dy = 1.
    !
    !allocate(model%temperature(model%n_y, model%n_x))
    !allocate(model%temperature_tmp(model%n_y, model%n_x))
    !
    !model%temperature = 0.
    !model%temperature_tmp = 0.
    !
    !call set_boundary_conditions(model%temperature)
    !call set_boundary_conditions(model%temperature_tmp)
    end subroutine initialize

    ! frees memory when program completes.
    subroutine cleanup(model)
    type (fastmech_model), target, intent (inout) :: model
    integer :: ier
    type(rivvar), pointer :: rvo
    type(calccond), pointer :: cco
    type(riv_w_var), pointer :: rvwo
    type(rivvartime), pointer :: rvto
    type(csed), pointer :: csedo
    rvo => model%t_rivvar
    cco => model%t_calccond
    rvwo => model%t_rivwvar
    rvto => model%t_rivvartime
    csedo => model%t_csed

    call dealloc_all(rvo, rvwo, cco, rvto, csedo)
    end subroutine cleanup

    ! steps the heat model forward in time.
    subroutine advance_in_time(model)
    type (fastmech_model), intent (inout) :: model

    call solve_fm(model)
    !model%temperature = model%temperature_tmp
    !model%t = model%t + model%dt
    end subroutine advance_in_time

    !    subroutine solve_fm(model)
    !    implicit none
    !    type (fastmech_model), target, intent (inout) :: model
    !    integer :: i, j, jswi, j1, j2, j3
    !    integer :: n, n2
    !    integer :: solindex
    !    integer :: ibccount
    !    integer :: ier
    !    integer :: iter
    !    real(kind=mp) :: va, ta, tb, tc, td, te, tf
    !    real(kind=mp) :: ba, bb, bc, be, bf, bg, bot, bot2
    !    real(kind=mp) :: be1, be2, bf1, bf2
    !    real(kind=mp) :: ua, uu, vv
    !    real(kind=mp) :: tmplev, tmpvardt
    !    real(kind=mp) :: nodechange
    !    real(kind=mp) :: ustar
    !    real(kind=mp) :: area
    !    real(kind=mp) :: dinc
    !    real(kind=mp) :: dube
    !    real(kind=mp) :: sumqpdiff
    !    type(calccond), pointer :: cco
    !    type(rivvartime), pointer :: rvto
    !    type(rivvar), pointer :: rvo
    !    type(rivvarvert), pointer :: rvvo
    !    type(riv_w_var), pointer :: rvwo
    !    type(csed), pointer :: csedo
    !    ! pointers
    !    integer, pointer :: nct, soltype, varstagetype, vardischtype, roughnesstype
    !    integer, pointer :: interitm, maxinteritermult, errorcode, iterplout, iplinc
    !    logical, pointer :: calccsed, calcquasi3d, calcquasi3drs, calcsedauto
    !    real(kind=mp), pointer :: g, rho, vkc, dn
    !    real(kind=mp), pointer :: vardt, olddisch, q, oldstage, newstage, newdisch, ptime
    !    real(kind=mp), pointer :: ddisch, dstage, dsstage!, elevoffset
    !    real(kind=mp), dimension(:,:), pointer :: e, hl, u, v, eta, totcd, cd, cdv, taus, taun, rn
    !    real(kind=mp), dimension(:,:), pointer :: dude, dvde, vout, uout, eka
    !    real(kind=mp), dimension(:,:), pointer :: qs, qn, con, x, y, harea
    !    real(kind=mp), dimension(:), pointer :: hav, r, w, xo, yo, phirotation
    !    real(kind=mp), dimension(:,:), pointer :: dsu, dne, dnq, dsv, dse, dsq, dnu, dnv, dqe
    !    real(kind=mp), dimension(:,:), pointer :: dqsde, dqnde, qv, qu, dq, de
    !    real(kind=mp), dimension(:), pointer :: dea, eav, am, bm, ccm, dm, em, qp
    !    real(kind=mp), dimension(:), pointer :: qprms
    !    integer, dimension(:), pointer :: isw
    !
    !    integer, dimension(:,:), pointer :: ibc, icon
    !    integer, pointer :: nclusters, nwetnodes, noldwetnodes
    !    integer, pointer :: ns, nn, nm
    !    !calculation condition pointers
    !    integer, pointer :: drytype, hiterinterval, hiterstop, imod
    !    integer, pointer :: levtype, levchangeiter, levbegiter, levenditer
    !    integer, pointer :: debugstop, dbgtimestep, dbgiternum, transeqtype, nsteps
    !    integer, pointer :: vbcds, vbc
    !    logical, pointer :: iwetdry, iterationout, io_3doutput
    !    real(kind=mp), pointer :: evc, startlev, endlev, urelax, erelax, arelax
    !    real(kind=mp), pointer :: hmin
    !    ! for hot start
    !    integer, pointer :: i_re_flag_i, i_re_flag_o, n_rest, i_tmp_count
    !    real(8), dimension(:),  pointer :: opt_tmp
    !    character(len = strmax), dimension(:), pointer :: tmp_file_o, tmp_caption
    !    character(len=strmax), pointer ::tmp_file_i, tmp_dummy, tmp_pass
    !
    !    tmpvardt = 0
    !
    !    cco => model%t_calccond
    !    rvto => model%t_rivvartime
    !    rvo => model%t_rivvar
    !    rvvo => model%t_rivvarvert
    !    rvwo => model%t_rivwvar
    !    csedo => model%t_csed
    !
    !    dbgtimestep => model%t_calccond%dbgtimestep
    !    dbgiternum => model%t_calccond%dbgiternum
    !    debugstop => model%t_calccond%debugstop
    !    iwetdry => model%t_calccond%iwetdry
    !    iterationout => model%t_calccond%iterationout
    !    io_3doutput => model%t_calccond%io_3doutput
    !    iplinc => model%t_calccond%iplinc
    !
    !    nsteps => model%t_rivvar%nsteps
    !    nm => model%t_rivvar%nm
    !    hav => model%t_rivvar%hav
    !    dn => model%t_rivvar%dn
    !    phirotation => model%t_rivvar%phirotation
    !    harea => model%t_rivvar%harea
    !    x => model%t_rivvar%x
    !    y => model%t_rivvar%y
    !    w => model%t_rivvar%w
    !    xo => model%t_rivvar%xo
    !    yo => model%t_rivvar%yo
    !    qs => model%t_rivvar%qs
    !    qn => model%t_rivvar%qn
    !    con => model%t_rivvar%con
    !    r => model%t_rivvar%r
    !    rn => model%t_rivvar%rn
    !    g => model%t_rivvar%g
    !    rho => model%t_rivvar%rho
    !    vkc => model%t_rivvar%vkc
    !    vbcds => model%t_calccond%vbcds
    !    transeqtype => model%t_calccond%transeqtype
    !    hmin => model%t_calccond%hmin
    !    vbc => model%t_calccond%vbc
    !    i_re_flag_i => model%t_calccond%i_re_flag_i
    !    i_re_flag_o => model%t_calccond%i_re_flag_o
    !    n_rest => model%t_calccond%n_rest
    !    i_tmp_count => model%t_calccond%i_tmp_count
    !    tmp_file_o => model%t_calccond%tmp_file_o
    !    tmp_caption => model%t_calccond%tmp_caption
    !    tmp_file_i => model%t_calccond%tmp_file_i
    !    tmp_dummy => model%t_calccond%tmp_dummy
    !    tmp_pass => model%t_calccond%tmp_pass
    !    opt_tmp => model%t_calccond%opt_tmp
    !
    !    iwetdry => model%t_calccond%iwetdry
    !    urelax => model%t_calccond%urelax
    !    erelax => model%t_calccond%erelax
    !    arelax => model%t_calccond%arelax
    !    levtype => model%t_calccond%levtype
    !    levchangeiter => model%t_calccond%levchangeiter
    !    levbegiter => model%t_calccond%levbegiter
    !    levenditer => model%t_calccond%levenditer
    !    evc => model%t_calccond%evc
    !    startlev => model%t_calccond%startlev
    !    endlev => model%t_calccond%endlev
    !
    !    imod => model%t_calccond%imod
    !    drytype => model%t_calccond%drytype
    !    hiterinterval => model%t_calccond%hiterinterval
    !    hiterstop => model%t_calccond%hiterstop
    !    qprms => model%t_rivwvar%qprms
    !    isw => model%t_rivwvar%isw
    !    dea => model%t_rivwvar%dea
    !    eav => model%t_rivwvar%eav
    !    am => model%t_rivwvar%am
    !    bm => model%t_rivwvar%bm
    !    ccm => model%t_rivwvar%ccm
    !    dm => model%t_rivwvar%dm
    !    em => model%t_rivwvar%em
    !    qp => model%t_rivwvar%qp
    !
    !    dqsde => model%t_rivwvar%dqsde
    !    dqnde => model%t_rivwvar%dqnde
    !    qv => model%t_rivwvar%qv
    !    qu => model%t_rivwvar%qu
    !    dq => model%t_rivwvar%dq
    !    de => model%t_rivwvar%de
    !    dnu => model%t_rivwvar%dnu
    !
    !    dnv => model%t_rivwvar%dnv
    !    dsu => model%t_rivwvar%dsu
    !    dne => model%t_rivwvar%dne
    !    dnq => model%t_rivwvar%dnq
    !    dsv => model%t_rivwvar%dsv
    !    dse => model%t_rivwvar%dse
    !    dsq => model%t_rivwvar%dsq
    !    dqe => model%t_rivwvar%dqe
    !
    !    dude => model%t_rivwvar%dude
    !    dvde => model%t_rivwvar%dvde
    !    eka => model%t_rivwvar%eka
    !    taus => model%t_rivvar%taus
    !    taun => model%t_rivvar%taun
    !
    !    totcd => model%t_rivvar%totcd
    !    cd => model%t_rivvar%cd
    !    cdv => model%t_rivvar%cdv
    !    ns => model%t_rivvar%ns
    !    nn => model%t_rivvar%nn
    !    errorcode => model%t_rivvar%errorcode
    !    ptime => model%t_rivvar%ptime
    !    nct => model%t_rivvar%nct
    !    vardt => model%t_rivvar%vardt
    !    calccsed => model%t_calccond%calccsed
    !    iterplout => model%t_calccond%iterplout
    !    calcquasi3d => model%t_calccond%calcquasi3d
    !    calcquasi3drs => model%t_calccond%calcquasi3drs
    !    calcsedauto => model%t_calccond%calcsedauto
    !    !elevoffset => model%t_rivvar%elevoffset
    !    olddisch => model%t_rivvar%olddisch
    !    oldstage => model%t_rivvar%oldstage
    !    newstage => model%t_rivvar%newstage
    !    newdisch => model%t_rivvar%newdisch
    !    ddisch => model%t_rivvar%ddisch
    !    dstage => model%t_rivvar%dstage
    !    dsstage => model%t_rivvar%dsstage
    !    varstagetype => model%t_rivvartime%varstagetype
    !    vardischtype => model%t_rivvartime%vardischtype
    !    roughnesstype => model%t_calccond%roughnesstype
    !    e => model%t_rivvar%e
    !    hl => model%t_rivvar%hl
    !    u => model%t_rivvar%u
    !    v => model%t_rivvar%v
    !    eta => model%t_rivvar%eta
    !    ibc => model%t_rivvar%ibc
    !    icon => model%t_rivvar%icon
    !    nclusters => model%t_rivvar%nclusters
    !    nwetnodes => model%t_rivvar%nwetnodes
    !    noldwetnodes => model%t_rivvar%noldwetnodes
    !    interitm => model%t_calccond%interitm
    !    maxinteritermult => model%t_calccond%maxinteritermult
    !    w => model%t_rivvar%w
    !    q => model%t_calccond%q
    !
    !    nct = nct+1
    !
    !    if(nct > 0.and.cco%soltype == 0) then
    !        return
    !    endif
    !    if(nct > 0.and.model%t_end.eq.rvto%vardischstarttime) then
    !        return
    !    endif
    !
    !    if(nct == 0) then
    !        vardt = model%t_calccond%fmdt
    !    else
    !        if(vardt > model%t_calccond%fmdt) then
    !            vardt = model%t_calccond%fmdt
    !        endif
    !    endif
    !
    !    if(calccsed.eqv..false.) then
    !        vardt = model%t_calccond%fmdt
    !    endif
    !
    !    olddisch = q
    !    oldstage = newstage
    !
    !    if(nct > 0) then
    !        model%t = model%t+vardt
    !        if(model%t > model%t_end) then
    !            write(*,*) 'model time %f is greater than end time %f', model%t, model%t_end
    !            return
    !        endif
    !    endif
    !
    !    if(vardischtype == 1) then !update discharge
    !        call getinterptimeseriesvalue(rvto, 1, model%t, newdisch)
    !        q = newdisch*1e6
    !    endif
    !    ddisch = q-olddisch
    !
    !    if(varstagetype == 1) then !stage time series
    !        call getinterpratingcurvevalue(rvto, 1, model%t, newstage)
    !        !            newstage = newstage * 100.
    !    else if (varstagetype == 2) then !stage rating curve
    !        call getinterpratingcurvevalue(rvto, 1, newdisch, newstage)
    !        !            newstage = newstage * 100.
    !    endif
    !
    !    dstage = newstage-oldstage
    !
    !    !        if(cdtype == 1) then !variable roughness
    !    !            call setroughness(ns, nn, roughnesstype, newdisch, cd, znaught)
    !    !        endif
    !
    !    if(vardischtype == 1) then
    !        !!  adjust water surface elevation on all wet nodes by dstage !!
    !        !!  check for minimum depth and turn nodes on or off.         !!
    !        call vardischeinit(rvo, cco, e, hl, u, v, eta, ibc, dstage, dsstage)
    !    endif
    !
    !    !! turn nodes on - check re-wetting !!
    !    call updatewetting(rvo, cco, e,hl,u,v,eta, ibc, dsstage)
    !
    !    !! set upstream boundary velocity  !!
    !    if(model%t_rivvartime%vardischtype == 1) then
    !        call vardischuinit(rvo, cco, u, v, w, eta, e, ibc, q, dstage)
    !    endif
    !
    !    !! finds isolated pockets of wet nodes and deletes them
    !    call label_clusters(rvo, -1, ibc, nclusters, icon)
    !    call delete_clusters(rvo, ibc, icon, nclusters, nwetnodes)
    !
    !    !!  check for ibc == 6 or ibc == 4  !!
    !    call updateibc(rvo)
    !
    !    !!  if grid extension then reset extension topography !!
    !    !        if(nct > 1) then
    !    !            call resetgridextension()
    !    !        endif
    !    if(roughnesstype == 1) then !roughness as z0
    !        !		call zotocdtwo() !calculate cd from log profile
    !        call z0tocd(rvo, rvvo, cco) !calculate cd from 2-part profile
    !    endif
    !
    !
    !    do i = 1,10
    !        write(*,*)
    !    enddo
    !
    !    write(*,*)
    !    write(*,*) 'variabledt', tmpvardt, 'dt', model%t_calccond%fmdt
    !    write(*,*) 'time step', nct, 'time', model%t, 'printtime', ptime
    !    write(*,*)
    !    if(calcquasi3d) then
    !        write(*,*) 'quasi-3d calculation is on'
    !        if(calcquasi3drs) then
    !            write(*,*) 'streamline curvature is on'
    !        endif
    !    else
    !        write(*,*) 'quasi-3d calculation is off'
    !    endif
    !    write(*,*)
    !    if(calccsed) then
    !        write(*,*) 'sediment-transport is on'
    !        if(calcsedauto) then
    !            write(*,*) 'time-stepping method is automatic'
    !        else
    !            write(*,*) 'time-stepping method is fixed'
    !        endif
    !        write(*,*) 'old wet nodes: ', noldwetnodes, ' new wet nodes: ', nwetnodes
    !
    !    else
    !        write(*,*) 'sediment-transport is off'
    !    endif
    !
    !    write(*,*)
    !    write(*,*) 'discharge', q/1e6, 'change discharge', ddisch/1e6
    !    write(*,*) 'stage', newstage/100., 'change stage', dstage/100.
    !    write(*,*)
    !
    !    if(nct.ge.1) then
    !        model%itermax=interitm
    !        nodechange = abs((real(noldwetnodes)-real(nwetnodes)))
    !        model%itermax = real(nodechange*interitm)
    !        if(model%itermax < interitm)then
    !            model%itermax = interitm
    !        else if (model%itermax >= 5*interitm) then
    !            model%itermax = maxinteritermult*interitm
    !        endif
    !        noldwetnodes = nwetnodes
    !    endif
    !
    !    if(model%itermax /= 0) then
    !        iter_loop: do iter = 1,model%itermax
    !            !       check if user cancelled simulation and if so exit
    !            call iric_check_cancel_f(errorcode)
    !            if(errorcode.eq.1) then
    !                call dealloc_all(rvo, rvwo, cco, rvto, csedo)
    !                return
    !            endif
    !
    !            if(iwetdry.and.mod(iter,imod).eq.0.and.model%t_rivvartime%vardischtype == 0.and.iter.lt.hiterstop.and.nct.lt.1)then
    !                call    updatewetting(rvo, cco, e, hl, u, v, eta, ibc, dsstage)
    !            endif
    !
    !            call updateibc(rvo)
    !
    !            do i = 1,ns
    !                do j = 1,nn
    !                    if(ibc(i,j).eq.1) then
    !                        ustar=(totcd(i,j)*(u(i,j)**2.))**.5
    !                    else
    !                        ustar=(totcd(i,j)*(u(i,j)**2.+v(i,j)**2.))**.5
    !                    endif
    !
    !                    !eka(i,j)=.005*200.*200.+vkc*hl(i,j)*ustar/6.
    !                    if(levtype == 0) then
    !                        eka(i,j) = evc+vkc*hl(i,j)*ustar/6.
    !                    else if(levtype.ne.0.and.nct.eq.0) then
    !                        if(iter < levbegiter) then
    !                            tmplev = startlev
    !                        elseif (iter >= levbegiter .and. iter < levenditer) then
    !                            tmplev = &
    !                                (((startlev-endlev)/(levbegiter-levenditer))*&
    !                                (iter-levbegiter))+startlev
    !                        elseif (iter > levenditer) then
    !                            tmplev = endlev
    !                        endif
    !                        eka(i,j) = tmplev+vkc*hl(i,j)*ustar/6.
    !
    !                    else
    !                        eka(i,j) = endlev+vkc*hl(i,j)*ustar/6.
    !                    endif
    !                enddo
    !            enddo
    !
    !
    !            do i=2,ns
    !                do jswi=1,2
    !                    if(jswi.eq.1) then
    !                        j1=1
    !                        j2=nn
    !                        j3=1
    !                    else
    !                        j1=nn
    !                        j2=1
    !                        j3=-1
    !                    endif
    !                    do j=j1,j2,j3
    !                        if(ibc(i,j).eq.2) then
    !                            va=(v(i-1,j)+v(i,j))/4.
    !                        else if(j.lt.nn) then
    !                            va=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4.
    !                        endif
    !                        if(ibc(i,j).eq.4.or.ibc(i,j).eq.0) then
    !                            dude(i,j)=0.
    !                            u(i,j)=0.
    !                            cycle
    !                        endif
    !
    !                        if(u(i,j).ge.0.or.i.eq.ns) then
    !                            ba=u(i,j)/dsu(i,j)
    !                            ta=ba*u(i-1,j)
    !                        else
    !                            ba=-1.*u(i,j)/dsu(i+1,j)
    !                            ta=ba*u(i+1,j)
    !                        endif
    !
    !                        if((va.ge.0.and.ibc(i,j).ne.1).or.ibc(i,j).eq.2) then
    !                            bb=va/dnu(i,j)
    !                            tb=bb*u(i,j-1)
    !                        else
    !                            bb=-1.*va/dnu(i,j+1)
    !                            tb=bb*u(i,j+1)
    !                        endif
    !
    !                        bc=-1.*va/(rn(i,j)*r(i))
    !                        if(ibc(i-1,j).eq.0) then
    !                            td = 0.
    !                        else
    !                            td=-1.*g*(e(i,j)-e(i-1,j))/dse(i,j)
    !                        endif
    !
    !                        if(i.eq.ns) then
    !                            be=(2./(dsu(i,j)))*eka(i,j)/dsu(i,j)
    !                            te=be*u(i-1,j)
    !                        else
    !                            be=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)/dsu(i,j))+&
    !                                (eka(i+1,j)/dsu(i+1,j)))
    !                            te=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)*u(i-1,j)/dsu(i,j))+&
    !                                (eka(i+1,j)*u(i+1,j)/dsu(i+1,j)))
    !                        endif
    !
    !                        if(ibc(i,j).eq.2) then
    !                            bf=(2./dnu(i,j))*(totcd(i,j)*abs(u(i,j))+eka(i,j)/dnu(i,j))
    !                            tf=(2./dnu(i,j))*(eka(i,j)*u(i,j-1)/dnu(i,j))
    !                        elseif(ibc(i,j).eq.1) then
    !                            bf=(4./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+&
    !                                totcd(i,j)*abs(u(i,j)))
    !                            tf=(4./(dnu(i,j)+dnu(i,j+1)))*(eka(i,j+1)*u(i,j+1)/dnu(i,j+1))
    !                        else
    !                            bf=-1.*eka(i,j)/(dnu(i,j)*rn(i,j)*r(i))+(2./(dnu(i,j)+&
    !                                dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+(eka(i,j)/dnu(i,j)))
    !                            tf=-1.*eka(i,j)*u(i,j-1)/(dnu(i,j-1)*rn(i,j-1)*r(i))+&
    !                                (2./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)*u(i,j+1)/dnu(i,j+1))+&
    !                                (eka(i,j)*u(i,j-1)/dnu(i,j)))
    !                        endif
    !
    !                        bot2=rho*(hl(i,j)+hl(i-1,j))
    !
    !                        if(bot2.eq.0.) then
    !                            bg=10.**12.
    !                        else
    !                            bg=2.*(totcd(i,j)*(u(i,j)**2.+va**2.)**.5)/bot2
    !                        endif
    !
    !                        bot=ba+bb+bc+be+bf+bg
    !
    !                        if(abs(bot).gt..000001) then
    !                            u(i,j)=u(i,j)*(1.-urelax)+urelax*(ta+tb+td+te+tf)/bot
    !                            if(i.eq.ns.and.u(i,j).lt.0) then
    !                                if(vbcds .ne. 0) then
    !                                    u(i,j) = 0.
    !                                endif
    !                            endif
    !                            dude(i,j)=-1.*g/(dse(i,j)*bot)
    !                        endif
    !                    enddo
    !                enddo
    !            enddo
    !
    !            do i=2,ns
    !                do jswi=1,2
    !                    if(jswi.eq.1) then
    !                        j1=1
    !                        j2=nn
    !                        j3=1
    !                    else
    !                        j1=nn
    !                        j2=1
    !                        j3=-1
    !                    endif
    !                    do j=j1,j2,j3
    !                        if(ibc(i,j).eq.0.or.ibc(i,j).eq.1) then
    !                            v(i,j)=0.
    !                            dvde(i,j)=0.
    !                            cycle ! will need to make sure this change is correct
    !                        endif
    !                        if(i.eq.ns) then
    !                            ua=(u(i,j)+u(i,j-1))/2.
    !                        else
    !                            ua=(u(i,j)+u(i,j-1)+u(i+1,j)+u(i+1,j-1))/4.
    !                        endif
    !                        if(ua.ge.0.or.i.eq.ns) then
    !                            ba=ua/dsv(i,j)
    !                            ta=ba*v(i-1,j)
    !                        else
    !                            ba=-1.*ua/dsv(i+1,j)
    !                            ta=ba*v(i+1,j)
    !                        endif
    !                        if(v(i,j).ge.0.or.ibc(i,j).eq.2) then
    !                            bb=v(i,j)/dnv(i,j)
    !                            tb=bb*v(i,j-1)
    !                        else
    !                            bb=-1.*v(i,j)/dnv(i,j+1)
    !                            tb=bb*v(i,j+1)
    !                        endif
    !
    !                        tc=-1.*(ua**2.)/(rn(i,j)*r(i))
    !                        if(ibc(i,j-1).eq.0) then
    !                            td = 0.
    !                        else
    !                            td=-1.*g*(e(i,j)-e(i,j-1))/dne(i,j)
    !                        endif
    !                        if(i.eq.ns) then
    !                            be=(1./dsv(i,j))*eka(i,j)/dsv(i,j)
    !                            te=be*v(i-1,j)
    !                        else
    !                            be1=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i+1,j)/dsv(i+1,j)
    !                            be2=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i,j)/dsv(i,j)
    !                            be=be1+be2
    !                            te=be1*v(i+1,j)+be2*v(i-1,j)
    !                        endif
    !                        if(ibc(i,j).eq.2) then
    !                            bf1=(1./dnv(i,j))*eka(i,j)/dnv(i,j)
    !                            bf2=bf1
    !                            bf=bf1+bf2
    !                            tf=bf1*v(i,j-1)
    !                        else
    !                            bf1=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j+1)/dnv(i,j+1)
    !                            bf2=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j)/dnv(i,j)
    !                            bf=bf1+bf2
    !                            tf=bf1*v(i,j+1)+bf2*v(i,j-1)
    !                        endif
    !                        bot2=rho*(hl(i,j)+hl(i,j-1))
    !                        if(bot2.eq.0.) then
    !                            bg=10.**12.
    !                        else
    !                            bg=2.*(totcd(i,j)*(ua**2.+v(i,j)**2.)**.5)/bot2
    !                        endif
    !                        bot=ba+bb+be+bf+bg
    !                        if(abs(bot).gt..000001) then
    !                            if(vbcds .eq. 0) then
    !                                v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
    !                            else
    !                                if(i.ne.ns) then
    !                                    v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
    !                                endif
    !                            endif
    !                            dvde(i,j)=-1.*g/(dn*bot)
    !                        endif
    !                    enddo
    !                enddo
    !            enddo
    !
    !            do i=1,ns
    !                do j=1,nn
    !                    if(i.eq.1) then
    !                        area=dnq(i,j)*hl(i,j)
    !                    else
    !                        area=dnq(i,j)*(hl(i,j)+hl(i-1,j))/2.
    !                    endif
    !                    if(j.eq.1.or.j.eq.nn) area=area/2.
    !                    qu(i,j)=u(i,j)*area
    !                    dqsde(i,j)=dude(i,j)*area
    !                    if(j.eq.1) then
    !                        area=dsq(i,j)*hl(i,j)
    !                    else
    !                        !        area=.5*(dsq(i,j-1)+dsq(i,j))*(hl(i,j)+hl(i,j-1))/2.
    !                        area=dsq(i,j)*(hl(i,j)+hl(i,j-1))/2.
    !
    !                    endif
    !                    qv(i,j)=v(i,j)*area
    !                    dqnde(i,j)=dvde(i,j)*area
    !                enddo
    !            enddo
    !
    !            do i=1,ns-1
    !                do j=1,nn
    !                    if(j.eq.nn) then
    !                        dq(i,j)=qu(i+1,j)-qu(i,j)-qv(i,j)
    !                    else
    !                        dq(i,j)=qu(i+1,j)-qu(i,j)+qv(i,j+1)-qv(i,j)
    !                    endif
    !                enddo
    !            enddo
    !
    !            do i=1,ns-1
    !                do j=1,nn
    !                    if(j.eq.nn) then
    !                        dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)
    !                    else
    !                        dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)+dqnde(i,j+1)
    !                    endif
    !                    if(abs(dqe(i,j)).lt..000001) then
    !                        de(i,j)=0.
    !                    else
    !                        de(i,j)=dq(i,j)/dqe(i,j)
    !                    endif
    !                enddo
    !            enddo
    !
    !            do i=ns-1,2,-1
    !                n=0.
    !                do j=1,nn
    !                    !c       if(abs(dqe(i,j)).gt..000001) then
    !                    if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
    !                        n=n+1
    !                        if(j.ne.1) then
    !                            am(n)=-1.*dqnde(i,j)
    !                        endif
    !                        bm(n)=dqe(i,j)
    !
    !                        if(j.eq.nn) then
    !                            ccm(n)=0.
    !                        else
    !                            ccm(n)=-1.*dqnde(i,j+1)
    !                        endif
    !                        dm(n)=dq(i,j)+dqsde(i+1,j)*de(i+1,j)+dqsde(i,j)*de(i-1,j)
    !                        isw(j)=1
    !                    else
    !                        isw(j)=0
    !                    endif
    !                enddo
    !                if(n.ge.1) then
    !                    call tridag(am,bm,ccm,dm,em,n,errorcode)
    !                    if(errorcode.eq.-1) then
    !                        !call write_error_cgns(str_in)
    !                        call dealloc_common2d(rvo)
    !                        if(vbc == 1) then
    !                            call dealloc_velbc(cco)
    !                        endif
    !
    !                        if(calcquasi3d) then
    !                            call dealloc_common3d(rvo)
    !                            if(transeqtype == 2) then
    !                                !call dealloc_csed3d_dt()
    !                            endif
    !                        endif
    !                        call dealloc_working(rvwo)
    !                        write(6,*) 'error at dqs', 'i = ', i
    !
    !                        !pause
    !                        return
    !                    endif
    !                endif
    !                n2=0
    !                do j=1,nn
    !                    if(isw(j).eq.1) then
    !                        n2=n2+1
    !                        de(i,j)=em(n2)
    !                    else
    !                        de(i,j)=0.
    !                    endif
    !                enddo
    !            enddo
    !
    !            do j=1,nn
    !                n=0
    !                do i=ns-1,2,-1
    !                    if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
    !                        n=n+1
    !                        am(n)=-1.*dqsde(i,j)
    !                        bm(n)=dqe(i,j)
    !                        ccm(n)=-1.*dqsde(i+1,j)
    !                        if(j.eq.1) then
    !                            dm(n)=dq(i,j)+dqnde(i,j+1)*de(i,j+1)
    !                        elseif(j.eq.nn) then
    !                            dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)
    !                        else
    !                            dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)+dqnde(i,j+1)*de(i,j+1)
    !                        endif
    !                        isw(i)=1
    !                    else
    !                        isw(i)=0
    !                    endif
    !                enddo
    !                if(n.ge.1) then
    !                    call tridag(am,bm,ccm,dm,em,n,errorcode)
    !                    if(errorcode.eq.-1) then
    !                        !call write_error_cgns(str_in)
    !                        call dealloc_common2d(rvo)
    !                        if(vbc == 1) then
    !                            call dealloc_velbc(cco)
    !                        endif
    !
    !                        if(calcquasi3d) then
    !                            call dealloc_common3d(rvo)
    !                            if(transeqtype == 2) then
    !                                !call dealloc_csed3d_dt()
    !                            endif
    !                        endif
    !                        call dealloc_working(rvwo)
    !                        write(6,*) 'error at dqn'
    !
    !                        !pause
    !                        return
    !                    endif
    !                endif
    !                n2=0
    !                do i=ns-1,2,-1
    !                    if(isw(i).eq.1) then
    !                        n2=n2+1
    !                        de(i,j)=em(n2)
    !                    else
    !                        de(i,j)=0.
    !                    endif
    !                enddo
    !            enddo
    !
    !            do i=1,ns
    !                qp(i)=0.
    !                eav(i)=0.
    !                ibccount = 0
    !                do j=1,nn
    !                    !        qp(i)=qp(i)+qu(i,j)
    !                    !        e(i,j)=e(i,j)+erelax*de(i,j)
    !                    if(ibc(i,j).ne.0) then
    !                        e(i,j)=e(i,j)+erelax*de(i,j)
    !                        ibccount = ibccount+1
    !                        eav(i)=eav(i)+e(i,j)
    !                    endif
    !                    qp(i)=qp(i)+qu(i,j)
    !                    hl(i,j)=e(i,j)-eta(i,j)
    !                    if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
    !                        if(drytype.eq.0) then
    !                            if(ibc(i,j).ne.0) then
    !                                u(i,j) = 0
    !                                v(i,j) = 0
    !                            endif
    !                            ibc(i,j)=0
    !                        endif
    !                        hl(i,j)=hmin
    !                    endif
    !                    !unsure of the effect of the code block below rmcd 12/23/05
    !                    !change rmcd made sure i>1 4/30/10
    !                    if(i.gt.1) then
    !                        if(abs(u(i,j)).lt..00001) then
    !                            u(i,j)=0.
    !                        endif
    !                        if(abs(v(i,j)).lt..00001) then
    !                            v(i,j)=0.
    !                        endif
    !                    endif
    !                    !!!!!!!!!!!!
    !                enddo
    !                eav(i)=eav(i)/ibccount
    !                do j=1,nn
    !                    if(ibc(i,j).eq.0)then
    !                        e(i,j) = eav(i)
    !                    endif
    !                enddo
    !            enddo
    !
    !            do  j=1,nn
    !                u(1,j)=u(1,j)*q/qp(1)
    !            enddo
    !
    !            dea(ns)=0.
    !            do  i=ns-1,1,-1
    !                dinc=arelax*erelax*abs(eav(i)-eav(i+1))*(1.-qp(i+1)/q)
    !                dea(i)=dea(i+1)+dinc
    !                do  j=1,nn
    !                    !! code change by rmcd !!
    !                    if(ibc(i,j).ne.0)then
    !600                     e(i,j)=e(i,j)+dea(i)
    !                        !                else
    !                        !                    e(i,j) = eav(i)+dea(i) !! added so that e(i,j) in dry nodes follows
    !                        !                                           !! average eav(i) during time-dependent runs
    !                    endif
    !                    !! end code change by rmcd !!
    !                enddo
    !            enddo
    !
    !            sumqpdiff = 0.
    !            do i=1,ns
    !                sumqpdiff = sumqpdiff + (((qp(i)-q)/q*100.))**2.
    !            enddo
    !
    !            qprms(iter) = sqrt(sumqpdiff/ns)
    !            !write(*,*) 'after qprms'
    !            if(isnan(qprms(iter))) then
    !                errorcode=-1
    !                !call write_error_cgns(str_in)
    !                !			call dealloc_all()
    !                call dealloc_common2d(rvo)
    !                if(vbc == 1) then
    !                    call dealloc_velbc(cco)
    !                endif
    !
    !                if(calcquasi3d) then
    !                    call dealloc_common3d(rvo)
    !                    if(transeqtype == 2) then
    !                        !call dealloc_csed3d_dt()
    !                    endif
    !                endif
    !                call dealloc_working(rvwo)
    !                write(6,*) 'error at qprms'
    !                !pause
    !                return
    !            endif
    !
    !            if(nct == dbgtimestep) then
    !                if(debugstop.eq.1.and.iter == dbgiternum) then
    !                    solindex = solindex+1
    !                    call cg_iric_write_sol_time_f(model%t, ier)
    !
    !                    call write_cgns2(rvo, cco, rvwo, model%t, q)
    !                    !			        call write_timestep_cgns(str_in, nct, tottime)
    !                    !			        call write_timeiter_cgns(str_in)
    !                    call dealloc_all(rvo, rvwo, cco, rvto, csedo)
    !                    !pause
    !
    !                    return
    !                endif
    !            endif
    !            if(levtype == 0) then
    !                write(*,*) 'iteration: ', iter,' mean error on discharge:', qprms(iter)
    !            else
    !                write(*,*)'iteration:', iter, ' lev:', tmplev,' mean error on discharge: ', qprms(iter)
    !            endif
    !
    !            if(iterationout.and.mod(iter,iterplout).eq.0) then
    !
    !                do i = 1,ns
    !                    do j = 1,nn
    !                        if(ibc(i,j).eq.0) then
    !                            v(i,j) = 0.
    !                            u(i,j) = 0.
    !                        endif
    !                    enddo
    !                enddo
    !
    !                do  i=1,ns
    !                    do  j=1,nn
    !                        !                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
    !                        if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
    !                            u(i,j)=0.
    !                            vout(i,j)=0.
    !                            !    go to 888
    !                            !endif
    !                        elseif(i.eq.1) then
    !                            vout(i,j)=v(i,j)
    !                            !    go to 888
    !                            !endif
    !
    !                        elseif(j.eq.nn) then
    !                            vout(i,j)=(v(i-1,j)+v(i,j))/4.
    !                        elseif(j.eq.1) then
    !                            vout(i,j) = 0.
    !                        elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
    !                            vout(i,j)=v(i,j)
    !                        elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
    !                            vout(i,j)=v(i,j+1)
    !                        else
    !                            vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
    !                        endif
    !
    !888                     if(i.eq.1) then
    !                            dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
    !                        else
    !                            dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
    !                        endif
    !                        taus(i,j)=dube*u(i,j)
    !                        taun(i,j)=dube*vout(i,j)
    !                    enddo
    !                enddo
    !
    !                solindex = solindex+1
    !                call cg_iric_write_sol_time_f(model%t, ier)
    !                call write_cgns2(rvo, cco, rvwo, model%t, q)
    !
    !            endif
    !
    !        end do iter_loop
    !
    !        if(i_re_flag_o.eq.1.and.nct.eq.0) then
    !            if(i_tmp_count <= n_rest) then
    !                if(iter.eq.opt_tmp(i_tmp_count)) then
    !                    tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
    !                    open(502,file=tmp_file_o(i_tmp_count) &
    !                        ,status='unknown',form='unformatted')
    !                    !
    !                    write(502) model%t,solindex,model%t_calccond%fmdt
    !                    write(502) ns,nn
    !                    !
    !                    write(502) ((eta(i,j),i=1,ns),j=1,nn)
    !                    write(502) ((u(i,j),i=1,ns),j=1,nn)
    !                    write(502) ((v(i,j),i=1,ns),j=1,nn)
    !                    write(502) ((ibc(i,j),i=1,ns),j=1,nn)
    !                    write(502) ((e(i,j),i=1,ns),j=1,nn)
    !                    write(502) ((hl(i,j),i=1,ns),j=1,nn)
    !                    write(502) (hav(i),i=1,ns)
    !
    !                    close(502)
    !
    !                    i_tmp_count = i_tmp_count +1
    !                endif
    !            endif
    !        endif
    !
    !64      format(25f7.3)
    !        !66      format(a15,i6,a25,f15.5)
    !66      format(a, a, i5,a, f15.5)
    !67      format(a, a, i5,a, f15.5,a, f15.5)
    !68      format(a)
    !75      format(25f7.0)
    !85      format(25f6.4)
    !86      format(25f8.3)
    !84      format(3i5,f6.2,2f12.1,2g12.4)
    !
    !        do i = 1,ns
    !            do j = 1,nn
    !                if(ibc(i,j).eq.0) then
    !                    v(i,j) = 0.
    !                    u(i,j) = 0.
    !                endif
    !            enddo
    !        enddo
    !    else
    !
    !        !601 continue
    !        do i=1,ns
    !            do j=1,nn
    !                !          if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
    !                if(ibc(i,j).ne.-1) then
    !                    u(i,j)=0.
    !                    vout(i,j)=0.
    !                    go to 88
    !                endif
    !                if(i.eq.1) then
    !                    vout(i,j)=v(i,j)
    !                    go to 88
    !                endif
    !
    !                if(j.eq.nn) then
    !                    vout(i,j)=(v(i-1,j)+v(i,j))/4.
    !                elseif(j.eq.1) then
    !                    vout(i,j) = 0.
    !                elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
    !                    vout(i,j)=v(i,j)
    !                elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
    !                    vout(i,j)=v(i,j+1)
    !                else
    !                    vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
    !                endif
    !88              if(i.eq.1) then
    !                    dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
    !                else
    !                    dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
    !                endif
    !                taus(i,j)=dube*u(i,j)
    !                taun(i,j)=dube*vout(i,j)
    !            enddo
    !        enddo
    !
    !        !c        call vert
    !        if(calcquasi3d) then
    !            call vert(rvo, rvvo, cco)
    !
    !        endif
    !        if(calccsed) then
    !            if(transeqtype == 2) then
    !                !call csed_dt(vardt,nct,nsteps, tmpvardt)
    !            else
    !                call csed2(csedo, rvo, cco, vardt,nct,nsteps, tmpvardt)
    !            endif
    !        else
    !            call stressdiv(rvo, cco, ibc, qs, qn, taus, taun, con, rn, r)
    !        endif
    !        if(calcsedauto) then
    !            vardt = tmpvardt
    !        endif
    !
    !4000    format(6f10.2)
    !        if(nct == 0) then
    !            solindex = solindex+1
    !            call calc_area(ns, nn, phirotation, x, y, xo, yo, nm, dn, harea)
    !            call cg_iric_write_sol_time_f(model%t, ier)
    !            call write_cgns2(rvo, cco, rvwo, model%t, q)
    !            !                if(calcquasi3d.and.io_3doutput) then
    !            !                    call write_cgns3d_grid()
    !            !                endif
    !            if(calcquasi3d.and.io_3doutput) then
    !                call write_cgns3d_grid(rvo)
    !                !                    call write_cgns3d_solgrid()
    !                call write_cgns3d_fixedbed(rvo, rvto, solindex, model%t, q)
    !            endif
    !
    !        else
    !            if(model%t >= ptime)then
    !                solindex = solindex+1
    !                call calc_area(ns, nn, phirotation, x, y, xo, yo, nm, dn, harea)
    !                call cg_iric_write_sol_time_f(model%t, ier)
    !                call write_cgns2(rvo, cco, rvwo, model%t, q)
    !                if(calcquasi3d.and.io_3doutput) then
    !                    !                        call write_cgns3d_solgrid()
    !                    call write_cgns3d_fixedbed(rvo, rvto, solindex, model%t, q)
    !                    !                    elseif(calcquasi3d.and.io_3doutput.and.calccsed) then
    !                    !!                        call write_cgns3d_moveablebed(tottime, q)
    !                endif
    !                ptime = ptime+(iplinc*model%t_calccond%fmdt)
    !            endif
    !            if(i_re_flag_o.eq.1.and.nct.ne.0) then
    !                if(i_tmp_count <= n_rest) then
    !                    if(model%t.ge.opt_tmp(i_tmp_count)) then
    !                        tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
    !                        open(502,file=tmp_file_o(i_tmp_count) &
    !                            ,status='unknown',form='unformatted')
    !                        !
    !                        write(502) model%t,solindex,model%t_calccond%fmdt
    !                        write(502) ns,nn
    !                        !
    !                        write(502) ((eta(i,j),i=1,ns),j=1,nn)
    !                        write(502) ((u(i,j),i=1,ns),j=1,nn)
    !                        write(502) ((v(i,j),i=1,ns),j=1,nn)
    !                        write(502) ((ibc(i,j),i=1,ns),j=1,nn)
    !                        write(502) ((e(i,j),i=1,ns),j=1,nn)
    !                        write(502) ((hl(i,j),i=1,ns),j=1,nn)
    !                        write(502) (hav(i),i=1,ns)
    !
    !                        close(502)
    !
    !                        i_tmp_count = i_tmp_count +1
    !                    endif
    !                endif
    !            endif
    !        endif
    !    endif
    !
    !    end subroutine solve_fm

    subroutine  solve_fm(model)
    implicit none
    type (fastmech_model), target, intent (inout) :: model
    integer :: i, j, jswi, j1, j2, j3
    integer :: n, n2
    integer :: solindex
    integer :: ibccount
    integer :: ier
    integer :: iter
    real(kind=mp) :: va, ta, tb, tc, td, te, tf
    real(kind=mp) :: ba, bb, bc, be, bf, bg, bot, bot2
    real(kind=mp) :: be1, be2, bf1, bf2
    real(kind=mp) :: ua, uu, vv
    real(kind=mp) :: tmplev, tmpvardt
    real(kind=mp) :: nodechange
    real(kind=mp) :: ustar
    real(kind=mp) :: area
    real(kind=mp) :: dinc
    real(kind=mp) :: dube
    real(kind=mp) :: sumqpdiff
    type(calccond), pointer :: cco
    type(rivvartime), pointer :: rvto
    type(rivvar), pointer :: rvo
    type(rivvarvert), pointer :: rvvo
    type(riv_w_var), pointer :: rvwo
    type(csed), pointer :: csedo

    tmpvardt = 0

    cco => model%t_calccond
    rvto => model%t_rivvartime
    rvo => model%t_rivvar
    rvvo => model%t_rivvarvert
    rvwo => model%t_rivwvar
    csedo => model%t_csed
    associate(  soltype => model%t_calccond%soltype, &
        dbgtimestep => model%t_calccond%dbgtimestep, &
        dbgiternum => model%t_calccond%dbgiternum, &
        debugstop => model%t_calccond%debugstop, &
        iwetdry => model%t_calccond%iwetdry, &
        iterationout => model%t_calccond%iterationout, &
        io_3doutput => model%t_calccond%io_3doutput, &
        iplinc => model%t_calccond%iplinc, &

        nsteps => model%t_rivvar%nsteps, &
        nm => model%t_rivvar%nm, &
        hav => model%t_rivvar%hav, &
        dn => model%t_rivvar%dn, &
        phirotation => model%t_rivvar%phirotation, &
        harea => model%t_rivvar%harea, &
        x => model%t_rivvar%x, &
        y => model%t_rivvar%y, &
        w => model%t_rivvar%w, &
        xo => model%t_rivvar%xo, &
        yo => model%t_rivvar%yo, &
        qs => model%t_rivvar%qs, &
        qn => model%t_rivvar%qn, &
        con => model%t_rivvar%con, &
        r => model%t_rivvar%r, &
        rn => model%t_rivvar%rn, &
        g => model%t_rivvar%g, &
        rho => model%t_rivvar%rho, &
        vkc => model%t_rivvar%vkc, &
        vbcds => model%t_calccond%vbcds, &
        transeqtype => model%t_calccond%transeqtype, &
        hmin => model%t_calccond%hmin, &
        vbc => model%t_calccond%vbc, &
        i_re_flag_i => model%t_calccond%i_re_flag_i, &
        i_re_flag_o => model%t_calccond%i_re_flag_o, &
        n_rest => model%t_calccond%n_rest, &
        i_tmp_count => model%t_calccond%i_tmp_count, &
        tmp_file_o => model%t_calccond%tmp_file_o, &
        tmp_caption => model%t_calccond%tmp_caption, &
        tmp_file_i => model%t_calccond%tmp_file_i, &
        tmp_dummy => model%t_calccond%tmp_dummy, &
        tmp_pass => model%t_calccond%tmp_pass, &
        opt_tmp => model%t_calccond%opt_tmp, &


        urelax => model%t_calccond%urelax, &
        erelax => model%t_calccond%erelax, &
        arelax => model%t_calccond%arelax, &
        levtype => model%t_calccond%levtype, &
        levchangeiter => model%t_calccond%levchangeiter, &
        levbegiter => model%t_calccond%levbegiter, &
        levenditer => model%t_calccond%levenditer, &
        evc => model%t_calccond%evc, &
        startlev => model%t_calccond%startlev, &
        endlev => model%t_calccond%endlev, &

        imod => model%t_calccond%imod, &
        drytype => model%t_calccond%drytype, &
        hiterinterval => model%t_calccond%hiterinterval, &
        hiterstop => model%t_calccond%hiterstop, &
        qprms => model%t_rivwvar%qprms, &
        isw => model%t_rivwvar%isw, &
        dea => model%t_rivwvar%dea, &
        eav => model%t_rivwvar%eav, &
        am => model%t_rivwvar%am, &
        bm => model%t_rivwvar%bm, &
        ccm => model%t_rivwvar%ccm, &
        dm => model%t_rivwvar%dm, &
        em => model%t_rivwvar%em, &
        qp => model%t_rivwvar%qp, &

        dqsde => model%t_rivwvar%dqsde, &
        dqnde => model%t_rivwvar%dqnde, &
        qv => model%t_rivwvar%qv, &
        qu => model%t_rivwvar%qu, &
        dq => model%t_rivwvar%dq, &
        de => model%t_rivwvar%de, &
        dnu => model%t_rivwvar%dnu, &

        dnv => model%t_rivwvar%dnv, &
        dsu => model%t_rivwvar%dsu, &
        dne => model%t_rivwvar%dne, &
        dnq => model%t_rivwvar%dnq, &
        dsv => model%t_rivwvar%dsv, &
        dse => model%t_rivwvar%dse, &
        dsq => model%t_rivwvar%dsq, &
        dqe => model%t_rivwvar%dqe, &

        dude => model%t_rivwvar%dude, &
        dvde => model%t_rivwvar%dvde, &
        eka => model%t_rivwvar%eka, &
        vout => model%t_rivwvar%vout, &
        uout => model%t_rivwvar%uout, &
        taus => model%t_rivvar%taus, &
        taun => model%t_rivvar%taun, &

        totcd => model%t_rivvar%totcd, &
        cd => model%t_rivvar%cd, &
        cdv => model%t_rivvar%cdv, &
        ns => model%t_rivvar%ns, &
        nn => model%t_rivvar%nn, &
        errorcode => model%t_rivvar%errorcode, &
        ptime => model%t_rivvar%ptime, &
        nct => model%t_rivvar%nct, &
        vardt => model%t_rivvar%vardt, &
        calccsed => model%t_calccond%calccsed, &
        iterplout => model%t_calccond%iterplout, &
        calcquasi3d => model%t_calccond%calcquasi3d, &
        calcquasi3drs => model%t_calccond%calcquasi3drs, &
        calcsedauto => model%t_calccond%calcsedauto, &
        !elevoffset => model%t_rivvar%elevoffset, &
        olddisch => model%t_rivvar%olddisch, &
        oldstage => model%t_rivvar%oldstage, &
        newstage => model%t_rivvar%newstage, &
        newdisch => model%t_rivvar%newdisch, &
        ddisch => model%t_rivvar%ddisch, &
        dstage => model%t_rivvar%dstage, &
        dsstage => model%t_rivvar%dsstage, &
        varstagetype => model%t_rivvartime%varstagetype, &
        vardischtype => model%t_rivvartime%vardischtype, &
        roughnesstype => model%t_calccond%roughnesstype, &
        e => model%t_rivvar%e, &
        hl => model%t_rivvar%hl, &
        u => model%t_rivvar%u, &
        v => model%t_rivvar%v, &
        eta => model%t_rivvar%eta, &
        ibc => model%t_rivvar%ibc, &
        icon => model%t_rivvar%icon, &
        nclusters => model%t_rivvar%nclusters, &
        nwetnodes => model%t_rivvar%nwetnodes, &
        noldwetnodes => model%t_rivvar%noldwetnodes, &
        interitm => model%t_calccond%interitm, &
        maxinteritermult => model%t_calccond%maxinteritermult, &
        q => model%t_calccond%q)

        nct = nct+1

        if(nct > 0.and.soltype == 0) then
            return
        endif
        if(nct > 0.and.model%t_end.eq.rvto%vardischstarttime) then
            return
        endif

        if(nct == 0) then
            vardt = model%t_calccond%fmdt
        else
            if(vardt > model%t_calccond%fmdt) then
                vardt = model%t_calccond%fmdt
            endif
        endif

        if(calccsed.eqv..false.) then
            vardt = model%t_calccond%fmdt
        endif

        olddisch = q
        oldstage = newstage

        if(nct > 0) then
            model%t = model%t+vardt
            if(model%t > model%t_end) then
                write(*,*) 'model time %f is greater than end time %f', model%t, model%t_end
                return
            endif
        endif

        if(vardischtype == 1) then !update discharge
            call getinterptimeseriesvalue(rvto, 1, model%t, newdisch)
            q = newdisch*1e6
        endif
        ddisch = q-olddisch

        if(varstagetype == 1) then !stage time series
            call getinterpratingcurvevalue(rvto, 1, model%t, newstage)
            !            newstage = newstage * 100.
        else if (varstagetype == 2) then !stage rating curve
            call getinterpratingcurvevalue(rvto, 1, newdisch, newstage)
            !            newstage = newstage * 100.
        endif

        dstage = newstage-oldstage

        !        if(cdtype == 1) then !variable roughness
        !            call setroughness(ns, nn, roughnesstype, newdisch, cd, znaught)
        !        endif

        if(vardischtype == 1) then
            !!  adjust water surface elevation on all wet nodes by dstage !!
            !!  check for minimum depth and turn nodes on or off.         !!
            call vardischeinit(rvo, cco, e, hl, u, v, eta, ibc, dstage, dsstage)
        endif

        !! turn nodes on - check re-wetting !!
        call updatewetting(rvo, cco, e,hl,u,v,eta, ibc, dsstage)

        !! set upstream boundary velocity  !!
        if(model%t_rivvartime%vardischtype == 1) then
            call vardischuinit(rvo, cco, u, v, w, eta, e, ibc, q, dstage)
        endif

        !! finds isolated pockets of wet nodes and deletes them
        call label_clusters(rvo, -1, ibc, nclusters, icon)
        call delete_clusters(rvo, ibc, icon, nclusters, nwetnodes)

        !!  check for ibc == 6 or ibc == 4  !!
        call updateibc(rvo)

        !!  if grid extension then reset extension topography !!
        !        if(nct > 1) then
        !            call resetgridextension()
        !        endif
        if(roughnesstype == 1) then !roughness as z0
            !		call zotocdtwo() !calculate cd from log profile
            call z0tocd(rvo, rvvo, cco) !calculate cd from 2-part profile
        endif


        do i = 1,10
            write(*,*)
        enddo

        write(*,*)
        write(*,*) 'variabledt', tmpvardt, 'dt', model%t_calccond%fmdt
        write(*,*) 'time step', nct, 'time', model%t, 'printtime', ptime
        write(*,*)
        if(calcquasi3d) then
            write(*,*) 'quasi-3d calculation is on'
            if(calcquasi3drs) then
                write(*,*) 'streamline curvature is on'
            endif
        else
            write(*,*) 'quasi-3d calculation is off'
        endif
        write(*,*)
        if(calccsed) then
            write(*,*) 'sediment-transport is on'
            if(calcsedauto) then
                write(*,*) 'time-stepping method is automatic'
            else
                write(*,*) 'time-stepping method is fixed'
            endif
            write(*,*) 'old wet nodes: ', noldwetnodes, ' new wet nodes: ', nwetnodes

        else
            write(*,*) 'sediment-transport is off'
        endif

        write(*,*)
        write(*,*) 'discharge', q/1e6, 'change discharge', ddisch/1e6
        write(*,*) 'stage', newstage/100., 'change stage', dstage/100.
        write(*,*)

        if(nct.ge.1) then
            model%itermax=interitm
            nodechange = abs((real(noldwetnodes)-real(nwetnodes)))
            model%itermax = real(nodechange*interitm)
            if(model%itermax < interitm)then
                model%itermax = interitm
            else if (model%itermax >= 5*interitm) then
                model%itermax = maxinteritermult*interitm
            endif
            noldwetnodes = nwetnodes
        endif

        if(model%itermax /= 0) then
            iter_loop: do iter = 1,model%itermax
                !       check if user cancelled simulation and if so exit
                call iric_check_cancel_f(errorcode)
                if(errorcode.eq.1) then
                    call dealloc_all(rvo, rvwo, cco, rvto, csedo)
                    return
                endif

                if(iwetdry.and.mod(iter,imod).eq.0.and.model%t_rivvartime%vardischtype == 0.and.iter.lt.hiterstop.and.nct.lt.1)then
                    call    updatewetting(rvo, cco, e, hl, u, v, eta, ibc, dsstage)
                endif

                call updateibc(rvo)

                do i = 1,ns
                    do j = 1,nn
                        if(ibc(i,j).eq.1) then
                            ustar=(totcd(i,j)*(u(i,j)**2.))**.5
                        else
                            ustar=(totcd(i,j)*(u(i,j)**2.+v(i,j)**2.))**.5
                        endif

                        !eka(i,j)=.005*200.*200.+vkc*hl(i,j)*ustar/6.
                        if(levtype == 0) then
                            eka(i,j) = evc+vkc*hl(i,j)*ustar/6.
                        else if(levtype.ne.0.and.nct.eq.0) then
                            if(iter < levbegiter) then
                                tmplev = startlev
                            elseif (iter >= levbegiter .and. iter < levenditer) then
                                tmplev = &
                                    (((startlev-endlev)/(levbegiter-levenditer))*&
                                    (iter-levbegiter))+startlev
                            elseif (iter > levenditer) then
                                tmplev = endlev
                            endif
                            eka(i,j) = tmplev+vkc*hl(i,j)*ustar/6.

                        else
                            eka(i,j) = endlev+vkc*hl(i,j)*ustar/6.
                        endif
                    enddo
                enddo


                do i=2,ns
                    do jswi=1,2
                        if(jswi.eq.1) then
                            j1=1
                            j2=nn
                            j3=1
                        else
                            j1=nn
                            j2=1
                            j3=-1
                        endif
                        do j=j1,j2,j3
                            if(ibc(i,j).eq.2) then
                                va=(v(i-1,j)+v(i,j))/4.
                            else if(j.lt.nn) then
                                va=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4.
                            endif
                            if(ibc(i,j).eq.4.or.ibc(i,j).eq.0) then
                                dude(i,j)=0.
                                u(i,j)=0.
                                cycle
                            endif

                            if(u(i,j).ge.0.or.i.eq.ns) then
                                ba=u(i,j)/dsu(i,j)
                                ta=ba*u(i-1,j)
                            else
                                ba=-1.*u(i,j)/dsu(i+1,j)
                                ta=ba*u(i+1,j)
                            endif

                            if((va.ge.0.and.ibc(i,j).ne.1).or.ibc(i,j).eq.2) then
                                bb=va/dnu(i,j)
                                tb=bb*u(i,j-1)
                            else
                                bb=-1.*va/dnu(i,j+1)
                                tb=bb*u(i,j+1)
                            endif

                            bc=-1.*va/(rn(i,j)*r(i))
                            if(ibc(i-1,j).eq.0) then
                                td = 0.
                            else
                                td=-1.*g*(e(i,j)-e(i-1,j))/dse(i,j)
                            endif

                            if(i.eq.ns) then
                                be=(2./(dsu(i,j)))*eka(i,j)/dsu(i,j)
                                te=be*u(i-1,j)
                            else
                                be=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)/dsu(i,j))+&
                                    (eka(i+1,j)/dsu(i+1,j)))
                                te=(4./(dsu(i,j)+dsu(i+1,j)))*((eka(i,j)*u(i-1,j)/dsu(i,j))+&
                                    (eka(i+1,j)*u(i+1,j)/dsu(i+1,j)))
                            endif

                            if(ibc(i,j).eq.2) then
                                bf=(2./dnu(i,j))*(totcd(i,j)*abs(u(i,j))+eka(i,j)/dnu(i,j))
                                tf=(2./dnu(i,j))*(eka(i,j)*u(i,j-1)/dnu(i,j))
                            elseif(ibc(i,j).eq.1) then
                                bf=(4./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+&
                                    totcd(i,j)*abs(u(i,j)))
                                tf=(4./(dnu(i,j)+dnu(i,j+1)))*(eka(i,j+1)*u(i,j+1)/dnu(i,j+1))
                            else
                                bf=-1.*eka(i,j)/(dnu(i,j)*rn(i,j)*r(i))+(2./(dnu(i,j)+&
                                    dnu(i,j+1)))*((eka(i,j+1)/dnu(i,j+1))+(eka(i,j)/dnu(i,j)))
                                tf=-1.*eka(i,j)*u(i,j-1)/(dnu(i,j-1)*rn(i,j-1)*r(i))+&
                                    (2./(dnu(i,j)+dnu(i,j+1)))*((eka(i,j+1)*u(i,j+1)/dnu(i,j+1))+&
                                    (eka(i,j)*u(i,j-1)/dnu(i,j)))
                            endif

                            bot2=rho*(hl(i,j)+hl(i-1,j))

                            if(bot2.eq.0.) then
                                bg=10.**12.
                            else
                                bg=2.*(totcd(i,j)*(u(i,j)**2.+va**2.)**.5)/bot2
                            endif

                            bot=ba+bb+bc+be+bf+bg

                            if(abs(bot).gt..000001) then
                                u(i,j)=u(i,j)*(1.-urelax)+urelax*(ta+tb+td+te+tf)/bot
                                if(i.eq.ns.and.u(i,j).lt.0) then
                                    if(vbcds .ne. 0) then
                                        u(i,j) = 0.
                                    endif
                                endif
                                dude(i,j)=-1.*g/(dse(i,j)*bot)
                            endif
                        enddo
                    enddo
                enddo

                do i=2,ns
                    do jswi=1,2
                        if(jswi.eq.1) then
                            j1=1
                            j2=nn
                            j3=1
                        else
                            j1=nn
                            j2=1
                            j3=-1
                        endif
                        do j=j1,j2,j3
                            if(ibc(i,j).eq.0.or.ibc(i,j).eq.1) then
                                v(i,j)=0.
                                dvde(i,j)=0.
                                cycle ! will need to make sure this change is correct
                            endif
                            if(i.eq.ns) then
                                ua=(u(i,j)+u(i,j-1))/2.
                            else
                                ua=(u(i,j)+u(i,j-1)+u(i+1,j)+u(i+1,j-1))/4.
                            endif
                            if(ua.ge.0.or.i.eq.ns) then
                                ba=ua/dsv(i,j)
                                ta=ba*v(i-1,j)
                            else
                                ba=-1.*ua/dsv(i+1,j)
                                ta=ba*v(i+1,j)
                            endif
                            if(v(i,j).ge.0.or.ibc(i,j).eq.2) then
                                bb=v(i,j)/dnv(i,j)
                                tb=bb*v(i,j-1)
                            else
                                bb=-1.*v(i,j)/dnv(i,j+1)
                                tb=bb*v(i,j+1)
                            endif

                            tc=-1.*(ua**2.)/(rn(i,j)*r(i))
                            if(ibc(i,j-1).eq.0) then
                                td = 0.
                            else
                                td=-1.*g*(e(i,j)-e(i,j-1))/dne(i,j)
                            endif
                            if(i.eq.ns) then
                                be=(1./dsv(i,j))*eka(i,j)/dsv(i,j)
                                te=be*v(i-1,j)
                            else
                                be1=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i+1,j)/dsv(i+1,j)
                                be2=(2./(dsv(i,j)+dsv(i+1,j)))*eka(i,j)/dsv(i,j)
                                be=be1+be2
                                te=be1*v(i+1,j)+be2*v(i-1,j)
                            endif
                            if(ibc(i,j).eq.2) then
                                bf1=(1./dnv(i,j))*eka(i,j)/dnv(i,j)
                                bf2=bf1
                                bf=bf1+bf2
                                tf=bf1*v(i,j-1)
                            else
                                bf1=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j+1)/dnv(i,j+1)
                                bf2=(2./(dnv(i,j)+dnv(i,j+1)))*eka(i,j)/dnv(i,j)
                                bf=bf1+bf2
                                tf=bf1*v(i,j+1)+bf2*v(i,j-1)
                            endif
                            bot2=rho*(hl(i,j)+hl(i,j-1))
                            if(bot2.eq.0.) then
                                bg=10.**12.
                            else
                                bg=2.*(totcd(i,j)*(ua**2.+v(i,j)**2.)**.5)/bot2
                            endif
                            bot=ba+bb+be+bf+bg
                            if(abs(bot).gt..000001) then
                                if(vbcds .eq. 0) then
                                    v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
                                else
                                    if(i.ne.ns) then
                                        v(i,j)=v(i,j)*(1.-urelax)+urelax*(ta+tb+tc+td+te+tf)/bot
                                    endif
                                endif
                                dvde(i,j)=-1.*g/(dn*bot)
                            endif
                        enddo
                    enddo
                enddo

                do i=1,ns
                    do j=1,nn
                        if(i.eq.1) then
                            area=dnq(i,j)*hl(i,j)
                        else
                            area=dnq(i,j)*(hl(i,j)+hl(i-1,j))/2.
                        endif
                        if(j.eq.1.or.j.eq.nn) area=area/2.
                        qu(i,j)=u(i,j)*area
                        dqsde(i,j)=dude(i,j)*area
                        if(j.eq.1) then
                            area=dsq(i,j)*hl(i,j)
                        else
                            !        area=.5*(dsq(i,j-1)+dsq(i,j))*(hl(i,j)+hl(i,j-1))/2.
                            area=dsq(i,j)*(hl(i,j)+hl(i,j-1))/2.

                        endif
                        qv(i,j)=v(i,j)*area
                        dqnde(i,j)=dvde(i,j)*area
                    enddo
                enddo

                do i=1,ns-1
                    do j=1,nn
                        if(j.eq.nn) then
                            dq(i,j)=qu(i+1,j)-qu(i,j)-qv(i,j)
                        else
                            dq(i,j)=qu(i+1,j)-qu(i,j)+qv(i,j+1)-qv(i,j)
                        endif
                    enddo
                enddo

                do i=1,ns-1
                    do j=1,nn
                        if(j.eq.nn) then
                            dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)
                        else
                            dqe(i,j)=dqsde(i,j)+dqsde(i+1,j)+dqnde(i,j)+dqnde(i,j+1)
                        endif
                        if(abs(dqe(i,j)).lt..000001) then
                            de(i,j)=0.
                        else
                            de(i,j)=dq(i,j)/dqe(i,j)
                        endif
                    enddo
                enddo

                do i=ns-1,2,-1
                    n=0.
                    do j=1,nn
                        !c       if(abs(dqe(i,j)).gt..000001) then
                        if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
                            n=n+1
                            if(j.ne.1) then
                                am(n)=-1.*dqnde(i,j)
                            endif
                            bm(n)=dqe(i,j)

                            if(j.eq.nn) then
                                ccm(n)=0.
                            else
                                ccm(n)=-1.*dqnde(i,j+1)
                            endif
                            dm(n)=dq(i,j)+dqsde(i+1,j)*de(i+1,j)+dqsde(i,j)*de(i-1,j)
                            isw(j)=1
                        else
                            isw(j)=0
                        endif
                    enddo
                    if(n.ge.1) then
                        call tridag(am,bm,ccm,dm,em,n,errorcode)
                        if(errorcode.eq.-1) then
                            !call write_error_cgns(str_in)
                            call dealloc_common2d(rvo)
                            if(vbc == 1) then
                                call dealloc_velbc(cco)
                            endif

                            if(calcquasi3d) then
                                call dealloc_common3d(rvo)
                                if(transeqtype == 2) then
                                    !call dealloc_csed3d_dt()
                                endif
                            endif
                            call dealloc_working(rvwo)
                            write(6,*) 'error at dqs', 'i = ', i

                            !pause
                            return
                        endif
                    endif
                    n2=0
                    do j=1,nn
                        if(isw(j).eq.1) then
                            n2=n2+1
                            de(i,j)=em(n2)
                        else
                            de(i,j)=0.
                        endif
                    enddo
                enddo

                do j=1,nn
                    n=0
                    do i=ns-1,2,-1
                        if(abs(dqe(i,j)).gt..000001.and.hl(i,j).gt.hmin) then
                            n=n+1
                            am(n)=-1.*dqsde(i,j)
                            bm(n)=dqe(i,j)
                            ccm(n)=-1.*dqsde(i+1,j)
                            if(j.eq.1) then
                                dm(n)=dq(i,j)+dqnde(i,j+1)*de(i,j+1)
                            elseif(j.eq.nn) then
                                dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)
                            else
                                dm(n)=dq(i,j)+dqnde(i,j)*de(i,j-1)+dqnde(i,j+1)*de(i,j+1)
                            endif
                            isw(i)=1
                        else
                            isw(i)=0
                        endif
                    enddo
                    if(n.ge.1) then
                        call tridag(am,bm,ccm,dm,em,n,errorcode)
                        if(errorcode.eq.-1) then
                            !call write_error_cgns(str_in)
                            call dealloc_common2d(rvo)
                            if(vbc == 1) then
                                call dealloc_velbc(cco)
                            endif

                            if(calcquasi3d) then
                                call dealloc_common3d(rvo)
                                if(transeqtype == 2) then
                                    !call dealloc_csed3d_dt()
                                endif
                            endif
                            call dealloc_working(rvwo)
                            write(6,*) 'error at dqn'

                            !pause
                            return
                        endif
                    endif
                    n2=0
                    do i=ns-1,2,-1
                        if(isw(i).eq.1) then
                            n2=n2+1
                            de(i,j)=em(n2)
                        else
                            de(i,j)=0.
                        endif
                    enddo
                enddo

                do i=1,ns
                    qp(i)=0.
                    eav(i)=0.
                    ibccount = 0
                    do j=1,nn
                        !        qp(i)=qp(i)+qu(i,j)
                        !        e(i,j)=e(i,j)+erelax*de(i,j)
                        if(ibc(i,j).ne.0) then
                            e(i,j)=e(i,j)+erelax*de(i,j)
                            ibccount = ibccount+1
                            eav(i)=eav(i)+e(i,j)
                        endif
                        qp(i)=qp(i)+qu(i,j)
                        hl(i,j)=e(i,j)-eta(i,j)
                        if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                            if(drytype.eq.0) then
                                if(ibc(i,j).ne.0) then
                                    u(i,j) = 0
                                    v(i,j) = 0
                                endif
                                ibc(i,j)=0
                            endif
                            hl(i,j)=hmin
                        endif
                        !unsure of the effect of the code block below rmcd 12/23/05
                        !change rmcd made sure i>1 4/30/10
                        if(i.gt.1) then
                            if(abs(u(i,j)).lt..00001) then
                                u(i,j)=0.
                            endif
                            if(abs(v(i,j)).lt..00001) then
                                v(i,j)=0.
                            endif
                        endif
                        !!!!!!!!!!!!
                    enddo
                    eav(i)=eav(i)/ibccount
                    do j=1,nn
                        if(ibc(i,j).eq.0)then
                            e(i,j) = eav(i)
                        endif
                    enddo
                enddo

                do  j=1,nn
                    u(1,j)=u(1,j)*q/qp(1)
                enddo

                dea(ns)=0.
                do  i=ns-1,1,-1
                    dinc=arelax*erelax*abs(eav(i)-eav(i+1))*(1.-qp(i+1)/q)
                    dea(i)=dea(i+1)+dinc
                    do  j=1,nn
                        !! code change by rmcd !!
                        if(ibc(i,j).ne.0)then
600                         e(i,j)=e(i,j)+dea(i)
                            !                else
                            !                    e(i,j) = eav(i)+dea(i) !! added so that e(i,j) in dry nodes follows
                            !                                           !! average eav(i) during time-dependent runs
                        endif
                        !! end code change by rmcd !!
                    enddo
                enddo

                sumqpdiff = 0.
                do i=1,ns
                    sumqpdiff = sumqpdiff + (((qp(i)-q)/q*100.))**2.
                enddo

                qprms(iter) = sqrt(sumqpdiff/ns)
                !write(*,*) 'after qprms'
                if(isnan(qprms(iter))) then
                    errorcode=-1
                    !call write_error_cgns(str_in)
                    !			call dealloc_all()
                    call dealloc_common2d(rvo)
                    if(vbc == 1) then
                        call dealloc_velbc(cco)
                    endif

                    if(calcquasi3d) then
                        call dealloc_common3d(rvo)
                        if(transeqtype == 2) then
                            !call dealloc_csed3d_dt()
                        endif
                    endif
                    call dealloc_working(rvwo)
                    write(6,*) 'error at qprms'
                    !pause
                    return
                endif

                if(nct == dbgtimestep) then
                    if(debugstop.eq.1.and.iter == dbgiternum) then
                        solindex = solindex+1
                        call cg_iric_write_sol_time_f(model%t, ier)

                        call write_cgns2(rvo, cco, rvwo, model%t, q)
                        !			        call write_timestep_cgns(str_in, nct, tottime)
                        !			        call write_timeiter_cgns(str_in)
                        call dealloc_all(rvo, rvwo, cco, rvto, csedo)
                        !pause

                        return
                    endif
                endif
                if(levtype == 0) then
                    write(*,*) 'iteration: ', iter,' mean error on discharge:', qprms(iter)
                else
                    write(*,*)'iteration:', iter, ' lev:', tmplev,' mean error on discharge: ', qprms(iter)
                endif

                if(iterationout.and.mod(iter,iterplout).eq.0) then

                    do i = 1,ns
                        do j = 1,nn
                            if(ibc(i,j).eq.0) then
                                v(i,j) = 0.
                                u(i,j) = 0.
                            endif
                        enddo
                    enddo

                    do  i=1,ns
                        do  j=1,nn
                            !                if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                            if(ibc(i,j).ne.-1.or.hl(i,j).le.hmin) then
                                u(i,j)=0.
                                vout(i,j)=0.
                                !    go to 888
                                !endif
                            elseif(i.eq.1) then
                                vout(i,j)=v(i,j)
                                !    go to 888
                                !endif

                            elseif(j.eq.nn) then
                                vout(i,j)=(v(i-1,j)+v(i,j))/4.
                            elseif(j.eq.1) then
                                vout(i,j) = 0.
                            elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
                                vout(i,j)=v(i,j)
                            elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
                                vout(i,j)=v(i,j+1)
                            else
                                vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
                            endif

888                         if(i.eq.1) then
                                dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
                            else
                                dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
                            endif
                            taus(i,j)=dube*u(i,j)
                            taun(i,j)=dube*vout(i,j)
                        enddo
                    enddo

                    solindex = solindex+1
                    call cg_iric_write_sol_time_f(model%t, ier)
                    call write_cgns2(rvo, cco, rvwo, model%t, q)

                endif

            end do iter_loop

            if(i_re_flag_o.eq.1.and.nct.eq.0) then
                if(i_tmp_count <= n_rest) then
                    if(iter.eq.opt_tmp(i_tmp_count)) then
                        tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
                        open(502,file=tmp_file_o(i_tmp_count) &
                            ,status='unknown',form='unformatted')
                        !
                        write(502) model%t,solindex,model%t_calccond%fmdt
                        write(502) ns,nn
                        !
                        write(502) ((eta(i,j),i=1,ns),j=1,nn)
                        write(502) ((u(i,j),i=1,ns),j=1,nn)
                        write(502) ((v(i,j),i=1,ns),j=1,nn)
                        write(502) ((ibc(i,j),i=1,ns),j=1,nn)
                        write(502) ((e(i,j),i=1,ns),j=1,nn)
                        write(502) ((hl(i,j),i=1,ns),j=1,nn)
                        write(502) (hav(i),i=1,ns)

                        close(502)

                        i_tmp_count = i_tmp_count +1
                    endif
                endif
            endif

64          format(25f7.3)
            !66      format(a15,i6,a25,f15.5)
66          format(a, a, i5,a, f15.5)
67          format(a, a, i5,a, f15.5,a, f15.5)
68          format(a)
75          format(25f7.0)
85          format(25f6.4)
86          format(25f8.3)
84          format(3i5,f6.2,2f12.1,2g12.4)

            do i = 1,ns
                do j = 1,nn
                    if(ibc(i,j).eq.0) then
                        v(i,j) = 0.
                        u(i,j) = 0.
                    endif
                enddo
            enddo
        else

            !601 continue
            do i=1,ns
                do j=1,nn
                    !          if(ibc(i,j).eq.0.or.hl(i,j).le.hmin) then
                    if(ibc(i,j).ne.-1) then
                        u(i,j)=0.
                        vout(i,j)=0.
                        go to 88
                    endif
                    if(i.eq.1) then
                        vout(i,j)=v(i,j)
                        go to 88
                    endif

                    if(j.eq.nn) then
                        vout(i,j)=(v(i-1,j)+v(i,j))/4.
                    elseif(j.eq.1) then
                        vout(i,j) = 0.
                    elseif(ibc(i,j).ne.0.and.ibc(i,j+1).eq.0) then
                        vout(i,j)=v(i,j)
                    elseif(ibc(i,j).ne.0.and.ibc(i,j-1).eq.0) then
                        vout(i,j)=v(i,j+1)
                    else
                        vout(i,j)=(v(i-1,j)+v(i,j)+v(i,j+1)+v(i-1,j+1))/4. !jmn
                    endif
88                  if(i.eq.1) then
                        dube=(cd(i,j))*sqrt(u(i,j)**2+vout(i,j)**2)
                    else
                        dube=(cd(i,j)+cd(i-1,j))*.5*sqrt(u(i,j)**2+vout(i,j)**2)
                    endif
                    taus(i,j)=dube*u(i,j)
                    taun(i,j)=dube*vout(i,j)
                enddo
            enddo

            !c        call vert
            if(calcquasi3d) then
                call vert(rvo, rvvo, cco)

            endif
            if(calccsed) then
                if(transeqtype == 2) then
                    !call csed_dt(vardt,nct,nsteps, tmpvardt)
                else
                    call csed2(csedo, rvo, cco, vardt,nct,nsteps, tmpvardt)
                endif
            else
                call stressdiv(rvo, cco, ibc, qs, qn, taus, taun, con, rn, r)
            endif
            if(calcsedauto) then
                vardt = tmpvardt
            endif

4000        format(6f10.2)
            if(nct == 0) then
                solindex = solindex+1
                call calc_area(ns, nn, phirotation, x, y, xo, yo, nm, dn, harea)
                call cg_iric_write_sol_time_f(model%t, ier)
                call write_cgns2(rvo, cco, rvwo, model%t, q)
                !                if(calcquasi3d.and.io_3doutput) then
                !                    call write_cgns3d_grid()
                !                endif
                if(calcquasi3d.and.io_3doutput) then
                    call write_cgns3d_grid(rvo)
                    !                    call write_cgns3d_solgrid()
                    call write_cgns3d_fixedbed(rvo, rvto, solindex, model%t, q)
                endif

            else
                if(model%t >= ptime)then
                    solindex = solindex+1
                    call calc_area(ns, nn, phirotation, x, y, xo, yo, nm, dn, harea)
                    call cg_iric_write_sol_time_f(model%t, ier)
                    call write_cgns2(rvo, cco, rvwo, model%t, q)
                    if(calcquasi3d.and.io_3doutput) then
                        !                        call write_cgns3d_solgrid()
                        call write_cgns3d_fixedbed(rvo, rvto, solindex, model%t, q)
                        !                    elseif(calcquasi3d.and.io_3doutput.and.calccsed) then
                        !!                        call write_cgns3d_moveablebed(tottime, q)
                    endif
                    ptime = ptime+(iplinc*model%t_calccond%fmdt)
                endif
                if(i_re_flag_o.eq.1.and.nct.ne.0) then
                    if(i_tmp_count <= n_rest) then
                        if(model%t.ge.opt_tmp(i_tmp_count)) then
                            tmp_file_o(i_tmp_count)=trim(tmp_pass)//tmp_file_o(i_tmp_count)  !i110419
                            open(502,file=tmp_file_o(i_tmp_count) &
                                ,status='unknown',form='unformatted')
                            !
                            write(502) model%t,solindex,model%t_calccond%fmdt
                            write(502) ns,nn
                            !
                            write(502) ((eta(i,j),i=1,ns),j=1,nn)
                            write(502) ((u(i,j),i=1,ns),j=1,nn)
                            write(502) ((v(i,j),i=1,ns),j=1,nn)
                            write(502) ((ibc(i,j),i=1,ns),j=1,nn)
                            write(502) ((e(i,j),i=1,ns),j=1,nn)
                            write(502) ((hl(i,j),i=1,ns),j=1,nn)
                            write(502) (hav(i),i=1,ns)

                            close(502)

                            i_tmp_count = i_tmp_count +1
                        endif
                    endif
                endif
            endif
        endif
    end associate
    end subroutine solve_fm

    subroutine get_grid_2d_ssvec(model, ctype, tmp)
    implicit none
    type(fastmech_model), intent(in) :: model
    type(rivvar) :: rvo
    type(riv_w_var) ::rwvo
    character(len=*) :: ctype
    real,dimension(:), intent(inout) :: tmp(:) !changed this to just real but could create problems
    real(kind=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    integer :: i, j, count, ier
    real(kind=mp), allocatable, dimension(:,:) :: tmpval
    rvo =  model%t_rivvar
    rwvo = model%t_rivwvar
    allocate(tmpval(rvo%ns2,rvo%nn), stat=ier)
    do i=1,rvo%ns2
        do j= 1,rvo%nn
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            ux = rvo%taus(i,j)*rcos - rvo%taun(i,j)*rsin
            uy = rvo%taus(i,j)*rsin + rvo%taun(i,j)*rcos

            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            select case (ctype)
            case('ShearStressX')
                tmpval(i,j) = (xx + rvo%xshift)/100.0d0
            case('ShearStressY')
                tmpval(i,j) = (yy + rvo%yshift)/100.0d0
            case default
                tmpval(i,j) = -9999
                return
            end select
        end do
    end do
    tmp = reshape(tmpval, (/rvo%ns2*rvo%nn/))
    deallocate(tmpval, stat=ier)
    end subroutine get_grid_2d_ssvec

    subroutine get_grid_2d_velvec(model, ctype, tmp)
    implicit none
    type(fastmech_model), intent(in) :: model
    type(rivvar) :: rvo
    type(riv_w_var) ::rwvo
    character(len=*) :: ctype
    real,dimension(:), intent(inout) :: tmp(:) !changed this to just real but could create problems
    real(kind=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    integer :: i, j, count, ier
    real(kind=mp), allocatable, dimension(:,:) :: tmpval
    rvo =  model%t_rivvar
    rwvo = model%t_rivwvar
    allocate(tmpval(rvo%ns2,rvo%nn), stat=ier)
    do i=1,rvo%ns2
        do j= 1,rvo%nn
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            ux = rvo%u(i,j)*rcos - rwvo%vout(i,j)*rsin
            uy = rvo%u(i,j)*rsin + rwvo%vout(i,j)*rcos

            xx = ux*rvo%fcos - uy*rvo%fsin
            yy = ux*rvo%fsin + uy*rvo%fcos
            select case (ctype)
            case('VelocityX')
                tmpval(i,j) = (xx + rvo%xshift)/100.0d0
            case('VelocityY')
                tmpval(i,j) = (yy + rvo%yshift)/100.0d0
            case default
                tmpval(i,j) = -9999
                return
            end select
        end do
    end do
    tmp = reshape(tmpval, (/rvo%ns2*rvo%nn/))
    deallocate(tmpval, stat=ier)
    end subroutine get_grid_2d_velvec

    subroutine get_grid_2d_coord(model, ctype, tmp)
    implicit none
    type(fastmech_model), intent(in) :: model
    type(rivvar) :: rvo
    character(len=*) :: ctype
    real,dimension(:), intent(inout) :: tmp(:) !changed this to just real but could create problems
    real(kind=mp) :: xx, yy, ux, uy
    double precision :: rcos, rsin
    integer :: i, j, count, ier
    real(kind=mp), allocatable, dimension(:,:) :: tmpval
    rvo =  model%t_rivvar
    allocate(tmpval(rvo%ns2,rvo%nn), stat=ier)
    do i=1,rvo%ns2
        do j= 1,rvo%nn
            rcos = cos(rvo%phirotation(i))
            rsin = sin(rvo%phirotation(i))
            ux = rvo%x(i,j)*rcos - rvo%y(i,j)*rsin
            uy = rvo%x(i,j)*rsin + rvo%y(i,j)*rcos

            xx = rvo%x(i,j)*rvo%fcos - rvo%y(i,j)*rvo%fsin
            yy = rvo%x(i,j)*rvo%fsin + rvo%y(i,j)*rvo%fcos
            select case (ctype)
            case('x')
                tmpval(i,j) = (xx + rvo%xshift)/100.0d0
            case('y')
                tmpval(i,j) = (yy + rvo%yshift)/100.0d0
                case default
                tmpval(i,j) = -9999
                return
            end select
        end do
    end do
    tmp = reshape(tmpval, (/rvo%ns2*rvo%nn/))
    deallocate(tmpval, stat=ier)
    end subroutine get_grid_2d_coord

    subroutine get_grid_3d_coord(model, ctype, tmp)
    implicit none
    type(fastmech_model), intent(in) :: model
    type(rivvar) :: rvo
    character(len=*) :: ctype
    real, dimension(:), intent(inout) :: tmp(:)
    integer :: i, j, k, ier
    double precision :: xx, yy, ux, uy
    double precision :: rcos, rsin
    real(kind=mp), allocatable, dimension(:,:,:) :: tmpval
    rvo = model%t_rivvar
    allocate(tmpval(rvo%ns2,rvo%nn,rvo%nz), stat=ier)
    do k = 1,rvo%nz
        do j=1,rvo%nn
            do i= 1,rvo%ns2
                rcos = cos(rvo%phirotation(i))
                rsin = sin(rvo%phirotation(i))

                ux = rvo%x(i,j)*rcos - rvo%y(i,j)*rsin
                uy = rvo%x(i,j)*rsin + rvo%y(i,j)*rcos

                xx = rvo%x(i,j)*rvo%fcos - rvo%y(i,j)*rvo%fsin
                yy = rvo%x(i,j)*rvo%fsin + rvo%y(i,j)*rvo%fcos
                select case (ctype)
                case('x')
                    tmpval(i,j,k) = (xx + rvo%xshift)/100.
                case('y')
                    tmpval(i,j,k) = (yy + rvo%yshift)/100.
                case('z')
                    tmpval(i,j,k) = (rvo%zz(i,j,k)/100.) !+ rvo%elevoffset
                end select
                !ty(i,j,k) = (yy + yshift)/100.
                !tz(i,j,k) = (zz(i,j,k)/100.) + elevoffset
            end do
        end do
    end do
    tmp = reshape(tmpval, (/rvo%ns2*rvo%nn*rvo%nz/))
    deallocate(tmpval, stat=ier)
    end subroutine get_grid_3d_coord

    subroutine print_info(model)
    type (fastmech_model), intent (in) :: model

    write(*,"(a10, i8)") "n_x:", model%n_x
    write(*,"(a10, i8)") "n_y:", model%n_y
    !write(*,"(a10, f8.2)") "dx:", model%dx
    !write(*,"(a10, f8.2)") "dy:", model%dy
    !write(*,"(a10, f8.2)") "alpha:", model%alpha
    write(*,"(a10, f8.2)") "dt:", model%dt
    write(*,"(a10, f8.2)") "t:", model%t
    write(*,"(a10, f8.2)") "t_end:", model%t_end
  end subroutine print_info
    end module fastmech
