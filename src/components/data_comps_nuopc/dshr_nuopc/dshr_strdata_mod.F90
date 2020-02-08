module dshr_strdata_mod

  ! -----------------------------------------------------
  ! holds data and methods to advance data models
  ! -----------------------------------------------------
  
  ! TODO: get the logic correct for doing the vector mapping
  ! TODO: make shr_strdata a linked list
  ! TODO: how to handle  z dimension?
  ! TODO: add scmlon and scmlat functionality
  ! TODO: implement fullfill stream functionality
  ! TODO: implement the interface that is fortran callable
  ! TODO: add functionality if the stream mesh needs to be created from a grid

  use ESMF
  use NUOPC
  use shr_kind_mod    , only : in=>shr_kind_in, r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : shr_const_pi
  use shr_mpi_mod     , only : shr_mpi_bcast
  use shr_cal_mod     , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod     , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod     , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod     , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_nl_mod      , only : shr_nl_find_group_name
  use shr_mpi_mod     , only : shr_mpi_bcast
  use shr_const_mod   , only : shr_const_cDay
  use shr_string_mod  , only : shr_string_listGetNum, shr_string_listGetName
  use shr_string_mod  , only : shr_string_listIsValid
  use shr_tinterp_mod , only : shr_tInterp_getAvgCosz
  use shr_tinterp_mod , only : shr_tInterp_getFactors
  use dshr_stream_mod , only : shr_stream_taxis_cycle
  use dshr_stream_mod , only : shr_stream_taxis_extend
  use dshr_stream_mod , only : shr_stream_taxis_limit
  use dshr_stream_mod , only : shr_stream_file_null
  use dshr_stream_mod , only : shr_stream_streamType
  use dshr_stream_mod , only : dshr_stream_findBounds
  use dshr_stream_mod , only : dshr_stream_getFilePath
  use dshr_stream_mod , only : dshr_stream_getPrevFileName
  use dshr_stream_mod , only : dshr_stream_getNextFileName
  use dshr_stream_mod , only : dshr_stream_getFileFieldName
  use dshr_stream_mod , only : dshr_stream_getModelFieldList
  use dshr_stream_mod , only : dshr_stream_getCurrFile
  use dshr_stream_mod , only : dshr_stream_setCurrFile
  use dshr_methods_mod, only : chkerr 
  use dshr_methods_mod, only : FB_getFieldN, FB_FldChk, FB_Diagnose
  use dshr_methods_mod, only : FB_Regrid, FB_getFldPtr 
  use perf_mod       
  use pio            

  implicit none
  private

  ! !PUBLIC TYPES:
  public shr_strdata_type

  ! !PUBLIC MEMBER FUNCTIONS:
  public ::  dshr_strdata_init
  public ::  dshr_strdata_advance
  public ::  dshr_strdata_setOrbs
  public ::  dshr_strdata_restRead
  public ::  dshr_strdata_restWrite

  private :: dshr_strdata_print
  private :: dshr_strdata_readLBUB
  private :: dshr_strdata_readstrm
  !private :: dshr_strdata_readstrm_fullfile 

  ! !PRIVATE:

  ! constants
  integer                              :: debug    = 0  ! local debug flag
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  real(R8)         ,parameter, private :: dtlimit_default = 1.5_R8

  ! mapping types
  integer , public, parameter :: nmappers       = 8
  integer , public, parameter :: mapunset       = 0
  integer , public, parameter :: mapbilnr       = 1
  integer , public, parameter :: mapconsf       = 2
  integer , public, parameter :: mapconsd       = 3
  integer , public, parameter :: mappatch       = 4
  integer , public, parameter :: mapfcopy       = 5
  integer , public, parameter :: mapnstod       = 6 ! nearest source to destination
  integer , public, parameter :: mapnstod_consd = 7 ! nearest source to destination followed by conservative dst
  integer , public, parameter :: mapnstod_consf = 8 ! nearest source to destination followed by conservative frac
  character(len=*) , public, parameter :: mapnames(nmappers) = &
       (/'bilinear   ','consf      ','consd      ',&
         'patch      ','fcopy      ','nstod      ',&
         'nstod_consd','nstod_consf'/)

  ! stream maximum sizes
  integer ,parameter          :: nStrMax = 30
  integer ,parameter          :: nVecMax = 30

  type shr_strdata_type
     character(CL)                  :: dataMode                                ! flags physics options wrt input data
     integer                        :: io_type                                 ! pio setting
     integer                        :: io_format                               ! pio setting
     type(iosystem_desc_t), pointer :: pio_subsystem => null()                 ! pio setting
     integer                        :: nvectors                                ! number of vectors
     type(ESMF_Mesh)                :: mesh_model                              ! model mesh
     integer                        :: ustrm (nVecMax) = 0                     ! vector interpolation stream  n -> model
     integer                        :: vstrm (nVecMax) = 0                     ! vector interpolation strearm n -> model
     integer                        :: ymd                                     ! model ymd
     integer                        :: tod                                     ! model tod
     character(CL)                  :: calendar                                ! model calendar for ymd,tod
     real(R8)                       :: eccen                                   ! at model time for cosz t-interp method
     real(R8)                       :: mvelpp                                  ! at model time for cosz t-interp method
     real(R8)                       :: lambm0                                  ! at model time for cosz t-interp method
     real(R8)                       :: obliqr                                  ! at model time for cosz t-interp method
     integer                        :: modeldt                                 ! model dt in seconds for cosz t-interp method
     character(CL)                  :: allocstring = 'strdata_allocated'
     !
     integer                        :: nstreams                                ! actual number of streams
     type(ESMF_Mesh)                :: mesh_streams(nstrMax)                   ! mesh for each stream
     character(CL)                  :: streamfiles(nStrMax)                    ! stream description file names
     character(CL)                  :: tintalgo(nStrMax)                       ! time interpolation algorithm
     character(CL)                  :: taxMode (nStrMax)                       ! time axis cycling mode
     real(R8)                       :: dtlimit (nStrMax)                       ! dt max/min limit
     character(CL)                  :: vector_names(nVecMax)                   ! define vectors to vector map
     character(CL)                  :: mapalgo (nStrMax)                       ! scalar map algorithm
     character(CL)                  :: readmode(nStrMax)                       ! file read mode
     type(shr_stream_streamType)    :: streams (nStrMax)                       ! stream info
     type(io_desc_t)                :: pio_iodesc(nStrMax)                     ! pio setting
     type(ESMF_FieldBundle)         :: FB_stream_lbound(nstrMax)               ! stream n field bundle for lb of time period (stream grid)
     type(ESMF_FieldBundle)         :: FB_stream_ubound(nstrMax)               ! stream n field bundle for ub of time period (stream grid)
     type(ESMF_FieldBundle)         :: FB_model_lbound(nstrMax)                ! stream n field bundle for lb of time period (model grid)
     type(ESMF_FieldBundle)         :: FB_model_ubound(nstrMax)                ! stream n field bundle for ub of time period (model grid)
     type(ESMF_FieldBundle)         :: FB_model(nStrMax)                       ! stream n field bundle for model time (model grid)
     type(ESMF_FieldBundle), pointer:: FB_stream_alltimes(:,:) => null()       ! field bundle for stream n for all time slices for stream
     type(ESMF_RouteHandle)         :: RH_stream2model(nMappers,nstrMax)       ! stream n -> model mesh mapping
     type(ESMF_Field)               :: field_coszen(nStrMax)                   ! needed for coszen time interp
     integer                        :: ymdLB(nStrMax)                          ! stream n ymd lower bound
     integer                        :: todLB(nStrMax)                          ! stream n time of day lower bound
     integer                        :: ymdUB(nStrMax)                          ! stream n ymd upper bound
     integer                        :: todUB(nStrMax)                          ! stream n time of day upper bound
     real(R8)                       :: dtmin(nStrMax)
     real(R8)                       :: dtmax(nStrMax)
  end type shr_strdata_type

  real(R8),parameter,private :: deg2rad = SHR_CONST_PI/180.0_R8

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  ! subroutine shr_strdata_create_from_inline(sdat, name, mpicom, compid, &
  !      yearFirst, yearLast, yearAlign, offset, &
  !      stream_meshfile, stream_filePath, stream_filenames, &
  !      stream_fieldnames, model_fieldnames, stream_readmode, &
  !      taxMode, dtlimit, tintalgo, mapalgo, calendar, rc)

  !      ! Set strdata and stream info from inline fortran interface.
  !      ! Note: When this is called, previous settings are reset to defaults
  !      ! and then the values passed are used.

  !   ! input/output arguments
  !   type(shr_strdata_type),intent(inout):: sdat                 ! strdata data data-type
  !   character(*)          ,intent(in)   :: name                 ! name of strdata
  !   integer               ,intent(in)   :: mpicom               ! mpi comm
  !   integer               ,intent(in)   :: compid               ! component id (needed for pio)
  !   integer               ,intent(in)   :: yearFirst            ! first year to use
  !   integer               ,intent(in)   :: yearLast             ! last  year to use
  !   integer               ,intent(in)   :: yearAlign            ! align yearFirst with this model year
  !   integer               ,intent(in)   :: offset               ! offset in seconds of stream data
  !   character(*)          ,intent(in)   :: stream_meshfile      ! stream mesh filename
  !   character(*)          ,intent(in)   :: stream_filepath      ! path to stream files
  !   character(*)          ,intent(in)   :: stream_filenames(:)  ! stream filenames in stream_filepath
  !   character(*)          ,intent(in)   :: stream_fieldnames(:) ! stream field names
  !   character(*)          ,intent(in)   :: model_fieldnames(:)  ! model field names corresponding to stream_fieldnames
  !   character(*),optional ,intent(in)   :: stream_readmode      ! file read mode (time slice or all times
  !   character(*),optional ,intent(in)   :: tintalgo             ! time interpolation algorithm
  !   character(*),optional ,intent(in)   :: mapalgo              ! scalar map algorithm
  !   character(*),optional ,intent(in)   :: taxMode
  !   real(R8)    ,optional ,intent(in)   :: dtlimit
  !   character(*),optional, intent(in)   :: calendar             ! model calendar

  !   ! local variables
  !   character(*),parameter :: subName = "(shr_strdata_create) "
  !   character(*),parameter ::   F00 = "('(shr_strdata_create) ',8a)"
  !   !-------------------------------------------------------------------------------

  !   ! The following only sets defaults - but does not read the namelist
  !   ! sdat%nstreams = 1
  !   ! call shr_strdata_pioinit(sdat, compid)

  !   ! call shr_strdata_readnml(sdat)
  !   ! if (present(taxMode)) then
  !   !    sdat%taxMode(1) = taxMode
  !   !    if (trim(sdat%taxMode(1)) == trim(shr_stream_taxis_extend)) sdat%dtlimit(1) = 1.0e30
  !   ! endif
  !   ! if (present(dtlimit))  sdat%dtlimit(1) = dtlimit
  !   ! if (present(mapalgo))  sdat%mapalgo(1) = mapalgo
  !   ! if (present(tintalgo)) sdat%tintalgo(1) = tintalgo
  !   ! if (present(calendar)) sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))

  !   ! call shr_stream_set(sdat%streamfiles(1), yearFirst, yearLast, yearAlign, offset, taxMode, &
  !   !      stream_meshfile, stream_fieldnames, fmodel_fieldnames, stream_filepath, filename) !???

  !   ! ! TODO: need to fill in information about the model mesh here

  !   ! !call shr_strdata_init(sdat, mpicom, compid, gsmap=gsmap, ggrid=ggrid, nxg=nxg, nyg=nyg, nzg=1)

  ! end subroutine shr_strdata_create_from_inline

  !===============================================================================

  subroutine dshr_strdata_create_from_cap(mesh_model, clock, nml_filename, compid, &
       mpicom, my_task, master_task, logunit, sdat, rc)

    ! This is ONLY relevant to data models caps - and not to calls from component models that
    ! leverage stream data functionality

    ! input/output arguments
    type(ESMF_Mesh)        ,intent(in)    :: mesh_model   ! data model mesh
    type(ESMF_Clock)       ,intent(in)    :: clock        ! data model clock
    character(*)           ,intent(in)    :: nml_filename ! shr_strdata namelist filename
    integer                ,intent(in)    :: compid       ! component id needed to initialize pio
    integer                ,intent(in)    :: mpicom       ! mpi communicator
    integer                ,intent(in)    :: my_task      ! my task number, 0 is default
    integer                ,intent(in)    :: master_task  ! master task number, 0 is default
    integer                ,intent(in)    :: logunit      ! output logunit
    type(shr_strdata_type) ,intent(inout) :: sdat         ! strdata data data-type
    integer, optional      ,intent(out)   :: rc           ! return code

    ! local variables
    type(ESMF_CalKind_Flag) :: esmf_caltype         ! esmf calendar type
    integer                 :: rCode                ! return code
    integer                 :: nUnit                ! fortran i/o unit number
    integer                 :: n                    ! generic loop index
    character(CL)           :: dataMode             ! flags physics options wrt input data
    character(CL)           :: streamfiles(nStrMax) ! stream description file names
    character(CL)           :: taxMode(nStrMax)     ! time axis cycling mode
    real(R8)                :: dtlimit(nStrMax)     ! delta time limiter
    character(CL)           :: vectors(nVecMax)     ! define vectors to vector map
    character(CL)           :: mapalgo(nStrMax)     ! scalar map algorithm
    character(CL)           :: tintalgo(nStrMax)    ! time interpolation algorithm
    character(CL)           :: readmode(nStrMax)    ! file read mode (single time slice or full)
    character(CL)           :: fileName             ! generic file name
    integer                 :: yearFirst            ! first year to use in data stream
    integer                 :: yearLast             ! last  year to use in data stream
    integer                 :: yearAlign            ! data year that aligns with yearFirst
    character(*),parameter  :: F00 = "('(shr_strdata_readnml) ',8a)"
    character(*),parameter  :: subName = "(shr_strdata_readnml) "
    !-------------------------------------------------------------------------------

    ! define namelist (NOTE: streams input has txtfile, yearFirst, yearLast, yearAlign as one entry)
    namelist / shr_strdata_nml / dataMode, streamfiles, taxMode, dtlimit, vectors, mapalgo, tintalgo, readmode

    rc = ESMF_SUCCESS

    ! initialize sdat model mesh
    sdat%mesh_model =  ESMF_MeshCreate(mesh_model, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initialize sdat calendar
    call ESMF_ClockGet(clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       sdat%calendar = trim(shr_cal_calendarName(trim(shr_cal_noleap)))
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       sdat%calendar = trim(shr_cal_calendarName(trim(shr_cal_gregorian)))
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(shr_cal_calendarName(trim(shr_cal_gregorian))))
    end if

    ! initialize shr_strdata pio
    call shr_strdata_pioinit(sdat, compid)

    ! read input namelist on master_task only
    if (my_task == master_task) then

       ! set default values for namelist vars
       dataMode        = 'NULL'
       streamfiles(:)  = trim(shr_strdata_nullstr)
       taxMode(:)      = trim(shr_stream_taxis_cycle)
       dtlimit(:)      = dtlimit_default
       vectors(:)      = trim(shr_strdata_nullstr)
       mapalgo(:)      = 'bilinear'
       tintalgo(:)     = 'linear'
       readmode(:)     = 'single'

       ! read input namelist
       write(logunit,F00) 'reading input namelist file: ',trim(nml_filename)
       open (newunit=nUnit,file=trim(nml_filename),status="old",action="read")
       call shr_nl_find_group_name(nUnit, 'shr_strdata_nml', status=rCode)
       if (rCode == 0) then
          read (nUnit, nml=shr_strdata_nml, iostat=rCode)
          if (rCode /= 0) then
             call shr_sys_abort(subName//": namelist read error "//trim(nml_filename))
          end if
       end if
       close(nUnit)

       ! copy temporary/local namelist vars into data structure
       sdat%nstreams = 0
       do n=1,nStrMax
          call dshr_stream_default(sdat%streamfiles(n))
       enddo
       sdat%dataMode        = dataMode
       sdat%streamfiles(:)  = streamfiles(:)
       sdat%taxMode(:)      = taxMode(:)
       sdat%dtlimit(:)      = dtlimit(:)
       sdat%vector_names(:) = vectors(:)
       sdat%mapalgo(:)      = mapalgo(:)
       sdat%tintalgo(:)     = tintalgo(:)
       sdat%readmode(:)     = readmode(:)
       do n=1,nStrMax
          if (trim(streamfiles(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nstreams = max(sdat%nstreams,n)
          end if
          if (trim(sdat%taxMode(n)) == trim(shr_stream_taxis_extend)) then
             sdat%dtlimit(n) = 1.0e30
          end if
       end do
       sdat%nvectors = 0
       do n=1,nVecMax
          if (trim(vectors(n)) /= trim(shr_strdata_nullstr)) sdat%nvectors = n
       end do

       ! initialize sdat%streams data types
       do n = 1,sdat%nstreams
          if (trim(sdat%streamfiles(n)) /= shr_strdata_nullstr) then
             ! extract fileName (stream description text file), yearAlign, yearFirst, yearLast from sdat%streamfiles(n)
             call dshr_stream_parseInput(sdat%streamfiles(n), fileName, yearAlign, yearFirst, yearLast)
             ! initialize stream datatype, read description text file
             call dshr_stream_init(sdat%streamfiles(n), fileName, yearFirst, yearLast, yearAlign, trim(sdat%taxMode(n)))
          end if
       enddo
       !call dshr_strdata_print(sdat, trim(file)//' NML_ONLY')

    endif   ! master_task

    call shr_mpi_bcast(sdat%dataMode      ,mpicom,'dataMode'    )
    call shr_mpi_bcast(sdat%calendar      ,mpicom,'calendar'    )
    call shr_mpi_bcast(sdat%nstreams      ,mpicom,'nstreams'    )
    call shr_mpi_bcast(sdat%nvectors      ,mpicom,'nvectors'    )
    call shr_mpi_bcast(sdat%vector_names  ,mpicom,'vectors'     )
    call shr_mpi_bcast(sdat%streamfiles   ,mpicom,'streamfiles' )
    call shr_mpi_bcast(sdat%taxMode       ,mpicom,'taxMode'     )
    call shr_mpi_bcast(sdat%dtlimit       ,mpicom,'dtlimit'     )
    call shr_mpi_bcast(sdat%mapalgo       ,mpicom,'mapalgo'     )
    call shr_mpi_bcast(sdat%tintalgo      ,mpicom,'tintalgo'    )
    call shr_mpi_bcast(sdat%readmode      ,mpicom,'readmode'    )

    sdat%ymdLB    = -1
    sdat%todLB    = -1
    sdat%ymdUB    = -1
    sdat%todUB    = -1
    sdat%dtmin    = 1.0e30
    sdat%dtmax    = 0.0
    sdat%eccen    = SHR_ORB_UNDEF_REAL
    sdat%mvelpp   = SHR_ORB_UNDEF_REAL
    sdat%lambm0   = SHR_ORB_UNDEF_REAL
    sdat%obliqr   = SHR_ORB_UNDEF_REAL
    sdat%modeldt  = 0

    if (my_task == master_task) then
       call shr_strdata_print(sdat,'SDAT data from '//trim(nml_filename))
    endif

  end subroutine dshr_strdata_create_from_cap

  !===============================================================================

  subroutine dshr_strdata_init(mpicom, my_task, master_task, logunit, sdat, rc)

    ! Note: for the NUOPC implementation the model MUST HAVE its own mesh file
    ! and the model domain CANNOT come from the first stream file

    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat

    ! input/output variables
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: my_task
    integer                , intent(in)    :: master_task
    integer                , intent(in)    :: logunit
    type(shr_strdata_type) , intent(inout) :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)         :: distgrid
    integer                     :: dimcount
    integer                     :: tilecount
    integer, allocatable        :: minIndexPTile(:,:)
    integer, allocatable        :: maxIndexPTile(:,:)
    integer                     :: lnx, lny   ! global mesh dimensions
    integer                     :: ns         ! number of local mesh elements
    integer                     :: n,m,k      ! generic index
    character(CL)               :: filePath   ! generic file path
    character(CL)               :: fileName   ! generic file name
    integer                     :: nfiles     ! number of data files for a given stream
    character(CXX)              :: fldList    ! list of fields
    character(CS)               :: uname      ! u vector field name
    character(CS)               :: vname      ! v vector field name
    integer                     :: nu, nv     ! vector indices
    integer                     :: nstream    ! loop stream index
    integer                     :: nvector    ! loop vector index
    integer                     :: nfld       ! loop stream field index 
    character(CS)               :: fldname    ! field name of field index nfld
    integer                     :: nflds      ! total number of fields in a given stream
    type(ESMF_Field)            :: lfield     ! temporary
    type(ESMF_Field)            :: lfield_src ! temporary
    type(ESMF_Field)            :: lfield_dst ! temporary
    integer                     :: srcTermProcessing_Value = 0 ! should this be a module variable?
    integer          ,pointer   :: dof(:)
    character(CS)               :: tmpstr
    character(*)     ,parameter :: F00 = "('(shr_strdata_init) ',8a)"
    character(len=*), parameter :: subname = "(shr_strdata_init) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! ------------------------------------
    ! Update sdat streams info - assumes that readnml has already been called for the shr_strdata namelist
    ! ------------------------------------

    if (my_task == master_task) then

       do n=1,nStrMax
          ! Determine the number of streams
          if (trim(sdat%streamfiles(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nstreams = max(sdat%nstreams,n)
          end if
          ! Determine number of streams - first check if a filename is defined in the stream
          call dshr_stream_getNFiles(sdat%streams(n), nfiles)
          if (nfiles > 0) then
             sdat%nstreams = max(sdat%nstreams,n)
          end if
          if (trim(sdat%taxMode(n)) == trim(shr_stream_taxis_extend)) then
             sdat%dtlimit(n) = 1.0e30
          end if
       end do
       ! Determine the number of vectors to map
       sdat%nvectors = 0
       do n=1,nVecMax
          if (trim(sdat%vector_names(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nvectors = n
          end if
       end do
    endif
    call shr_mpi_bcast(sdat%nstreams  ,mpicom,'nstreams')
    call shr_mpi_bcast(sdat%nvectors  ,mpicom,'nvectors')
    call shr_mpi_bcast(sdat%dtlimit   ,mpicom,'dtlimit')

    ! ------------------------------------
    ! Loop through the streams
    ! ------------------------------------

    do n = 1,sdat%nstreams

       ! ------------------------------------
       ! Create the target stream mesh from the stream mesh file
       ! ------------------------------------
       call dshr_stream_getMeshFileName (sdat%streams(n), filename)
       sdat%mesh_streams(n) = ESMF_MeshCreate(trim(filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! TODO: add functionality if the stream mesh needs to be created from a grid

       ! ------------------------------------
       ! Initialize the stream pio decomp
       ! ------------------------------------
       ! Determine lnx and lny - since are assuming a mesh, lny=1 and lnx is the total number of mesh elements
       call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       lnx = maxval(maxIndexPTile)
       lny = 1

       ! determine the dof needed to initialize pio
       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns, size(dof), dof(1), dof(ns)  !minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       ! Initialize the stream pio decomp
       call pio_initdecomp(sdat%pio_subsystem, pio_double, (/lnx, lny/), dof, sdat%pio_iodesc(n))
       deallocate(dof)

       ! ------------------------------------
       ! broadcast the stream calendar
       ! ------------------------------------
       call shr_mpi_bcast(sdat%streams(n)%calendar, mpicom)

       ! ------------------------------------
       ! Create sdat field bundles and coszen field
       ! ------------------------------------

       ! create empty field bundles
       sdat%FB_stream_lbound(n) = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at lower time bound
       sdat%FB_stream_ubound(n) = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at upper time bound
       sdat%FB_model_lbound(n)  = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%FB_model_lbound(n)  = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%FB_model(n)         = ESMF_FieldBundleCreate(rc=rc) ! time interpolation on model mesh

       ! create colon deliminted string of the model names corresponding to the stream data field names
       call dshr_stream_getModelFieldList(sdat%streams(n), fldList)

       ! loop over field names in fldList
       nflds = shr_string_listGetNum(fldList)
       do nfld = 1,nflds
          ! get nth fldname in colon delimited string
          call shr_string_listGetName(fldlist, nfld, fldname)

          ! create temporary field with name fldname on stream mesh
          ! add the field to the lbound and ubound field bundle
          lfield = ESMF_FieldCreate(sdat%mesh_streams(n), ESMF_TYPEKIND_R8, name=trim(fldname), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%FB_stream_lbound(n), (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%FB_stream_ubound(n), (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! create temporary field with name fldname on model mesh
          ! add the field to the lbound and ubound field bundle as well as the time interpolated field bundle
          lfield = ESMF_FieldCreate(sdat%mesh_model, ESMF_TYPEKIND_R8, name=trim(fldname), &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%FB_model_lbound(n), (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%FB_model_ubound(n), (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldBundleAdd(sdat%FB_model(n)       , (/lfield/), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do

       ! create a field for coszen time interpolation for this stream if needed
       if (trim(sdat%tintalgo(n)) == 'coszen') then
          sdat%field_coszen = ESMF_FieldCreate(sdat%mesh_streams(n), ESMF_TYPEKIND_R8, name='tavCosz', &
               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       endif

    end do

    ! ------------------------------------
    ! Create the sdat route handles
    ! ------------------------------------
    ! create the source and destination fields needed for the route handles
    ! these fields will be used to create the route handles
    ! since all fields in a stream share the same mesh and there is only a unique model mesh
    ! can do this outside of a stream loop by just using the first stream index

    call FB_getFieldN(sdat%FB_stream_lbound(1), 1, lfield_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFieldN(sdat%FB_model_lbound(1), 1, lfield_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Loop over all active streams and create bilinear route handle for stream(n) -> model mapping
    do n = 1,sdat%nstreams
       if (sdat%mapalgo(n) == 'bilinear') then
          call ESMF_FieldRegridStore(lfield_src, lfield_dst, routehandle=sdat%RH_stream2model(mapfcopy,n), &
               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, polemethod=ESMF_POLEMETHOD_ALLAVG, &
               extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
              !srcTermProcessing=srcTermProcessing_Value, ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
               srcTermProcessing=srcTermProcessing_Value, ignoreDegenerate=.true., rc=rc)
       end if
    end do ! end of loop over streams

    ! ------------------------------------
    ! create the sdat ustrm and vstrm fields
    ! ------------------------------------
    ! check vectors and compute ustrm,vstrm
    ! determine if vector field names match any field names in the input stream

    do m = 1,sdat%nvectors
       ! check that vector field list is a valid colon delimited string
       if (.not. shr_string_listIsValid(sdat%vector_names(m))) then
          write(logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(sdat%vector_names(m))
          call shr_sys_abort(subname//': vec fldlist invalid:'//trim(sdat%vector_names(m)))
       endif
       ! check that only 2 fields are contained for any vector pairing
       if (shr_string_listGetNum(sdat%vector_names(m)) /= 2) then
          write(logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(sdat%vector_names(m))
          call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(sdat%vector_names(m)))
       endif

       ! Determine stream indices that vector_names
       nu = 0
       nv = 0
       ! get name of the first and second field in the colon delimited string
       call shr_string_listGetName(sdat%vector_names(m),1,uname)
       call shr_string_listGetName(sdat%vector_names(m),2,vname)

       ! loop through the streams and find which stream(s) contain the vector field pair
       ! normally both fields in the pair will be on one stream
       do n = 1,sdat%nstreams
          if (FB_fldchk(sdat%FB_stream_lbound(n), trim(uname), rc=rc)) nu = n
          if (FB_fldchk(sdat%FB_stream_lbound(n), trim(vname), rc=rc)) nv = n
       enddo
       if (nu == 0  .or. nv == 0) then
          ! if the input fields are not contained - then abort
          write(logunit,*) trim(subname),' vec flds not found  m=',m,trim(sdat%vector_names(m))
          call shr_sys_abort(subname//': vec flds not found:'//trim(sdat%vector_names(m)))
       else if (nu /= nv) then
          ! TODO: a grid comparison was made for the vector fields before - is this necessary?
          write(logunit,*) trim(subname),' vec fld doms not same m=',m,trim(sdat%vector_names(m))
          call shr_sys_abort(subname//': vec fld doms not same:'//trim(sdat%vector_names(m)))
       else
          sdat%ustrm(m)= nu
          sdat%vstrm(m)= nv
       end if
    enddo

  end subroutine dshr_strdata_init

  !===============================================================================

  subroutine dshr_strdata_advance(sdat, ymd, tod, mpicom, logunit, istr, timers, my_task, master_task, rc)

    ! -------------------------------------------------------
    ! Mismatching calendars: 4 cases
    ! (0) The stream calendar and model calendar are identical
    ! (1) The stream is a no leap calendar and the model is gregorian
    ! (2) The stream is a gregorian calendar and the model is a noleap calendar
    ! (3) The calendars mismatch and none of the above
    ! -------------------------------------------------------
    !
    ! ymdmod and todmod are the ymd and tod to time interpolate to.
    ! Generally, these are just the model date and time.  Also, always
    ! use the stream calendar for time interpolation for reasons
    ! described below.  When there is a calendar mismatch, support Feb
    ! 29 in a special way as needed to get reasonable values.  Note
    ! that when Feb 29 needs to be treated specially, a discontinuity
    ! will be introduced.  The size of that discontinuity will depend
    ! on the time series input data.
    !
    ! (0) The stream calendar and model calendar are identical:
    ! Proceed in the standard way.
    !
    ! (1) The stream is a no leap calendar and the model is gregorian:
    ! Time interpolate on the noleap calendar.  Then if the model date
    ! is Feb 29, compute stream data for Feb 28 by setting ymdmod and
    ! todmod to Feb 28.  This results in duplicate stream data on Feb
    ! 28 and Feb 29 and a discontinuity at the start of Feb 29.  This
    ! could be potentially updated by using the gregorian calendar for
    ! time interpolation when the input data is relatively infrequent
    ! (say greater than daily) with the following concerns.
    !   - The forcing will not be reproduced identically on
    !     the same day with climatological inputs data
    !   - Input data with variable input frequency might behave funny
    !   - An arbitrary discontinuity will be introduced in the time 
    !     interpolation method based upon the logic chosen to transition
    !     from reproducing Feb 28 on Feb 29 and interpolating to Feb 29.
    !   - The time gradient of data will change by adding a day arbitrarily.
    !
    ! (2) The stream is a gregorian calendar and the model is a noleap calendar:
    ! Then just time interpolate on the gregorian calendar. This
    ! causes Feb 29 stream data to be skipped and will lead to a
    ! discontinuity at the start of March 1.
    !
    ! (3) If the calendars mismatch and neither of the three cases above
    ! are recognized, then abort.
    ! -------------------------------------------------------

    ! input/output variables
    type(shr_strdata_type) ,intent(inout)       :: sdat
    integer                ,intent(in)          :: ymd    ! current model date
    integer                ,intent(in)          :: tod    ! current model date
    integer                ,intent(in)          :: mpicom
    integer                ,intent(in)          :: logunit
    integer                ,intent(in)          :: my_task 
    integer                ,intent(in)          :: master_task 
    character(len=*)       ,intent(in),optional :: istr
    logical                ,intent(in),optional :: timers
    integer                ,intent(out)         :: rc 

    ! local variables
    integer                     :: n,m,i,kf             ! generic index
    logical ,allocatable        :: newData(:)
    integer                     :: ierr
    integer                     :: nu,nv
    integer                     :: lsize
    integer ,allocatable        :: ymdmod(:)            ! modified model dates to handle Feb 29
    integer                     :: todmod               ! modified model dates to handle Feb 29
    character(len=32)           :: lstr
    logical                     :: ltimers
    real(R8)                    :: flb,fub              ! factor for lb and ub
    real(R8) ,pointer           :: lon_model(:)         ! lon radians
    real(R8) ,pointer           :: lat_model(:)         ! lat radians
    real(R8) ,pointer           :: cosz(:)              ! cosz
    real(R8) ,pointer           :: tavCosz(:)           ! cosz, time avg over [LB,UB]
    real(R8) ,pointer           :: dataptr(:)           ! pointer into field bundle
    real(R8) ,pointer           :: dataptr_lbound(:)    ! pointer into field bundle
    type(ESMF_Time)             :: timeLB, timeUB       ! lb and ub times
    type(ESMF_TimeInterval)     :: timeint              ! delta time
    integer                     :: dday                 ! delta days
    real(R8)                    :: dtime                ! delta time
    integer                     :: uvar,vvar
    integer                     :: year,month,day       ! date year month day
    integer                     :: spatialDim
    integer                     :: numOwnedElements     ! local size of mesh
    character(CS)               :: uname                ! u vector field name
    character(CS)               :: vname                ! v vector field name
    real(r8), pointer           :: nu_coords(:)         ! local element mesh lat and lons
    real(r8), pointer           :: nv_coords(:)         ! local element mesh lat and lons
    real(r8), pointer           :: model_coords(:)      ! local element mesh lat and lons
    type(ESMF_Field)            :: field_src
    type(ESMF_Field)            :: field_dst
    real(r8), pointer           :: data2d_src(:,:)
    real(r8), pointer           :: data2d_dst(:,:)
    real(r8), pointer           :: data_u_src(:)
    real(r8), pointer           :: data_v_src(:)
    real(r8), pointer           :: data_u_dst(:)
    real(r8), pointer           :: data_v_dst(:)
    real(r8)                    :: lon, lat
    real(r8)                    :: sinlon, sinlat
    real(r8)                    :: coslon, coslat
    real(r8)                    :: ux, uy, uz 
    logical                     :: checkflag = .false.
    integer                     :: fieldCount
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    real(R8)         ,parameter :: solZenMin = 0.001_R8 ! minimum solar zenith angle
    integer          ,parameter :: tadj = 2
    character(len=*) ,parameter :: timname = "_strd_adv"
    character(*)     ,parameter :: subname = "(shr_strdata_advance) "
    !-------------------------------------------------------------------------------

    if (sdat%nstreams < 1) return

    lstr = ''
    if (present(istr)) lstr = trim(istr)
    ltimers = .true.
    if (present(timers)) ltimers = timers

    if (.not.ltimers) call t_adj_detailf(tadj)
    call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
    call t_startf(trim(lstr)//trim(timname)//'_total')

    sdat%ymd = ymd
    sdat%tod = tod

    if (sdat%nstreams > 0) then
       allocate(newData(sdat%nstreams))
       allocate(ymdmod(sdat%nstreams))

       do n = 1,sdat%nstreams

          ! ---------------------------------------------------------
          ! Consistency checks
          ! ---------------------------------------------------------

          ! case(0)
          ymdmod(n) = ymd
          todmod    = tod
          if (trim(sdat%calendar) /= trim(sdat%streams(n)%calendar)) then
             if ( (trim(sdat%calendar) == trim(shr_cal_gregorian)) .and. &
                  (trim(sdat%streams(n)%calendar) == trim(shr_cal_noleap))) then

                ! case (1), set feb 29 = feb 28
                call shr_cal_date2ymd (ymd,year,month,day)
                if (month == 2 .and. day == 29) then
                   call shr_cal_ymd2date(year,2,28,ymdmod(n))
                endif

             else if ((trim(sdat%calendar) == trim(shr_cal_noleap)) .and. &
                      (trim(sdat%streams(n)%calendar) == trim(shr_cal_gregorian))) then
                ! case (2), feb 29 input data will be skipped automatically

             else
                ! case (3), abort
                write(logunit,*) trim(subname),' ERROR: mismatch calendar ', &
                     trim(sdat%calendar),':',trim(sdat%streams(n)%calendar)
                call shr_sys_abort(trim(subname)//' ERROR: mismatch calendar ')
             endif
          endif

          ! ---------------------------------------------------------
          ! Determine if new data is read in - if so
          ! Copy FB_stream_ubound to FB_stream_lbound
          ! Read in new FB_stream_ubound data
          ! ---------------------------------------------------------

          call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
          call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

          call dshr_strdata_readLBUB(sdat%streams(n),                  &
               sdat%pio_subsystem, sdat%io_type, sdat%pio_iodesc(n),   &
               ymdmod(n), todmod, mpicom, my_task, master_task,        &
               sdat%FB_stream_lbound(n), sdat%ymdLB(n), sdat%todLB(n), &
               sdat%FB_stream_ubound(n), sdat%ymdUB(n), sdat%todUB(n), &
               sdat%FB_stream_alltimes(:,n),                           &
               trim(sdat%readmode(n)), newData(n), istr=trim(lstr)//'_readLBUB')

          if (debug > 0) then
             write(logunit,*) trim(subname),' newData flag = ',n,newData(n)
             write(logunit,*) trim(subname),' LB ymd,tod = ',n,sdat%ymdLB(n),sdat%todLB(n)
             write(logunit,*) trim(subname),' UB ymd,tod = ',n,sdat%ymdUB(n),sdat%todUB(n)
          endif

          ! ---------------------------------------------------------
          ! Reset time bounds if newdata read in
          ! ---------------------------------------------------------

          ! If new data was read in - then set sdat%ymdLB, sdat%ymdUB, sdat%dtmin and sdat%dtmax
          if (newData(n)) then
             if (debug > 0) then
                call FB_diagnose(sdat%FB_stream_lbound(n), subname//':FB_stream_lbound', rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call FB_diagnose(sdat%FB_stream_ubound(n), subname//':FB_stream_ubound', rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif

             call shr_cal_date2ymd(sdat%ymdLB(n),year,month,day)
             call shr_cal_timeSet(timeLB, sdat%ymdLB(n), 0, sdat%streams(n)%calendar)
             call shr_cal_timeSet(timeUB, sdat%ymdUB(n), 0, sdat%streams(n)%calendar)
             timeint = timeUB-timeLB
             call ESMF_TimeIntervalGet(timeint, StartTimeIn=timeLB, d=dday)
             dtime = abs(real(dday, R8) + real(sdat%todUB(n)-sdat%todLB(n), R8)/shr_const_cDay)

             sdat%dtmin(n) = min(sdat%dtmin(n), dtime)
             sdat%dtmax(n) = max(sdat%dtmax(n), dtime)
             if ((sdat%dtmax(n)/sdat%dtmin(n)) > sdat%dtlimit(n)) then
                write(logunit,*) trim(subname),' ERROR: for stream ',n
                write(logunit,*) trim(subName),' ERROR: dt limit1 ',sdat%dtmax(n),sdat%dtmin(n),sdat%dtlimit(n)
                write(logunit,*) trim(subName),' ERROR: dt limit2 ',sdat%ymdLB(n),sdat%todLB(n),sdat%ymdUB(n),sdat%todUB(n)
                call shr_sys_abort(trim(subName)//' ERROR dt limit for stream')
             endif
          endif
          call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')

          ! ---------------------------------------------------------
          ! Do spatial interpolation if newdata read in
          ! ---------------------------------------------------------
          ! TODO: generalize the bilinear mapping below

          ! If new data was read in, interpolate the lower and upper bound data to the model grid
          if (newData(n)) then
             call t_startf(trim(lstr)//trim(timname)//'_map')
             call FB_Regrid(sdat%FB_stream_lbound(n), sdat%FB_model_lbound(n), sdat%rh_stream2model(mapbilnr,n), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_Regrid(sdat%FB_stream_ubound(n), sdat%FB_model_ubound(n), sdat%rh_stream2model(mapbilnr,n), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call t_stopf(trim(lstr)//trim(timname)//'_map')
             if (debug > 0) then
                call FB_diagnose(sdat%FB_model_lbound(n), subname//':FB_model_lbound',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
          endif

       enddo ! end of loop over streams(n)

       ! ---------------------------------------------------------
       ! remap with vectors if needed
       ! ---------------------------------------------------------

       do m = 1,sdat%nvectors
          nu = sdat%ustrm(m) ! nu is the stream index that contains the u vector
          nv = sdat%vstrm(m) ! nv is the stream index that contains the v vector

          ! TODO: this is not correct logic - need to change it
          if ((nu > 0 .or. nv > 0) .and. (newdata(nu) .or. newdata(nv))) then

             call t_startf(trim(lstr)//trim(timname)//'_vect')

             ! get nu coords
             call ESMF_MeshGet(sdat%mesh_streams(nu), spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(nu_coords(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%mesh_streams(nu), ownedElemCoords=nu_coords)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get nv coords
             call ESMF_MeshGet(sdat%mesh_streams(nv), spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(nv_coords(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%mesh_streams(nv), ownedElemCoords=nv_coords)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get model coords
             call ESMF_MeshGet(sdat%mesh_model, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(model_coords(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%mesh_model, ownedElemCoords=model_coords)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! create source field and destination fields
             ! TODO: assume that the two meshes are idential - but need to confirm this
             field_src = ESMF_FieldCreate(sdat%mesh_streams(nu), ESMF_TYPEKIND_R8, name='field_src', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             field_dst = ESMF_FieldCreate(sdat%mesh_model, ESMF_TYPEKIND_R8, name='field_dst', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! set pointers to source and destination data that will be filled in with rotation to cart3d
             call ESMF_FieldGet(field_src, farrayPtr=data2d_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(field_dst, farrayPtr=data2d_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get names of vector field pairs
             call shr_string_listGetName(sdat%vector_names(m), 1, uname)
             call shr_string_listGetName(sdat%vector_names(m), 2, vname)

             ! map lower bounds: rotate source data, regrid, then rotate back
             call FB_getFldPtr(sdat%FB_stream_lbound(nu), trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound(nu), trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound(nv), trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound(nv), trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_src)
                lon = nu_coords(2*n-1)
                lat = nu_coords(2*n)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                data2d_src(1,n) = coslon * data_u_src(n) - sinlon * data_v_src(n)
                data2d_src(2,n) = sinlon * data_u_src(n) + coslon * data_v_src(n)
             enddo
             call ESMF_FieldRegrid(field_src, field_dst, sdat%RH_stream2model(mapbilnr,nu), &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_dst)
                lon = model_coords(2*n-1)
                lat = model_coords(2*n)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,n)
                uy = data2d_dst(2,n)
                uz = data2d_dst(3,n)
                data_u_dst(n) =  coslon * data_u_dst(n) + sinlon * data_v_dst(n)
                data_v_dst(n) = -sinlon * data_u_dst(n) + coslon * data_v_dst(n)
             enddo

             ! map upper bounds: rotate source data, regrid, then rotate back
             call FB_getFldPtr(sdat%FB_stream_ubound(nu), trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound(nu), trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound(nv), trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound(nv), trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_src)
                lon = nu_coords(2*n-1)
                lat = nu_coords(2*n)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                data2d_src(1,n) = coslon * data_u_src(n) - sinlon * data_v_src(n)
                data2d_src(2,n) = sinlon * data_u_src(n) + coslon * data_v_src(n)
             enddo
             call ESMF_FieldRegrid(field_src, field_dst, sdat%RH_stream2model(mapbilnr,nu), &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_dst)
                lon = model_coords(2*n-1)
                lat = model_coords(2*n)
                sinlon = sin(lon*deg2rad); coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad); coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,n)
                uy = data2d_dst(2,n)
                uz = data2d_dst(3,n)
                data_u_dst(n) =  coslon * data_u_dst(n) + sinlon * data_v_dst(n)
                data_v_dst(n) = -sinlon * data_u_dst(n) + coslon * data_v_dst(n)
             enddo

             call t_stopf(trim(lstr)//trim(timname)//'_vect')
          endif
       enddo

       ! ---------------------------------------------------------
       ! Do time interpolation to create FB_model
       ! ---------------------------------------------------------

          ! get model mesh lat/lon data
          call ESMF_MeshGet(sdat%mesh_model, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          call ESMF_MeshGet(mesh_model, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(lon_model(numOwnedElements))
          allocate(lat_model(numOwnedElements))
          do n = 1, numOwnedElements
             lon_model(n) = ownedElemCoords(2*n-1)
             lat_model(n) = ownedElemCoords(2*n)
          end do
          deallocate(ownedElementCoords)

       do n = 1,sdat%nstreams

          ! Get field namelist
          call ESMF_FieldBundleGet(sdat%FB_model(n), fieldCount=fieldCount, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          allocate(lfieldNameList(fieldCount))
          call ESMF_FieldBundleGet(sdat%FB_model, fieldNameList=lfieldNameList, rc=rc)


          ! allocate cosz and avgcosz arrays
          allocate(tavCosz(numOwnedElements),cosz(numOwnedElements))

          if (trim(sdat%tintalgo(n)) == 'coszen') then

             ! ------------------------------------------
             ! time interpolation method is coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_coszen')

             ! make sure orb info has been set
             if (sdat%eccen == SHR_ORB_UNDEF_REAL) then
                call shr_sys_abort(subname//' ERROR in orb params for coszen tinterp')
             else if (sdat%modeldt < 1) then
                call shr_sys_abort(subname//' ERROR: model dt < 1 for coszen tinterp')
             endif

             ! get cosz factor
             call t_startf(trim(lstr)//trim(timname)//'_coszenC')
             cosz(:) = 0.0_r8
             call dshr_tInterp_getCosz(cosz, lon_model, lat_model, ymdmod(n), todmod, &
                  sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%streams(n)%calendar)
             call t_stopf(trim(lstr)//trim(timname)//'_coszenC')

             ! get avg cosz factor
             if (newdata(n)) then
                ! compute a new avg cosz
                call t_startf(trim(lstr)//trim(timname)//'_coszenN')
                call shr_tInterp_getAvgCosz(tavCosz, lon_model, lat_model,  &
                     sdat%ymdLB(n), sdat%todLB(n),  sdat%ymdUB(n), sdat%todUB(n),  &
                     sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%modeldt, &
                     sdat%streams(n)%calendar)
                sdat%avgCoszen(:) = tavgCosz(:)
                call t_stopf(trim(lstr)//trim(timname)//'_coszenN')
             else
                ! reuse existing avg cosz
                tavgCosz(:) = sdat%avgCoszen(:)
             endif

             ! compute time interperpolate value - LB data normalized with this factor: cosz/tavCosz
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                call FB_getfldptr(sdat%FB_model_lbound, fieldNameList(n), dataptr_lb, rc=rc)
                do i = 1,lsize
                   if (cosz(i) > solZenMin) then
                      dataptr(i) = dataptr_lb(i)*cosz(i)/sdat%tavCosz(i)
                   else
                      dataptr(i) = 0._r8
                   endif
                end do
             end do

             ! deallocate memory
             call t_stopf(trim(lstr)//trim(timname)//'_coszen')

          elseif (trim(sdat%tintalgo(n)) /= trim(shr_strdata_nullstr)) then

             ! ------------------------------------------
             ! time interpolation method is not coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_tint')
             call shr_tInterp_getFactors(&
                  sdat%ymdlb(n), sdat%todlb(n), sdat%ymdub(n), sdat%todub(n), &
                  ymdmod(n),todmod,flb,fub, &
                  calendar=sdat%streams(n)%calendar,algo=trim(sdat%tintalgo(n)))
             if (debug > 0) then
                write(logunit,*) trim(subname),' interp = ',n,flb,fub
             endif
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                call FB_getfldptr(sdat%FB_model_lbound, fieldNameList(n), dataptr_lb, rc=rc)
                call FB_getfldptr(sdat%FB_model_ubound, fieldNameList(n), dataptr_ub, rc=rc)
                dataptr(:) = dataptr_lb(:) * flb + dataptr_ub(:) * fub
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_tint')

          else

             ! ------------------------------------------
             ! zero out stream data for this field
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_zero')
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                dataptr(:) = 0._r8
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_zero')

          endif
          deallocate(lfieldNameList(fieldCount))

          if (debug > 0) then
             ! TODO: call field bundle diagnose here
          endif

       enddo
       deallocate(newData)
       deallocate(ymdmod)

    endif    ! nstreams > 0

    deallocate(tavCosz,cosz,lonr,latr)

    call t_stopf(trim(lstr)//trim(timname)//'_total')
    if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine dshr_strdata_advance

  !===============================================================================

  subroutine dshr_strdata_restWrite(filename,sdat,mpicom,str1,str2)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: sdat
    integer               ,intent(in)    :: mpicom
    character(len=*)      ,intent(in)    :: str1
    character(len=*)      ,intent(in)    :: str2

    !--- local ----
    integer     :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restWrite) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    if (my_task == 0) then
       call dshr_stream_restWrite(sdat%streams, trim(filename), trim(str1), trim(str2), sdat%nstreams)
    endif

  end subroutine dshr_strdata_restWrite

  !===============================================================================

  subroutine dshr_strdata_restRead(filename,sdat,mpicom)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: sdat
    integer               ,intent(in)    :: mpicom

    !--- local ----
    integer     :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restRead) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    if (my_task == 0) then
       call dshr_stream_restRead(sdat%streams, trim(filename), sdat%nstreams)
    endif

  end subroutine dshr_strdata_restRead

  !===============================================================================

  subroutine dshr_strdata_setOrbs(sdat,eccen,mvelpp,lambm0,obliqr,modeldt)

    ! input/output variables
    type(shr_strdata_type) ,intent(inout) :: sdat
    real(R8)               ,intent(in)    :: eccen
    real(R8)               ,intent(in)    :: mvelpp
    real(R8)               ,intent(in)    :: lambm0
    real(R8)               ,intent(in)    :: obliqr
    integer                ,intent(in)    :: modeldt
    !-------------------------------------------------------------------------------

    sdat%eccen   = eccen
    sdat%mvelpp  = mvelpp
    sdat%lambm0  = lambm0
    sdat%obliqr  = obliqr
    sdat%modeldt = modeldt

  end subroutine dshr_strdata_setOrbs

  !===============================================================================

  subroutine dshr_strdata_print(sdat, name, logunit)

    !  Print strdata common to all data models

    ! !input/output parameters
    type(shr_strdata_type)    ,intent(in) :: sdat  ! strdata data data-type
    character(len=*),optional ,intent(in) :: name  ! just a name for tracking
    integer                   ,intent(in) :: logunit

    ! local variables
    integer                :: n
    character(CL)          :: lname
    character(*),parameter :: F00 = "('(shr_strdata_print) ',8a)"
    character(*),parameter :: F01 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter :: F02 = "('(shr_strdata_print) ',a,es13.6)"
    character(*),parameter :: F03 = "('(shr_strdata_print) ',a,l6)"
    character(*),parameter :: F04 = "('(shr_strdata_print) ',a,i2,a,a)"
    character(*),parameter :: F05 = "('(shr_strdata_print) ',a,i2,a,i6)"
    character(*),parameter :: F06 = "('(shr_strdata_print) ',a,i2,a,l2)"
    character(*),parameter :: F07 = "('(shr_strdata_print) ',a,i2,a,es13.6)"
    character(*),parameter :: F20 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter :: F90 = "('(shr_strdata_print) ',58('-'))"
    character(*),parameter :: subName = "(shr_strdata_print) "
    !-------------------------------------------------------------------------------

    lname = 'unknown'
    if (present(name)) then
       lname = trim(name)
    endif

    !----------------------------------------------------------------------------
    ! document datatype settings
    !----------------------------------------------------------------------------
    write(logunit,F90)
    write(logunit,F00) "name        = ",trim(lname)
    write(logunit,F00) "dataMode    = ",trim(sdat%dataMode)
    write(logunit,F00) "calendar    = ",trim(sdat%calendar)
    write(logunit,F01) "io_type     = ",sdat%io_type
    write(logunit,F02) "eccen       = ",sdat%eccen
    write(logunit,F02) "mvelpp      = ",sdat%mvelpp
    write(logunit,F02) "lambm0      = ",sdat%lambm0
    write(logunit,F02) "obliqr      = ",sdat%obliqr
    write(logunit,F01) "nstreams    = ",sdat%nstreams
    write(logunit,F01) "pio_iotype  = ",sdat%io_type

    do n=1, sdat%nstreams
       write(logunit,F04) "  streamfiles (",n,") = ",trim(sdat%streamfiles(n))
       write(logunit,F04) "  taxMode (",n,")     = ",trim(sdat%taxMode(n))
       write(logunit,F07) "  dtlimit (",n,")     = ",sdat%dtlimit(n)
       write(logunit,F06) "  domaps  (",n,")     = ",sdat%domaps(n)
       write(logunit,F04) "  mapalgo (",n,")     = ",trim(sdat%mapalgo(n))
       write(logunit,F04) "  tintalgo(",n,")     = ",trim(sdat%tintalgo(n))
       write(logunit,F04) "  readmode(",n,")     = ",trim(sdat%readmode(n))
       write(logunit,F01) " "
    end do
    write(logunit,F01) "nvectors    = ",sdat%nvectors
    do n=1, sdat%nvectors
       write(logunit,F04) "  vectors (",n,") = ",trim(sdat%vectors(n))
    end do
    write(logunit,F90)
    call shr_sys_flush(logunit)

  end subroutine dshr_strdata_print

  !===============================================================================

  subroutine dshr_strdata_readLBUB(stream, &
       pio_subsystem, pio_iotype, pio_iodesc, &
       mDate, mSec, mpicom, my_task, master_task, &
       FB_stream_lbound, mDateLB, mSecLB, &
       FB_stream_ubound, mDateUB, mSecUB, &
       FB_stream_alltimes, &
       readMode, newData, rmOldFile, istr)


    ! input/output variables
    type(shr_stream_streamType)   ,intent(inout) :: stream
    type(iosystem_desc_t), target ,intent(inout) :: pio_subsystem
    integer                       ,intent(in)    :: pio_iotype
    type(io_desc_t)               ,intent(inout) :: pio_iodesc
    integer                       ,intent(in)    :: mDate  ,mSec
    integer                       ,intent(in)    :: mpicom
    integer                       ,intent(in)    :: my_task
    integer                       ,intent(in)    :: master_task
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_lbound
    integer                       ,intent(inout) :: mDateLB,mSecLB
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_ubound
    integer                       ,intent(inout) :: mDateUB,mSecUB
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_alltimes(:)
    character(len=*)              ,intent(in)    :: readMode
    logical                       ,intent(out)   :: newData
    logical          ,optional    ,intent(in)    :: rmOldFile
    character(len=*) ,optional    ,intent(in)    :: istr

    ! local variables
    integer           :: ierr       ! error code
    integer           :: rCode      ! return code
    logical           :: localCopy,fileexists
    integer           :: ivals(6)   ! bcast buffer
    integer           :: oDateLB,oSecLB,dDateLB
    integer           :: oDateUB,oSecUB,dDateUB
    real(R8)          :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer           :: n_lb, n_ub
    character(CL)     :: fn_lb,fn_ub,fn_next,fn_prev
    character(CL)     :: path
    character(len=32) :: lstr
    real(R8)          :: spd

    character(*), parameter :: subname = '(dshr_strdata_readLBUB) '
    character(*), parameter :: F00   = "('(dshr_strdata_readLBUB) ',8a)"
    character(*), parameter :: F01   = "('(dshr_strdata_readLBUB) ',a,5i8)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Read LB and UB stream data
    !----------------------------------------------------------------------------

    lstr = 'dshr_strdata_readLBUB'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    call t_startf(trim(lstr)//'_setup')

    spd = shr_const_cday

    newData = .false.
    n_lb = -1
    n_ub = -1
    fn_lb = 'undefinedlb'
    fn_ub = 'undefinedub'

    oDateLB = mDateLB
    oSecLB  = mSecLB
    oDateUB = mDateUB
    oSecUB  = mSecUB

    rDateM  = real(mDate  ,R8) + real(mSec  ,R8)/spd
    rDateLB = real(mDateLB,R8) + real(mSecLB,R8)/spd
    rDateUB = real(mDateUB,R8) + real(mSecUB,R8)/spd
    call t_stopf(trim(lstr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(lstr)//'_fbound')
       if (my_task == master_task) then
          call dshr_stream_findBounds(stream,mDate,mSec, &
               ivals(1),dDateLB,ivals(2),ivals(5),fn_lb, &
               ivals(3),dDateUB,ivals(4),ivals(6),fn_ub  )
          call dshr_stream_getFilePath(stream,path)
       endif
       call t_stopf(trim(lstr)//'_fbound')

       call t_startf(trim(lstr)//'_bcast')
       call shr_mpi_bcast(stream%calendar,mpicom)
       call shr_mpi_bcast(ivals,mpicom)
       mDateLB = ivals(1)
       mSecLB  = ivals(2)
       mDateUB = ivals(3)
       mSecUB  = ivals(4)
       n_lb    = ivals(5)
       n_ub    = ivals(6)
       call t_stopf(trim(lstr)//'_bcast')
    endif

    if (mDateLB /= oDateLB .or. mSecLB /= oSecLB) then
       newdata = .true.

       if (mDateLB == oDateUB .and. mSecLB == oSecUB) then

          ! copy FB_stream_ubound to FB_stream_lbound
          call t_startf(trim(lstr)//'_LB_copy')
          call ESMF_FieldBundleGet(FB_stream_ubound, fieldCount=fieldCount, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          allocate(lfieldNameList(fieldCount))
          call ESMF_FieldBundleGet(FB_stream_ubound, itemNameList=lfieldNameList, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do n = 1,fieldCount
             fldname = trim(lfieldnamelist(n))
             call ESMF_FieldBundleGet(FB_stream_ubound, fieldName=fldname, field=lfield_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(FB_stream_lbound, fieldName=fldname, field=lfield_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_ub, farrayPtr=dataptr_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_lb, farrayPtr=dataptr_lb, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr_lb(:) = dataptr_ub(:)
          end do
          call t_stopf(trim(lstr)//'_LB_copy')

       else

          select case(readMode)
          case ('single')
             call dshr_strdata_readstrm(stream, FB_stream, &
                  pio_subsystem, pio_iotype, pio_iodesc, mpicom, path, &
                  fn_lb, n_lb, istr=trim(lstr)//'_UB', boundstr = 'ub')
          case ('full_file')
             ! TODO: implement full file read
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
          end select
       endif
    endif

    if (mDateUB /= oDateUB .or. mSecUB /= oSecUB) then
       newdata = .true.

       select case(readMode)
       case ('single')
          call dshr_strdata_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, avUB, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case ('full_file')
          call dshr_strdata_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
               gsMap, avUB, avFile, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case default
          write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
          call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
       end select

    endif

    call t_startf(trim(lstr)//'_filemgt')
    !--- determine previous & next data files in list of files ---
    if (my_task == master_task .and. newdata) then
       call dshr_stream_getFilePath(stream,path)
    endif
    call t_stopf(trim(lstr)//'_filemgt')

  end subroutine dshr_strdata_readLBUB

  !===============================================================================

  subroutine dshr_strdata_readstrm(stream, FB_stream, &
       pio_subsystem, pio_iotype, pio_iodesc, mpicom, my_task, master_task, &
       path, fn, nt, istr, boundstr)

    ! input/output variables
    type(shr_stream_streamType) ,intent(inout)         :: stream
    type(ESMF_FieldBundle)      ,intent(inout)         :: FB_stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)            :: pio_iotype
    type(io_desc_t)             ,intent(inout)         :: pio_iodesc
    integer                     ,intent(in)            :: mpicom
    integer                     ,intent(in)            :: my_task
    integer                     ,intent(in)            :: master_task
    character(len=*)            ,intent(in)            :: path
    character(len=*)            ,intent(in)            :: fn
    integer                     ,intent(in)            :: nt
    character(len=*),optional   ,intent(in)            :: istr
    character(len=*),optional   ,intent(in)            :: boundstr

    ! local variables
    type(ESMF_Field)              :: lfield
    character(CS)                 :: fldname
    character(CL)                 :: fileName
    character(CL)                 :: currfile
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    integer                       :: k, n
    integer                       :: ierr
    logical                       :: fileexists
    integer                       :: rCode      ! return code
    character(len=32)             :: lstr
    logical                       :: fileopen
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(*), parameter :: subname = '(dshr_strdata_readstrm) '
    character(*), parameter :: F00   = "('(dshr_strdata_readstrm) ',8a)"
    character(*), parameter :: F02   = "('(dshr_strdata_readstrm) ',2a,i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_strdata_readstrm'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    ! Set up file to read from
    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    if (my_task == master_task) then
       fileName = trim(path)//fn
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif
    call shr_mpi_bcast(filename,mpicom,'filename')
    call t_stopf(trim(lstr)//'_setup')

    ! Get current file and determine if it is open
    call dshr_stream_getCurrFile(stream, fileopen=fileopen, currfile=currfile, currpioid=pioid)
    if (fileopen .and. currfile==filename) then
       ! don't reopen file, all good
    else
       ! otherwise close the old file if open and open new file
       if (fileopen) then
          if (my_task == master_task) then
             write(logunit,F00) 'close  : ',trim(currfile)
             call shr_sys_flush(logunit)
          endif
          call pio_closefile(pioid)
       endif
       if (my_task == master_task) then
          write(logunit,F00) 'open   : ',trim(filename)
          call shr_sys_flush(logunit)
       endif
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call dshr_stream_setCurrFile(stream, fileopen=.true., currfile=trim(filename), currpioid=pioid)
    endif
    if (my_task == master_task) then
       write(logunit,*) 'file '// trim(filename),nt
    endif

    ! Always use pio to read in stream data
    call t_startf(trim(lstr)//'_readpio')
    if (my_task == master_task) then
       write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename), nt
    endif
    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

    ! Loop over all stream fields in FB_stream
    call ESMF_FieldBundleGet(FB_stream, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FB_stream, itemNameList=lfieldNameList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,fieldCount
       fldname = trim(lfieldnamelist(n))
       ! get varid of field n
       rcode = pio_inq_varid(pioid, fldname, varid)
       ! set frame to time index
       frame = nt
       call pio_setframe(pioid,varid,frame)
       ! set pointer rdata to field data
       call ESMF_FieldBundleGet(FB_stream, fieldName=fldname, field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, farrayPtr=dataptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! read dataptr - which sets FB_stream value for stream field fldname
       call pio_read_darray(pioid, varid, pio_iodesc, dataptr, rcode)
    enddo
    call t_stopf(trim(lstr)//'_readpio')

  end subroutine dshr_strdata_readstrm

end module dshr_strdata_mod
