!*******************************************************************
!*                        PENVOX                                   *
!*                                                                 *
!* Short description:                                              *
!*   Performs the rectilinear transport of particles in a geometry *
!*   that combines quadric surfaces and voxels. Routines from      *
!*   PENGEOM are used to handle the quadric part of the geometry.  *
!*                                                                 *
!* Updates:                                                        *
!*     Adapted by                                                  *
!*       J Sempau                                                  *
!*     Created by                                                  *
!*       A Badal & J Sempau, 2008                                  *
!*       Universitat Politecnica de Catalunya, Barcelona, Spain    *
!* Rev:                                                            *
!*   2014-12-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module geovoxMod
!*******************************************************************
!*    PENVOX vars                                                  *
!*******************************************************************
      implicit none
      save
      character*80 voxfilen
      logical isQuad,isOverlay            ! Indicate: .true. if quadrics defined; .true. if quadrics overlay voxels
      integer granul
      integer, parameter :: maxGranul=100 ! Max num. points in voxel side to compute mass
      integer*4 nx,ny,nz,nxy,nvox,bodymask,matmask
      real*8 dx,dy,dz,volvox,idx,idy,idz,vbb(3)
      real*8 :: memvox=0.0                ! Memory (bytes) used by voxel arrays
      integer, allocatable :: matvox(:)
      real, allocatable :: densvox(:),massvox(:)
      real*8, allocatable :: massvoxd(:)  ! Temporary array to perform integration in DPrec
      end


      module partvoxMod
!*******************************************************************
!*    Particle state vars in the voxelized geo.                    *
!*******************************************************************
      implicit none
      save
      integer*4 xvox,yvox,zvox,absvox
      real*8 uold,vold,ivx,ivy,ivz,sdx,sdy,sdz
      end


      subroutine setIsQuad(value)
!*******************************************************************
!*    Sets quadrics flag in PENVOX module, used to access voxels   *
!*    module without actually declaring its USE.                   *
!*                                                                 *
!*    Input:                                                       *
!*      value: .true. to set to true, else set to .false.          *
!*******************************************************************
      use geovoxMod, only : isQuad
      implicit none
      logical value

      isQuad = value
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine stepx(ds,dsef,ncross)
!*******************************************************************
!*    Substitute for PENGEOM's STEP that transports a particle     *
!*    across a geometry that combines quadric surfaces and voxels. *
!*    The particle will travel across multiple voxels in a single  *
!*    step if the material does not change. The flight stops if:   *
!*    (i) mat changes; (ii) the particle enters or exits the voxels*
!*    bounding box (VBB); (iii) similarly to STEP, when a vacuum   *
!*    gap is found, and the particle stops on the next non-vacuum  *
!*    mat; (iv) the particle escapes to infinity; and (v) the      *
!*    distance ds has been traversed.                              *
!*                                                                 *
!*    Input:                                                       *
!*      ds: distance to travel (cm)                                *
!*    Output:                                                      *
!*      dsef: effective distance travelled (cm)                    *
!*      ncross: 0 if ds has been completed, >0 else                *
!*    Comments:                                                    *
!*      - Values of ncross different from 0 may NOT have the same  *
!*        meaning as in the original STEP.                         *
!*******************************************************************
      use geovoxMod
      use partvoxMod
      use TRACK_mod
      implicit none
      integer*4 ncross
      real*8 ds,dsef

      logical ismix
      integer*4 mat0,ncquad
      real*8 x0,y0,z0,dsvox,dsquad,dsvbb,fuel
      real*8 infty,eps,oneeps,distVBB
      parameter (infty=1.0d30)
      parameter (eps=1.0d-10,oneeps=1.0d0+1.0d-12)  ! From PENGEOM

      if (absvox.eq.0) then             ! Particle in quadric geometry
        if (ibody.ne.bodymask) then     !  and not in transparent body
          call step(ds,dsef,ncross)     !  so interrogate quadrics only
          call locatevox                ! Locate in voxels
          return                        ! Stop flight
        endif

        ! Else, transparent mat but out of VBB:
        dsvbb = distVBB()*oneeps+eps    ! Distance (enlarged) to VBB
        x0 = x                          ! Save position before STEP
        y0 = y
        z0 = z
        call step(ds,dsef,ncross)       ! Distance to quadrics
        if (dsef.le.dsvbb) then         ! Quadrics are closer
          if (ncross.ne.0) call locatevox ! Vacuum gap
        else                            ! Entered voxels
          ibody = bodymask              ! So, restore transparent mat
          ncross = 1                    ! Return value
          dsef = dsvbb                  ! Ditto
          x = x0+dsvbb*u                ! Move particle up to VBB
          y = y0+dsvbb*v
          z = z0+dsvbb*w
          call locatevox                ! Set indices and mat
          if (absvox.eq.0) then         ! Consistency check, VBB missed
            call locate                 ! Reset ibody and mat
            !! write(*,'(a,3(1x,es12.5))')
            !!  'stepx:ERROR: inconsistency found at x,y,z:',x,y,z
            !! write(*,'(a)') '  Expected inside VBB, but it is not.'
            !! stop
          endif
        endif
        return                          ! Stop flight
      endif
      ! ...else particle already in voxels

      ! Get ready for the flight:
      fuel = ds                         ! 'Fuel' is distance in nominal mat density
      dsef = 0.0d0                      ! Reset travelled distance
      mat0 = mat                        ! Save original vox material
      ismix = matvox(absvox).lt.0       ! True if in an overlaid (-) vox

      ! Move particle across voxels of same mat:
      do                                ! Up to next (-) voxel or quadric
        call stepvox(fuel,dsvox,ncross) ! Find distance to (-) voxel walls
        if (ismix) then                 ! It's a (-) voxel
          x0 = x                        ! Save original position
          y0 = y                        !   (STEPVOX does not change it,
          z0 = z                        !    but STEP does)
          mat = matmask                 ! Set mat to quad geom for STEP to work
          call step(+infty,dsquad,ncquad) ! Find distance to quads
          if (dsquad.lt.dsvox+eps) then ! Entering quad geom
            absvox = 0                  ! Mark as quadrics
            dsef = dsef+dsquad          ! Accumulate travelled distance
            ncross = ncquad             ! Set ncross as found by STEP
            call locatevox              ! Re-entered after vac gap
            return                      ! xyz,ibody,mat set by STEP
          endif
          x = x0                        ! Restore original position ...
          y = y0
          z = z0
          ibody = bodymask              ! ... and body ...
          mat = mat0                    ! ... and material
        endif

        ! Voxels only, quadrics play no role:
        dsef = dsef+dsvox               ! Accumulate travelled distance
        if (ncross.eq.-9999) dsvox = dsvox*oneeps+eps ! Particle left the VBB, push out
        x = x+dsvox*u                   ! Move--no need to enlarge in voxs
        y = y+dsvox*v
        z = z+dsvox*w
        if (ncross.eq.0) return         ! Interaction inside current mat
        if (ncross.gt.0) then           ! Found voxel with new mat
          mat = ncross                  ! Update mat
          return
        endif ! ncross.lt.0             ! overlaid (-) voxel or out of VBB
        ncross = -ncross                ! Prevent returning neg values
        call locate                     ! Reset ibody and mat
        if (ibody.ne.bodymask) then     ! Entered quadrics
          absvox = 0                    ! Mark as quadrics
          if (mat.ne.0) return          ! Not in vacuum, so we're done
          call step(+infty,dsquad,ncquad) ! Skip vacuum
          call locatevox                ! Locate in voxels
          return                        ! Truncate flight after a vacuum
        endif
        if (ncross.eq.9999) return      ! In bodymask, but out of VBB
        mat = ncross                    ! Entered a (-) voxel, restore mat
        if (mat.ne.mat0) return         ! Material changed
        ismix = .true.                  ! From 2nd cycle on, overlaid vox
      enddo
      end


      subroutine stepvox(fuel,dsef,ncross)
!*******************************************************************
!*    Starting from a point inside a voxel, this routine moves     *
!*    the particle across pure (+) voxels made of the same MAT.    *
!*    The flight stops when: (i) voxel material has changed; (ii)  *
!*    an interaction happens; (iii) entered an overlaid (-) voxel  *
!*    or (iv) the particle has exited the voxels bounding box (VBB)*
!*                                                                 *
!*    Input:                                                       *
!*      fuel: max flight distance if in nominal mat density        *
!*      Voxel indices {xyz}vox,absvox                              *
!*    Output:                                                      *
!*      fuel: remaining fuel (if new mat found, else undef)        *
!*      dsef: actual distance travelled (cm)                       *
!*      ncross: 0 if interaction, new material (+ or -), or -9999  *
!*      Voxel indices {xyz}vox,absvox                              *
!*    Comments:                                                    *
!*      - Only voxel indices are updated; neither xyz nor mat are. *
!*      - If particle leaves the VBB, absvox is set to 0 and       *
!*        ncross is set to -9999.                                  *
!*******************************************************************
      use geovoxMod
      use partvoxMod
      use TRACK_mod
      implicit none
      integer*4 ncross
      real*8 fuel,dsef

      integer dxvox,dyvox,dzvox,mat0
      real*8 sx,sy,sz,dfuel,infty
      parameter (infty=1.0d30)

      ! Check if direction has changed since last call and store
      !   info to save time when going across multiple voxels:
      if (u.ne.uold.or.v.ne.vold) then
        if (u.ne.0.0d0) then
          ivx = 1.0d0/u                 ! Inverse of direction cosines
          sdx = abs(dx*ivx)             ! Distance to travel to span dx
        else
          ivx = +infty
          sdx = +infty
        endif
        if (v.ne.0.0d0) then
          ivy = 1.0d0/v
          sdy = abs(dy*ivy)
        else
          ivy = +infty
          sdy = +infty
        endif
        if (w.ne.0.0d0) then
          ivz = 1.0d0/w
          sdz = abs(dz*ivz)
        else
          ivz = +infty
          sdz = +infty
        endif
        uold = u                        ! Save new particle's direction
        vold = v
      endif

      ! Determine min distance to the three relevant voxel walls:
      if (ivx.gt.0.0d0) then
        sx = (xvox*dx-x)*ivx
        dxvox = +1                          ! Will move one voxel forward
      else
        sx = ((xvox-1)*dx-x)*ivx
        dxvox = -1                          ! Will move one voxel backwards
      endif
      if (ivy.gt.0.0d0) then
        sy = (yvox*dy-y)*ivy
        dyvox = +1
      else
        sy = ((yvox-1)*dy-y)*ivy
        dyvox = -1
      endif
      if (ivz.gt.0.0d0) then
        sz = (zvox*dz-z)*ivz
        dzvox = +1
      else
        sz = ((zvox-1)*dz-z)*ivz
        dzvox = -1
      endif

      dsef = 0.0d0                          ! Reset traveled distance
      mat0 = mat                            ! Store original material
      do                                    ! One voxel per cycle
        if (sx.lt.sy) then
          if (sx.lt.sz) then                ! Intersection with x wall

            sx = max(sx,0.0d0)              ! Never go backwards
            dfuel = sx*densvox(absvox)      ! Dist2Surf*dens(vox)/dens(nominal)
            if (dfuel.gt.fuel) then         ! Run out of fuel => interaction
              ncross = 0                    ! Signals an interaction
              dsef = dsef+fuel/densvox(absvox) ! Distance up to interaction
              return                        ! Stop flight
            endif
            dsef = dsef+sx                  ! Update travelled distance
            xvox = xvox+dxvox               ! Update coordinate voxel index
            if (xvox.lt.1.or.xvox.gt.nx) then  ! Particle is gone
              ncross = -9999
              absvox = 0                    ! Mark as quadrics
              return
            endif
            fuel = fuel-dfuel               ! Remove spent fuel from tank
            absvox = absvox+dxvox           ! Update abs voxel index
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return      ! MAT changed or (-) voxel
            sy = sy-sx                      ! Update distances to walls
            sz = sz-sx
            sx = sdx                        ! Reset to full voxel width
            cycle                           ! Process next voxel

          else                              ! Intersection with z wall
                                            ! Repeat code above for z
            sz = max(sz,0.0d0)
            dfuel = sz*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel/densvox(absvox)
              return
            endif
            dsef = dsef+sz
            zvox = zvox+dzvox
            if (zvox.lt.1.or.zvox.gt.nz) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dzvox*nxy
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sx = sx-sz
            sy = sy-sz
            sz = sdz
            cycle

          endif
        else
          if (sy.lt.sz) then                ! Intersection with y wall
                                            ! Repeat code above for y
            sy = max(sy,0.0d0)
            dfuel = sy*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel/densvox(absvox)
              return
            endif
            dsef = dsef+sy
            yvox = yvox+dyvox
            if (yvox.lt.1.or.yvox.gt.ny) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dyvox*nx
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sz = sz-sy
            sx = sx-sy
            sy = sdy
            cycle

          else                              ! Intersection with z wall
                                            ! Repeat code above for z
            sz = max(sz,0.0d0)
            dfuel = sz*densvox(absvox)
            if (dfuel.gt.fuel) then
              ncross = 0
              dsef = dsef+fuel/densvox(absvox)
              return
            endif
            dsef = dsef+sz
            zvox = zvox+dzvox
            if (zvox.lt.1.or.zvox.gt.nz) then
              ncross = -9999
              absvox = 0
              return
            endif
            fuel = fuel-dfuel
            absvox = absvox+dzvox*nxy
            ncross = matvox(absvox)
            if (ncross.ne.mat0) return
            sx = sx-sz
            sy = sy-sz
            sz = sdz
            cycle

          endif
        endif
      enddo
      end


      real*8 function distVBB()
!*******************************************************************
!*    Returns distance to intersection between the particle's      *
!*    trajectory and the VBB. If there is no intersection, +infty  *
!*    is returned.                                                 *
!*                                                                 *
!*    Input:                                                       *
!*      - xyz, uvw, vbb needs to be defined.                       *
!*    Comments:                                                    *
!*      - Max considered distance is infty=1e90 cm.                *
!*******************************************************************
      use geovoxMod
      use TRACK_mod
      implicit none

      integer i
      real*8 r(3),d(3),infty,slo,shi,invd
      parameter (infty=1.0d90)

      r = (/ x,y,z /)    ! Starting point
      d = (/ u,v,w /)    ! Direction
      slo = 0.0d0        ! Position 'low pointer' at start of segment
      shi = +infty       ! Position 'high pointer' at end of segment
      distVBB = +infty   ! Infinity signals 'no intersection'

      do i=1,3                         ! For x,y,z
        if (d(i).eq.0.0d0) then        ! Going parallel to surface
          if (r(i).lt.0.0d0.or.r(i).gt.vbb(i)) return ! No VBB inters.
          cycle                        ! Otherwise, no need to clip
        endif
        invd = -1.0d0/d(i)             ! Avoids repeating the division
        if (invd.lt.0.0d0) then        ! Going forward
          if (r(i).gt.vbb(i)) return ! No intersection
          shi = min(shi,(r(i)-vbb(i))*invd)         ! Clip high end
          if (r(i).lt.0.0d0) slo = max(slo,r(i)*invd) ! Clip low end
        else                           ! Going backwards
          if (r(i).lt.0.0d0) return    ! No intersection
          shi = min(shi,r(i)*invd)     ! Clip high end
          if (r(i).gt.vbb(i)) slo = max(slo,(r(i)-vbb(i))*invd) ! Low end
        endif
        if (shi.lt.slo) return         ! No intersection
      enddo
      distVBB = slo                    ! Distance to VBB
      end


      subroutine locatevox
!*******************************************************************
!*    Determines the voxel that contains the point with coordinates*
!*    {x,y,z} (taken from TRACK_MOD) and sets the voxel indices and*
!*    the material. If the position is not in a voxel, the absolute*
!*    voxel index ABSVOX is set to 0---taking into account the     *
!*    consideration of being overlaid with a quadric.              *
!*******************************************************************
      use geovoxMod
      use partvoxMod
      use TRACK_mod
      implicit none

      if (ibody.ne.bodymask) then         ! Not in quadric transparent body, i.e., not in a voxel
        absvox = 0                        ! Mark as quadrics
        return                            ! Done
      endif

      xvox = x*idx+1.0d0
      if (xvox.lt.1.or.xvox.gt.nx) then   ! Outside voxel's bounding box (VBB)
        absvox = 0                        ! Set mark to quadrics
        return                            ! Done
      endif
      yvox = y*idy+1.0d0
      if (yvox.lt.1.or.yvox.gt.ny) then
        absvox = 0
        return
      endif
      zvox = z*idz+1.0d0
      if (zvox.lt.1.or.zvox.gt.nz) then
        absvox = 0
        return
      endif
      absvox = xvox+(yvox-1)*nx+(zvox-1)*nxy  ! Abs voxel index
      mat = abs(matvox(absvox))
      end


      subroutine inivox(nmatvox,mem)
!*******************************************************************
!*    Initializes the voxel geometry package, setting material and *
!*    density arrays.                                              *
!*                                                                 *
!*    Output:                                                      *
!*      nmatvox: no. of materials in VOX file.                     *
!*      mem: memory used by voxel arrays.                          *
!*    Comments:                                                    *
!*      - Computation of voxels mass and material sign is not      *
!*        performed here since PENELOPE has not been initialized,  *
!*        and so the mass density for each material has yet to be  *
!*        defined. It will be determined later.                    *
!*******************************************************************
      use geovoxMod
      use partvoxMod, only : uold,vold,ivx,ivy,ivz
      use PENGEOM_mod, only : mater,nbody
      implicit none
      integer*4 nmatvox
      real*8 mem,nowcpu,cputime

      character*80 buffer
      integer nmask
      integer*4 i

      ! Read voxels file, if there is one:
      read(*,'(a80)') voxfilen ! Quadrics filename
      voxfilen = adjustl(voxfilen)
      voxfilen = voxfilen(1:scan(voxfilen,' ')) ! Clip at 1st blank
      write(*,*) ''
      if (voxfilen.eq.'-') then
        write(*,'(a)') 'No voxels geometry file defined.'
        read(*,'(a80)') buffer   ! Dummy (unused) lines from input
        read(*,'(a80)') buffer
        nmatvox = 0              ! No voxel materials
        bodymask = -1            ! Impossible transparent body signals no voxelized geo
        return                   ! Nothing else to do
      endif                      ! Else, a vox file exists
      write(*,'(a)')
     &  '>>>> Initializing voxelized geometry >>>>'
      write(*,'(a)')
     &  'NOTICE: when voxels are active, some tallies may not yield'
      write(*,'(a)')
     &  '  correct values if evaluated inside the voxelized region'
      write(*,'(a)')
     &  '  (e.g. fluence, cylindrical and spherical dose distrib).'
      write(*,'(a)')
     &  '  Read supplied documentation for more info.'

      nowcpu = cputime()         ! Start stopwatch
      call readvox(nmatvox,mem)

      if (isQuad) then           ! Quadrics geom present
        write(*,'(a)') 'Quadric material transparent to voxels:'
        read(*,*) matmask
        write(*,'(1x,i0)') matmask
        if (matmask.lt.1) then
          write(*,'(a)') 'inivox:ERROR: invalid transparent mat.'
          stop  ! Vacuum forbidden as transparent mat (but allowed elsewhere)
        endif   !   because it would be missed by STEP
        write(*,'(a)')
     &    'Searching for bodies containing the transparent material...'
        nmask = 0                ! No. of bodies with matmask
        do i=1,nbody
          if (mater(i).eq.matmask) then  ! This body contains the transparent mat
            bodymask = i
            nmask = nmask+1
            write(*,'(a,i0)')
     &        '  found in body no. (PENGEOM internal numbering): ',i
          endif
        enddo
        if (nmask.ne.1) then
          write(*,'(a)') 'inivox:ERROR: '//
     &      'there must be ONE body with transparent material.'
          stop
        endif

        write(*,'(a)') 'Voxel scan granularity:'
        read(*,*) granul  ! No. rays along voxel side (granul**2 per voxel)
        write(*,'(1x,i0)') granul
        if (granul.lt.2.or.granul.gt.maxGranul) then
          write(*,'(a,i0)')
     &      'inivox:ERROR: granularity must be between 2 and ',
     &      maxGranul
          stop
        endif

      else                       ! No quadric geometry
        write(*,'(a)') 'No quadrics present, voxels assumed in vacuum.'
        read(*,*) matmask        ! Read a dummy line
        matmask = 0              ! Unused, but just in case
        if (mater(1).eq.0) then  ! Only one body is defined in default GEO
          bodymask = 1
        else
          write(*,'(a,i0)')
     &      'inivox:ERROR: internal inconsistency; expecting '//
     &      'vacuum, found instead: ',mater(1)
          stop
        endif

        read(*,*) granul         ! Read a dummy line
        granul = 2               ! Unused, but just in case
      endif

      ! Initialize particle state variables with dummy values:
      uold = 1.0d0
      vold = 1.0d0
      ivx  = 1.0d0
      ivy  = 1.0d0
      ivz  = 1.0d0

      ! Note that voxels mass and sign are not computed here;
      !   see header for comments on this.

      ! Report time usage:
      write(*,'(a)') 'CPU time used for voxels init (s):'
      write(*,'(1x,f9.2)') cputime()-nowcpu

      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Voxelized geometry initialization finished >>>>'
      end


      subroutine readvox(nmatvox,mem)
!*******************************************************************
!*    Reads the voxel geometry file and allocates arrays.          *
!*                                                                 *
!*    Output:                                                      *
!*      nmatvox: no. of materials in VOX file.                     *
!*      mem: memory used by voxel arrays to load mat & dens.       *
!*******************************************************************
      use geovoxMod
      use dataTypesMod
      implicit none
      integer*4 nmatvox
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION VOXELS HEADER v.2008-04-13]')
      parameter (eos='[END OF VXH SECTION]')
      character*500 buffer
      integer errorcode,finduf,ufile,mat,isblank
      integer colMat,colDens,maxCol,ncol,n
      parameter (maxCol=20)
      integer*4 i,j,k,line,vox
      integer*4, parameter :: nvoxmax=int(2.14748e9,kind=kindInt4)  ! ~max int*4
      real zero,dens,col(maxCol)
      real*8 minvoxSide,eps
      parameter (eps=1.0d-10,minvoxSide=eps*1.0d4,zero=1.0e-30)

      write(*,*) ''
      write(*,'(a)') '>>>> Reading voxels file >>>>'
      write(*,'(a)') 'Opening file:'
      write(*,'(1x,a30)') voxfilen
      ufile = finduf()  ! Find a valid unit file
      open(ufile,file=voxfilen,status='old',iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'readvox:ERROR: unable to open voxels file'
        stop
      endif

      line = 0                           ! Clear file line number
      call getvoxline                    ! Get a non-comment line into buffer
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'readvox:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)') 'No. of voxels in x,y,z and total:'
      line = line+1                     ! Update line counter
      read(ufile,*) nx,ny,nz            ! Note that no comments are allowed inside section
      if (min(nx,ny,nz).lt.1) then
        write(*,'(a)') 'readvox:ERROR: invalid no. voxels.'
        stop
      endif
      nxy = nx*ny
      nvox = nxy*nz
      write(*,'(4(1x,i0))') nx,ny,nz,nvox
      if (dble(nx)*dble(ny)*dble(nz).gt.dble(nvoxmax)) then  ! Check overflow in double prec.
        write(*,'(a,i0)') 'readvox:ERROR: No. voxels exceeds ',
     &    nvoxmax
        stop
      endif
      ! Allocate arrays:
      mem = (sizeofInt+sizeofReal)*real(nvox)
      write(*,'(a)') 'Memory required to allocate voxel arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6
      allocate(matvox(nvox),densvox(nvox),stat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'readvox:ERROR: not enough memory.'
        stop
      endif
      memvox = memvox+mem  ! Store for later use by voxels dose report

      write(*,'(a)') 'Voxel dimensions in x,y,z (cm):'
      line = line+1
      read(ufile,*) dx,dy,dz
      write(*,'(3(1x,es12.5))') dx,dy,dz
      if (min(dx,dy,dz).lt.minvoxSide) then
        write(*,'(a)')
     &    'readvox:ERROR: voxel side too small, tracking algorithm'
        write(*,'(a,es12.5)')
     &    '  requires voxel sides to be larger than (cm):',minvoxSide
        stop
      endif
      volvox  = dx*dy*dz
      write(*,'(a)') 'Voxels volume (cm^3):'
      write(*,'(1x,es12.5)') volvox
      idx = 1.0d0/dx
      idy = 1.0d0/dy
      idz = 1.0d0/dz
      vbb = (/dx*nx,dy*ny,dz*nz/)
      write(*,'(a)') 'Size of voxels bounding box, VBB (cm):'
      write(*,'(3(1x,es12.5))') (vbb(i), i=1,3)
      if (maxval(vbb).gt.1.0d5) then
        write(*,'(a)') 'readvox:ERROR: VBB too large.'
        stop
      endif

      write(*,'(a)') 'Material ID expected in column:'
      line = line+1
      read(ufile,*) colMat
      write(*,'(1x,i0)') colMat
      write(*,'(a)') 'Mass density expected in column:'
      line = line+1
      read(ufile,*) colDens
      write(*,'(1x,i0)') colDens
      ncol = max(colMat,colDens)
      if (min(colMat,colDens).lt.1.or.ncol.gt.maxCol)
     &  then
        write(*,'(a,i0)')
     &    'readvox:ERROR: column numbers must be between 1 and ',
     &     maxCol
        stop
      endif

      line = line+1
      read(ufile,*) isblank
      if (isblank.eq.1) then
        write(*,'(a)')
     &    'Expecting blank lines to separate y,z-cycles'
      else
        write(*,'(a)')
     &    'NOT expecting blank lines to separate y,z-cycles'
      endif

      ! Check section integrity:
      line = line+1
      read(ufile,'(a500)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'readvox:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)') 'Now loading mat and dens for each voxel...'
      nmatvox = 0                             ! Init num. mats in vox geometry
      do k=1,nz
        do j=1,ny

          do i=1,nx
            call getvoxline
            read(buffer,*,iostat=errorcode) (col(n), n=1,ncol)
            mat  = int(col(colMat)+0.1)  ! Round to integer
            dens = col(colDens)
            if (errorcode.ne.0.or.mat.lt.1.or.dens.lt.0.0) then
              write(*,'(a,i0)')
     &          'readvox:ERROR: invalid entry at line ',line
              write(*,'(a)') '  Line contents:'
              write(*,'(a)') buffer
              stop
            endif
            vox = i+(j-1)*nx+(k-1)*nxy        ! Absolute voxel index
            matvox(vox)  = mat
            densvox(vox) = max(zero,dens)     ! Min dens is +zero
            ! Note: densvox will be redefined later to normalize
            !   to the nominal mass density, not known here before PENELOPE init
          enddo

          if (isblank.eq.1.and.(j.lt.ny.or.k.lt.nz)) then
            call getvoxline                   ! Read one blank line except in the last line
            if (buffer.ne.' ') then
              write(*,'(a,i0)')
     &          'readvox:ERROR: Line should be blank, line ',line
              stop
            endif
          endif
        enddo

        if (isblank.eq.1.and.k.lt.nz) then     ! Read one blank line except in the last line
          call getvoxline
          if (buffer.ne.' ') then
            write(*,'(a,i0)')
     &        'readvox:ERROR: Line should be blank, line ',line
            stop
          endif
        endif

        ! Report progress every 10%:
        if (mod(10*k,nz).lt.10)
     &    write(*,'(a,i0,a,$)') ',',10*((10*k)/nz),'%'
      enddo
      write(*,*) ''                            ! Terminate progress report line
      close(ufile)

      nmatvox = maxval(matvox)                 ! Largest mat index in file
      write(*,'(a,i0,a)') 'Done. ',line,' lines processed.'

      contains
        subroutine getvoxline
        !*******************************************************************
        !*    Reads a new line from voxels file. The line is returned      *
        !*    only if it is not a comment line.                            *
        !*                                                                 *
        !*    Input (from parent routine):                                 *
        !*      ufile: file unit to be read                                *
        !*      line: last line number read from that file                 *
        !*    Output (from parent routine):                                *
        !*      line: updated line number                                  *
        !*      buffer: line read from file                                *
        !*    Comments:                                                    *
        !*      - comment lines start with a '#'                           *
        !*      - blank lines are NOT comments                             *
        !*******************************************************************
        implicit none
        do
          line = line+1
          read(ufile,'(a)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a)')
     &        'getvoxline:ERROR: unable to read vox file at line:'
            write(*,'(1x,i0)') line
            stop
          endif
          if (buffer(1:1).ne.'#') exit  ! Else a comment line, read another
        enddo
        end subroutine
      end


      subroutine inimassvox(mem)
!*******************************************************************
!*    Performs three main tasks:                                   *
!*    1. Compute the mass of all voxels by integrating density over*
!*       volume. In case of voxels overlaid with quadrics, voxel   *
!*       masses exclude the portion covered by quadrics.           *
!*    2. Change the sign of a voxel material to negative if the    *
!*       voxel is overlaid with a quadric.                         *
!*    3. Renormalize voxel densities by dividing by their          *
!*       materials' nominal density.                               *
!*                                                                 *
!*    Output:                                                      *
!*      mem: memory used by permanent arrays in this routine.      *
!*    Comments:                                                    *
!*      - Must be called after PENELOPE init because requires      *
!*        knowledge of the materials density (PEINIT).             *
!*******************************************************************
      use geovoxMod
      use dataTypesMod
      use TRACK_mod
      use COMPOScommonMod, only : rho
      implicit none
      real*8 mem

      integer errorcode
      integer*4 xvox,yvox,zvox,i,j,vox,nOverlay
      integer*4 voxIni,zvoxIni,zvoxFin,ncross,matold
      real*8 delta,eps,oneeps,nowcpu,cputime
      parameter (delta=0.01d0)                      ! A small fraction of unity
      parameter (eps=1.0d-10,oneeps=1.0d0+1.0d-12)  ! From PENGEOM
      real*8 darea,zinf,xIni,yIni,zIni,dist,dsef,shift
      real*8, allocatable :: da(:,:),xg(:),yg(:)

      if (bodymask.eq.-1.or..not.isQuad) then ! No voxelized geometry or no quadrics
        isOverlay = .false.                   ! There are no overlaid voxels
        do vox=1,nvox                         ! Renormalize voxels density to nominal mat density
          densvox(vox) = densvox(vox)/rho(matvox(vox))
        enddo
        return                                ! Nothing else to do, mass is simply volvox*densvox
      endif

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') '[QUADRICS&VOXELS OVERLAYS PROCESSING]'
      write(*,'(a)') 'Memory required to allocate temp arrays (MB):'
      write(*,'(1x,es12.5)')
     &  sizeofReal8*real(nvox+granul**2+granul*2)*1.0e-6
      allocate(massvoxd(nvox),da(granul,granul),xg(granul),yg(granul),
     &  stat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'inimassvox:ERROR: '//
     &    'not enough memory for temp arrays.'
        stop
      endif

      write(*,'(a)') 'Evaluating overlays and voxels mass...'
      nowcpu = cputime()                      ! Start stopwatch

      ! Set grid points and weights to sweep a voxel according to GRANUL:
      da = dx*dy/(granul-1)**2                ! Area of cells in which voxels area subdivided
      da(:,1)      = 0.5*da(:,1)
      da(:,granul) = 0.5*da(:,granul)         ! Halve weights of 1st and last columns
      da(1,:)      = 0.5*da(1,:)
      da(granul,:) = 0.5*da(granul,:)         ! Halve weights of 1st and last rows; note that corners were halved twice
      do i=1,granul
        xg(i) = (i-1)*dx/(granul-1)
        yg(i) = (i-1)*dy/(granul-1)
      enddo

      ! Shift sides:
      shift = delta*dx/(granul-1)
      xg(1) = shift
      xg(granul) = dx-shift
      shift = delta*dy/(granul-1)
      yg(1) = shift
      yg(granul) = dy-shift
      ! ...Note that corners are shifted both in x and y

      ! Compute voxels mass and sign:
      zinf = 1.01*nz*dz                       ! Practical infinity (>VBB in z) avoids overflows
      u = 0.0
      v = 0.0
      w = 1.0                                 ! Pointing upwards
      massvoxd = 0.0                          ! Init mass integrator
      do yvox=1,ny                            ! For each voxel
        yIni = (yvox-1)*dy                    ! Coordinates of the voxel corner
        do xvox=1,nx
          xIni = (xvox-1)*dx
          voxIni = xvox+(yvox-1)*nx           ! Absolute voxel index for this column

          do j=1,granul                       ! For each ray (i.e. each cell) in a voxel
            do i=1,granul
              darea = da(i,j)                 ! dArea to compute voxel mass
              x = xIni+xg(i)                  ! Coordinates of the {x,y} ray
              y = yIni+yg(j)
              z = 0.0                         ! Initial z on the VBB
              call locate                     ! Set ibody and mat
              zvoxIni = 1                     ! Initial z index for this ray
              vox = voxIni                    ! Init voxel index for this ray
              zIni = z                        ! First raystep position

              raystep: do                     ! For each step in a ray
                matold = mat                  ! Save current mat to compute vox mass
                z = zIni                      ! Set start of next step for STEP to work
                call step(zinf,dsef,ncross)   ! Compute next step of current ray
                if (mat.eq.0) then            ! Gone, last step
                  if (matold.eq.0) then       ! From vac to vac, no matter
                    z = +zinf                 ! Set z to practical infinity
                  else
                    z = zIni+oneeps*dsef+eps  ! Reposition in close vacuum
                  endif
                else if (ncross.gt.1) then    ! Matter after vaccum gap
                  z = zIni+oneeps*dsef+eps    ! Reposition inside vac gap
                endif
                zvoxFin = z*idz+1.0           ! Final voxel of the current step

                do zvox=max(1,zvoxIni),min(nz,zvoxFin)           ! For each visited voxel along z
                  ! Compute distance travelled inside this voxel:
                  dist = dz                                      ! For a whole voxel
                  if (zvox.eq.zvoxIni) dist = dist-zIni+(zvox-1)*dz  ! Correct start-of-step
                  if (zvox.eq.zvoxFin) dist = dist+z-zvox*dz     ! Correct end-of-step
                  ! Compute mass and mark overlaid voxels:
                  if (matold.eq.matmask) then                    ! Moved inside voxels
                    massvoxd(vox) = massvoxd(vox)+
     &                darea*dist*densvox(vox)                    ! Add contribution to voxel mass
                  else                                           ! Moved in quadric part of an overlaid voxel
                    matvox(vox) = -abs(matvox(vox))              ! Mark voxel as overlaid
                  endif
                  vox = vox+nxy                                  ! Move one voxel up along z
                enddo
                vox = vox-nxy                                    ! Move 'pointer' back to last voxel

                zIni = z                         ! Update next step initial position
                zvoxIni = zvoxFin                ! Update next step initial index
                if (zvoxIni.gt.nz) exit raystep  ! All voxels along this ray have been completed
              enddo raystep                      ! Cycle next step
            enddo
          enddo
        enddo

        ! Report progress every 10%:
        if (mod(10*yvox,ny).lt.10)
     &    write(*,'(a,i0,a,$)') ',',10*((10*yvox)/ny),'%'
      enddo
      write(*,*) ''                            ! Terminate progress report line
      deallocate(da,xg,yg)

      ! Renormalize voxels density to nominal mat density:
      do vox=1,nvox
        densvox(vox) = densvox(vox)/rho(abs(matvox(vox)))  ! ABS() because overlaid voxels have neg material
      enddo

      ! Evaluate overlays:
      nOverlay = count(matvox.lt.0,kind=kindInt4)
      if (nOverlay.gt.0) then                    ! Overlaying quadrics means voxels mass can be variable
        isOverlay = .true.                       ! Set flag
        write(*,'(a)') 'Fraction of voxels overlaid with quadrics (%):'
        write(*,'(1x,es12.5)') nOverlay/real(nvox)*100.0
        mem = sizeofReal*real(nvox)
        write(*,'(a)')
     &    'Done. Memory required to allocate overlay arrays (MB):'
        write(*,'(1x,es12.5)') mem*1.0e-6
        allocate(massvox(nvox),stat=errorcode)   ! Array to allocate voxel mass in 4 bytes
        if (errorcode.ne.0) then
          write(*,'(a)') 'inimassvox:ERROR: not enough memory.'
          stop
        endif
        massvox = massvoxd                       ! Transfer data from real8 to real4
      else
        isOverlay = .false.                      ! Set flag
        mem = 0.0
        write(*,'(a)') 'Done. None of the voxels is overlaid '//
     &    'with quadrics; no permanent memory allocation required.'
      endif
      deallocate(massvoxd)                       ! Free real8 memory
      ! Memory usage:
      memvox = memvox+mem

      ! Report time usage:
      write(*,'(a)') 'CPU time used for overlay processing (s):'
      write(*,'(1x,f9.2)') cputime()-nowcpu

      !! call writeMassvox                       ! Write voxels mass and sign to a file (for future versions)
      write(*,'(a)') '>>>> Voxel overlays processing finished >>>>'
      end


      real*8 function fmassvox(vox)
!*******************************************************************
!*    Returns the mass of a voxel.                                 *
!*                                                                 *
!*    Input:                                                       *
!*      vox: absolute voxel index.                                 *
!*******************************************************************
      use geovoxMod, only : isOverlay,densvox,matvox,volvox,massvox
      use COMPOScommonMod, only : rho
      implicit none
      integer*4 vox

      if (isOverlay) then         ! Quadrics overlay voxels
        fmassvox = dble(massvox(vox))
      else                        ! All voxels are pure and matvox>0
        fmassvox = volvox*densvox(vox)*rho(matvox(vox))
      endif
      end


      real*8 function flocaldens()
!*******************************************************************
!*    Returns the local mass density, taking care of whether the   *
!*    point lies in a quadric or in a voxel.                       *
!*******************************************************************
      use geovoxMod, only : densvox,matvox
      use partvoxMod, only : absvox
      use TRACK_mod, only : mat
      use COMPOScommonMod, only : rho
      implicit none

      call locate            ! Sets mat
      call locatevox         ! Sets absvox
      if (absvox.eq.0) then  ! Point inside quadrics
        if (mat.gt.0) then
          flocaldens = rho(mat)
        else
          flocaldens = 0.0
        endif
      else                   ! Point in a voxel
        flocaldens = densvox(absvox)*rho(abs(matvox(absvox)))  ! Overlaid voxels have matvox<0
      endif
      end


      subroutine writeMassvox
!*******************************************************************
!*    Writes voxels mass and sign to an external file for reference*
!*                                                                 *
!*    Comments:                                                    *
!*    -> Currently unused, awaiting for future versions that would *
!*       allow to load the mass file instead of recalculating it.  *
!*******************************************************************
      use geovoxMod
      implicit none
      character*90 buffer
      integer ufile,finduf,errorcode
      integer*4 vox,i,j,k,dk,djk
      real*8 fmassvox

      buffer = trim(voxfilen)//'.mass'
      ufile = finduf()  ! Find a valid unit file
      open(ufile,file=buffer,iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a,a)') 'writeMassvox:ERROR: unable to open file: ',
     &    buffer
        stop
      endif

      write(ufile,'(a,a)')
     &  '# This file was automatically created by penEasy from '//
     &  'the voxels file named: ',trim(voxfilen)
      write(ufile,'(a)')
     &  '# It contains the list of voxels mass (in g) with sign,'
      write(ufile,'(a)')
     &  '#   (+) for pure voxels and (-) for voxels overlaid '//
     &  'by a quadric body.'
      write(ufile,'(a)')
     &  '# The voxels mass excludes the volume covered by quadrics.'
      write(ufile,'(a,i0)')
     &  '# The granularity used to compute masses was: ',granul
      write(ufile,'(a)') '#'

      ! Write data to file:
      do k=1,nz
        write(ufile,'(a,i0)') '# zVoxIndex=',k
        dk = (k-1)*nxy
        do j=1,ny
          write(ufile,'(a,i0)') '# yVoxIndex=',j
          djk = (j-1)*nx+dk
          do i=1,nx
            vox = i+djk
            write(ufile,'(es12.5)') fmassvox(vox)*sign(1,matvox(vox))
          enddo
          if (nx.gt.1) write(ufile,*) '' ! Separate 2D data blocks
        enddo
        if (ny.gt.1) write(ufile,*) ''   ! Separate 3D data blocks
      enddo
      close(ufile)
      end


      subroutine autogeo(filename)
!*******************************************************************
!*    Creates a PENGEOM file with an empty (filled with vacuum)    *
!*    geometry. Used when no GEO file has been provided.           *
!*                                                                 *
!*    Input:                                                       *
!*      filename: name of the file to be created                   *
!*******************************************************************
      implicit none
      character*(*) filename

      character*80 linezero,lineend,linesurf,lineindex,linescale
      character*80 lineshift,linebody,linemat,linepoint
      parameter (linezero='0000000000000000000000000000000000000'//
     &                    '000000000000000000000000000')
      parameter (lineend ='END      0000000000000000000000000000'//
     &                    '000000000000000000000000000')
      parameter (linesurf  ='SURFACE (')
      parameter (lineindex ='INDICES=(')
      parameter (linescale ='-SCALE=(+')
      parameter (lineshift ='-SHIFT=(+')
      parameter (linebody  ='BODY    (')
      parameter (linemat   ='MATERIAL(')
      parameter (linepoint =', SIDE POINTER=(')
      integer ufile,finduf

      ufile = finduf()
      open(ufile,file=filename)
      write(ufile,'(a)')
     & '# This PENGEOM file was automatically created by penEasy.'
      write(ufile,'(a)')
     & '# It can be removed only after the execution is completed.'
      write(ufile,'(a)')
     & '# Do not edit this file while the code is running.'
      write(ufile,*) ''
      write(ufile,'(a80)') linezero         ! Start with a line of zeroes
      write(ufile,'(a9,i4,")")') linebody,1 ! Define a body
      write(ufile,'(a9,i4,")")') linemat,0  ! Filled with vacuum (no surfs)
      write(ufile,'(a80)') linezero         ! Another line of zeroes
      write(ufile,'(a80)') lineend          ! And the end line
      close(ufile)

      ! Other possible elements:
      !!! Pair of X planes:
      !!write(ufile,'(a9,i4,")")') linesurf,1
      !!write(ufile,'(a9,4(i2,","),i2,")")') lineindex,1,0,0,0,-1
      !!write(ufile,'(a1,a9,es21.15,",",i4,")")') 'X',
     &!!  linescale,xscale,0
      !!write(ufile,'(a1,a9,es21.15,",",i4,")")') 'X',
     &!!  lineshift,xscale,0
      !!write(ufile,'(a80)') linezero
      !!! Limiting quadrics:
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,1,linepoint,-1
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,2,linepoint,-1
      !!write(ufile,'(a9,i4,")",a16,i2,")")') linesurf,3,linepoint,-1
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
