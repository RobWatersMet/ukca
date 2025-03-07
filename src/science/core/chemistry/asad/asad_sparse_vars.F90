module asad_sparse_vars
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_SPARSE_VARS'
  public

  ! Overloaded interface for both integer and real arrays
  interface array_xor_checksum
    module procedure array_xor_checksum_int
    module procedure array_xor_checksum_real
    module procedure array_xor_checksum_logical
  end interface

  contains

  function array_xor_checksum_int(arr) result(checksum)
    implicit none
    integer, intent(in)  :: arr(..)
    integer              :: checksum
    integer              :: i

    checksum = 0
    do i = 1, size(arr)
      checksum = ieor(checksum, arr(i))
    end do
  end function array_xor_checksum_int

  function array_xor_checksum_real(arr) result(checksum)
    implicit none
    real, intent(in)     :: arr(..)
    integer              :: checksum
    integer              :: i, int_value

    checksum = 0
    do i = 1, size(arr)
        int_value = int(arr(i) * 1e6)  ! Convert to int with scaling
        checksum = ieor(checksum, int_value)
    end do
  end function array_xor_checksum_real

  function array_xor_checksum_logical(arr) result(checksum)
    implicit none
    logical, intent(in)     :: arr(..)
    integer                 :: checksum
    integer                 :: i, int_value

    checksum = 0
    do i = 1, size(arr)
        int_value = merge(1, 0, arr(i))
        checksum = ieor(checksum, int_value)
    end do
  end function array_xor_checksum_logical

  subroutine setup_spfuljac()
    use asad_mod, only : base_tracer, ffrac, fracterms, frpx, jpcspf, jpfrpx, jpmsp, jpspec, madvtr, maxfterms, maxterms, modified_map, ndepd, ndepw, negterms, nfracterms, nfrpx, njcoth, nltrf, nmsjac, nmzjac, nnegterms, nonzero_map, nonzero_map_unordered, npdfr, nposterms, nsjac1, nstst, ntabpd, ntrf, ntro3, nzjac1, posterms, reorder, specf, spfjsize_max, total, ztabpd
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_normal, prstatus_oper, ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use ukca_um_legacy_mod, only : mype
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: total1
    integer :: errcode
    integer :: irj
    integer :: itrd
    integer :: i
    integer :: j
    integer :: jc
    integer :: itrcr
    integer :: j3
    integer :: jn
    integer :: js
    integer :: i1
    integer :: i2
    integer :: isx
    integer :: kr
    integer :: ikr
    integer :: krj
    integer :: ij1
    integer :: itemp1
    integer, dimension(jpmsp) :: ij
    integer, dimension(jpcspf) :: activity
    integer, allocatable, dimension(:,:) :: permute
    integer, allocatable, dimension(:,:) :: map
    logical, dimension(jpmsp) :: lj
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'SETUP_SPFULJAC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19

    #ifdef ARRAYCHECKSUM
      integer :: checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (.NOT.ALLOCATED(map)) then
      ALLOCATE(map(1:jpcspf,1:jpcspf))
    end if

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(map)
      print*, 'Checksum - map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(ij)
      print*, 'Checksum - ij in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(lj)
      print*, 'Checksum - lj in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nonzero_map)
      print*, 'Checksum - nonzero_map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(base_tracer)
      print*, 'Checksum - base_tracer in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(activity)
      print*, 'Checksum - activity in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(reorder)
      print*, 'Checksum - reorder in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(permute)
      print*, 'Checksum - permute in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(modified_map)
      print*, 'Checksum - modified_map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nposterms)
      print*, 'Checksum - nposterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nnegterms)
      print*, 'Checksum - nnegterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nfracterms)
      print*, 'Checksum - nfracterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(posterms)
      print*, 'Checksum - posterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(negterms)
      print*, 'Checksum - negterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(fracterms)
      print*, 'Checksum - fracterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(ffrac)
      print*, 'Checksum - ffrac in ', RoutineName, ' = ', checksum_result
    #endif

    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(map, dim=2), UBOUND(map, dim=2), 1
      do idx_1 = LBOUND(map, dim=1), UBOUND(map, dim=1), 1
        map(idx_1,idx) = 0
      enddo
    enddo
    !$acc loop gang vector independent
    do i = 1, jpcspf, 1
      map(i,i) = 1
    enddo
    !$acc loop seq
    do jc = 1, ntrf, 1
      itrcr = nltrf(jc)
      !$acc loop seq
      do j3 = 1, nmzjac(itrcr), 1
        irj = nzjac1(j3,itrcr)
        !$acc loop gang vector independent
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        !$acc loop seq
        do jn = 1, jpmsp, 1
          if (lj(jn)) then
            map(ij(jn),itrcr) = 1
          end if
        enddo
        if (npdfr(irj,1) /= 0) then
          i1 = npdfr(irj,1)
          i2 = npdfr(irj,2)
          !$acc loop seq
          do jn = i1, i2, 1
            isx = ntabpd(jn,1)
            map(isx,itrcr) = 1
          enddo
        end if
      enddo
    enddo
    !$acc loop seq
    do jc = 1, nstst, 1
      !$acc loop seq
      do j3 = 1, nmsjac(jc), 1
        irj = nsjac1(j3,jc)
        !$acc loop gang vector independent
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        !$acc loop seq
        do jn = 1, jpmsp, 1
          if (lj(jn)) then
            map(ij(jn),ntro3) = 1
          end if
        enddo
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(map)
      print*, 'Checksum - map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(ij)
      print*, 'Checksum - ij in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(lj)
      print*, 'Checksum - lj in ', RoutineName, ' = ', checksum_result
    #endif

    if (ndepw /= 0 .OR. ndepd /= 0) then
      !$acc kernels
      !$acc loop seq
      do js = 1, jpspec, 1
        itrd = madvtr(js)
        if (itrd /= 0) then
          map(itrd,itrd) = 1
        end if
      enddo
      !$acc end kernels
      #ifdef ARRAYCHECKSUM
        checksum_result = array_xor_checksum(map)
        print*, 'Checksum - map in ', RoutineName, ' = ', checksum_result
      #endif
    end if
    total = SUM(map)
    if (mype == 0 .AND. printstatus >= prstatus_oper) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I6)') 'TOTAL NUMBER OF NONZERO ENTRIES IN JACOBIAN: ', total
      call umprint(ummessage, src='asad_sparse_vars')
    end if
    !$acc kernels
    isx = 0
    !$acc loop gang vector independent collapse(2)
    do idx_2 = LBOUND(nonzero_map, dim=2), UBOUND(nonzero_map, dim=2), 1
      do idx_3 = LBOUND(nonzero_map, dim=1), UBOUND(nonzero_map, dim=1), 1
        nonzero_map(idx_3,idx_2) = 0
      enddo
    enddo
    !$acc loop seq
    do i = 1, jpcspf, 1
      !$acc loop seq
      do j = 1, jpcspf, 1
        if (map(i,j) == 1) then
          isx = isx + 1
          nonzero_map(i,j) = isx
        end if
      enddo
    enddo
    !$acc loop gang vector independent
    do idx_4 = LBOUND(base_tracer, dim=1), UBOUND(base_tracer, dim=1), 1
      base_tracer(idx_4) = 0
    enddo
    !$acc loop seq
    do i = 1, jpcspf, 1
      !$acc loop seq
      do j = 1, jpcspf, 1
        i1 = nonzero_map(j,i)
        if (i1 > 0) then
          base_tracer(i1) = i
        end if
      enddo
    enddo
    !$acc loop gang vector independent
    do i = 1, jpcspf, 1
      activity(i) = SUM(map(:,i)) + SUM(map(i,:))
      reorder(i) = i
    enddo
    !$acc loop seq
    do i = 1, jpcspf - 1, 1
      !$acc loop seq
      do j = i + 1, jpcspf, 1
        if (activity(i) > activity(j)) then
          itemp1 = reorder(i)
          reorder(i) = reorder(j)
          reorder(j) = itemp1
          itemp1 = activity(i)
          activity(i) = activity(j)
          activity(j) = itemp1
        end if
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(nonzero_map)
      print*, 'Checksum - nonzero_map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(base_tracer)
      print*, 'Checksum - base_tracer in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(activity)
      print*, 'Checksum - activity in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(reorder)
      print*, 'Checksum - reorder in ', RoutineName, ' = ', checksum_result
    #endif

    print *, "RW - value of jpcspf = ", jpcspf
    print *, "RW - part 0 - min value of nonzero_map = ", minval(nonzero_map)
    print *, "RW - lbounds nonzero dim 2", lbound(nonzero_map, 2)
    print *, "RW - ubounds nonzero dim 2", ubound(nonzero_map, 2)
    print *, "RW - lbounds nonzero dim 1", lbound(nonzero_map, 1)
    print *, "RW - ubounds nonzero dim 1", ubound(nonzero_map, 1)
    do i = 1, jpcspf, 1
      if (mype == 0 .AND. printstatus >= prstatus_oper) then
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(A,I3,1X,I3,1X,A10,1X,I3)') 'IN ASAD SETUP_SPFULJAC: ', i, reorder(i), specf(reorder(i)), activity(i)
        call umprint(ummessage, src='asad_sparse_vars')
      end if
    enddo
    if (.NOT.ALLOCATED(permute)) then
      ALLOCATE(permute(1:jpcspf,1:jpcspf))
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do idx_5 = LBOUND(permute, dim=2), UBOUND(permute, dim=2), 1
      do idx_6 = LBOUND(permute, dim=1), UBOUND(permute, dim=1), 1
        permute(idx_6,idx_5) = 0
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(permute)
      print*, 'Checksum - permute in ', RoutineName, ' = ', checksum_result
    #endif

    !$acc kernels
    !$acc loop gang vector independent
    do i = 1, jpcspf, 1
      permute(i,reorder(i)) = 1
    enddo
    !$acc end kernels
    nonzero_map_unordered(:,:) = MATMUL(MATMUL(permute, nonzero_map), TRANSPOSE(permute))
    if (ALLOCATED(permute)) then
      DEALLOCATE(permute)
    end if
    if (ALLOCATED(map)) then
      DEALLOCATE(map)
    end if
    !$acc kernels
    total1 = total
    !$acc loop gang vector independent collapse(2)
    do idx_7 = LBOUND(modified_map, dim=2), UBOUND(modified_map, dim=2), 1
      do idx_8 = LBOUND(modified_map, dim=1), UBOUND(modified_map, dim=1), 1
        modified_map(idx_8,idx_7) = nonzero_map_unordered(idx_8 + (LBOUND(nonzero_map_unordered, dim=1) - LBOUND(modified_map, dim=1)),idx_7 + (LBOUND(nonzero_map_unordered, dim=2) - LBOUND(modified_map, dim=2)))
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(modified_map)
      print*, 'Checksum - modified_map in ', RoutineName, ' = ', checksum_result
    #endif

    do kr = 1, jpcspf, 1
      do i = kr + 1, jpcspf, 1
        ikr = modified_map(i,kr)
        if (ikr > 0) then
          do j = kr + 1, jpcspf, 1
            krj = modified_map(kr,j)
            if (krj > 0) then
              ij1 = modified_map(i,j)
              if (ij1 <= 0) then
                total1 = total1 + 1
                if (total1 > spfjsize_max) then
                  errcode = total1
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, '(A,2I4)') 'Total1 exceeded spfjsize_max: ', total1, spfjsize_max
                  call umprint(ummessage, src='asad_sparse_vars')
                  cmessage(:) = 'Increase spfjsize_max'
                  call ereport('SETUP_SPFULJAC', errcode, cmessage)
                end if
                modified_map(i,j) = total1
              end if
            end if
          enddo
        end if
      enddo
    enddo
    if (mype == 0 .AND. printstatus >= prstatus_normal) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,2I4)') 'Initial and final fill-in: ', total, total1
      call umprint(ummessage, src='asad_sparse_vars')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,I4)') 'spfjsize_max is currently set to: ', spfjsize_max
      call umprint(ummessage, src='asad_sparse_vars')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A)') 'Reduce spfjsize_max if much greater than total1 '
      call umprint(ummessage, src='asad_sparse_vars')
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx_9 = LBOUND(nposterms, dim=1), UBOUND(nposterms, dim=1), 1
      nposterms(idx_9) = 0
    enddo
    !$acc loop gang vector independent
    do idx_10 = LBOUND(nnegterms, dim=1), UBOUND(nnegterms, dim=1), 1
      nnegterms(idx_10) = 0
    enddo
    !$acc loop gang vector independent
    do idx_11 = LBOUND(nfracterms, dim=1), UBOUND(nfracterms, dim=1), 1
      nfracterms(idx_11) = 0
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_12 = LBOUND(posterms, dim=2), UBOUND(posterms, dim=2), 1
      do idx_13 = LBOUND(posterms, dim=1), UBOUND(posterms, dim=1), 1
        posterms(idx_13,idx_12) = 0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_14 = LBOUND(negterms, dim=2), UBOUND(negterms, dim=2), 1
      do idx_15 = LBOUND(negterms, dim=1), UBOUND(negterms, dim=1), 1
        negterms(idx_15,idx_14) = 0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_16 = LBOUND(fracterms, dim=2), UBOUND(fracterms, dim=2), 1
      do idx_17 = LBOUND(fracterms, dim=1), UBOUND(fracterms, dim=1), 1
        fracterms(idx_17,idx_16) = 0
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_18 = LBOUND(ffrac, dim=2), UBOUND(ffrac, dim=2), 1
      do idx_19 = LBOUND(ffrac, dim=1), UBOUND(ffrac, dim=1), 1
        ffrac(idx_19,idx_18) = 0.0
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(nposterms)
      print*, 'Checksum - nposterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nnegterms)
      print*, 'Checksum - nnegterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(nfracterms)
      print*, 'Checksum - nfracterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(posterms)
      print*, 'Checksum - posterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(negterms)
      print*, 'Checksum - negterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(fracterms)
      print*, 'Checksum - fracterms in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(ffrac)
      print*, 'Checksum - ffrac in ', RoutineName, ' = ', checksum_result
    #endif

    do jc = 1, ntrf, 1
      itrcr = nltrf(jc)
      do j3 = 1, nmzjac(itrcr), 1
        !$acc kernels
        irj = nzjac1(j3,itrcr)
        !$acc loop gang vector independent
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        !$acc end kernels

        #ifdef ARRAYCHECKSUM
          checksum_result = array_xor_checksum(ij)
          print*, 'Checksum - ij in ', RoutineName, ' = ', checksum_result
          checksum_result = array_xor_checksum(lj)
          print*, 'Checksum - lj in ', RoutineName, ' = ', checksum_result
        #endif

        do jn = 1, 2, 1
          if (lj(jn)) then
            i = nonzero_map(ij(jn),itrcr)
            nnegterms(i) = nnegterms(i) + 1
            if (nnegterms(i) > maxterms) then
              errcode = 1
              cmessage(:) = ' Increase maxterms'
              call ereport('SETUP_SPFULJAC', errcode, cmessage)
            end if
            negterms(i,nnegterms(i)) = irj
          end if
        enddo
        if (nfrpx(irj) == 0) then
          do jn = 3, jpmsp, 1
            if (lj(jn)) then
              i = nonzero_map(ij(jn),itrcr)
              nposterms(i) = nposterms(i) + 1
              if (nposterms(i) > maxterms) then
                errcode = 2
                cmessage(:) = ' Increase maxterms'
                call ereport('SETUP_SPFULJAC', errcode, cmessage)
              end if
              posterms(i,nposterms(i)) = irj
            end if
          enddo
        else
          do jn = 3, jpmsp, 1
            if (lj(jn)) then
              i = nonzero_map(ij(jn),itrcr)
              nfracterms(i) = nfracterms(i) + 1
              if (nfracterms(i) > maxfterms) then
                errcode = 3
                cmessage(:) = ' Increase maxfterms'
                call ereport('SETUP_SPFULJAC', errcode, cmessage)
              end if
              if (nfrpx(irj) + jn - 3 > jpfrpx) then
                cmessage(:) = ' frpx array index > jpfrpx'
                errcode = 3
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(A)') cmessage
                call umprint(ummessage, src='asad_sparse_vars')
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(A,I4,A,I4)') 'irj: ', irj, ' nfrpx(irj): ', nfrpx(irj)
                call umprint(ummessage, src='asad_sparse_vars')
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Write_Stmt
                WRITE(umMessage, '(A,I4,A,I4)') 'i: ', i, ' jn: ', jn
                call umprint(ummessage, src='asad_sparse_vars')
                call ereport('SETUP_SPFULJAC', errcode, cmessage)
              end if
              fracterms(i,nfracterms(i)) = irj
              ffrac(i,nfracterms(i)) = frpx(nfrpx(irj) + jn - 3)
            end if
          enddo
        end if
        if (npdfr(irj,1) /= 0) then
          i1 = npdfr(irj,1)
          i2 = npdfr(irj,2)
          do jn = i1, i2, 1
            isx = ntabpd(jn,1)
            i = nonzero_map(isx,itrcr)
            nfracterms(i) = nfracterms(i) + 1
            if (nfracterms(i) > maxfterms) then
              errcode = 3
              cmessage(:) = ' Increase maxfterms'
              call ereport('SETUP_SPFULJAC', errcode, cmessage)
            end if
            fracterms(i,nfracterms(i)) = irj
            ffrac(i,nfracterms(i)) = ztabpd(jn,1)
          enddo
        end if
      enddo
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine setup_spfuljac
  subroutine spfuljac(n_points, cdt, min_pivot, nonzero_map, spfj)
    use asad_mod, only : base_tracer, ctype, deriv, dpd, dpw, f, ffrac, fracterms, jpcspf, jpfm, jpif, jpmsp, jpspec, linfam, madvtr, moffam, ndepd, ndepw, negterms, nfracterms, njcoth, nmsjac, nnegterms, nodd, nposterms, nsjac1, nstst, ntro3, posterms, prk, spfjsize_max, total, y
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use umPrintMgr, only : printstatus, prstatus_normal, prstatus_oper, ummessage, umprint
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, intent(in) :: cdt
    real, intent(in) :: min_pivot
    integer, dimension(jpcspf,jpcspf), intent(in) :: nonzero_map
    real, dimension(n_points,spfjsize_max), intent(out) :: spfj
    real :: deltt
    CHARACTER(LEN = 2) :: ityped
    integer :: p
    integer :: irj
    integer :: ifamd
    integer :: itrd
    integer :: i
    integer :: jc
    integer :: j3
    integer :: jn
    integer :: js
    integer :: jl
    integer, dimension(jpmsp) :: ij
    logical, dimension(jpmsp) :: lj
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'SPFULJAC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9

    #ifdef ARRAYCHECKSUM
      integer :: checksum_result
    #endif

    print *, "RW - part 4 - min value of nonzero_map = ", minval(nonzero_map)

    print *, "RW - Size of lj:", LBOUND(lj,1), UBOUND(lj,1)
    call umprint(ummessage, src='asad_sparse_vars')

    print *, "RW - jpmsp:", jpmsp
    call umprint(ummessage, src='asad_sparse_vars')

    if (ALLOCATED(deriv)) then
      print *, 'RW - derive allocated'
      call umprint(ummessage, src='asad_mod')
    else
      print *, 'RW - derive not allocated'
      call umprint(ummessage, src='asad_mod')
    end if

    print *, "RW - lbounds prk:", LBOUND(prk, dim=1)
    call umprint(ummessage, src='asad_sparse_vars')
    print *, "RW - lbounds spfj:", LBOUND(spfj, dim=1)
    call umprint(ummessage, src='asad_sparse_vars')
    print *, "RW - lbounds deriv:", LBOUND(deriv, dim=1)
    call umprint(ummessage, src='asad_sparse_vars')


    do jc = 1, nstst, 1
      do j3 = 1, nmsjac(jc), 1
        irj = nsjac1(j3,jc)
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        do jn = 1, 2, 1
          if (lj(jn)) then
            p = nonzero_map(ij(jn),ntro3)
            if (p < LBOUND(spfj,2) .or. p > UBOUND(spfj,2)) then
              print *, "RW ERROR: Out of bounds p = ", p
            endif
            if (irj < LBOUND(prk,2) .or. irj > UBOUND(prk,2)) then
              print *, "RW ERROR: Out of bounds irj = ", irj
            endif
          end if
        enddo
        do jn = 3, jpmsp, 1
          if (lj(jn)) then
            p = nonzero_map(ij(jn),ntro3)
            if (p < LBOUND(spfj,2) .or. p > UBOUND(spfj,2)) then
              print *, "RW ERROR: Out of bounds p = ", p
            endif
            if (irj < LBOUND(prk,2) .or. irj > UBOUND(prk,2)) then
              print *, "RW ERROR: Out of bounds irj = ", irj
            endif
          end if
        enddo
      enddo
    enddo

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(f)
      print*, 'Checksum - f in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(spfj)
      print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
    #endif

    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do jl = 1, jpcspf, 1
      do i = 1, n_points, 1
        if (f(i,jl) < min_pivot) then
          f(i,jl) = min_pivot
        end if
      enddo
    enddo
    deltt = 1.0 / cdt
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(spfj, dim=2), UBOUND(spfj, dim=2), 1
      do idx_1 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
        spfj(idx_1,idx) = 0.0
      enddo
    enddo
    !$acc loop seq
    do i = 1, jpcspf, 1
      !$acc loop gang vector independent
      do idx_2 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
        spfj(idx_2,nonzero_map(i,i)) = -deltt * f(idx_2 + (LBOUND(f, dim=1) - LBOUND(spfj, dim=1)),i)
      enddo
    enddo
    !$acc loop seq
    do p = 1, total, 1
      !$acc loop seq
      do i = nposterms(p), 1, -1
        !$acc loop gang vector independent
        do idx_3 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
          spfj(idx_3,p) = spfj(idx_3,p) + prk(idx_3 + (LBOUND(prk, dim=1) - LBOUND(spfj, dim=1)),posterms(p,i))
        enddo
      enddo
      !$acc loop seq
      do i = nnegterms(p), 1, -1
        !$acc loop gang vector independent
        do idx_4 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
          spfj(idx_4,p) = spfj(idx_4,p) - prk(idx_4 + (LBOUND(prk, dim=1) - LBOUND(spfj, dim=1)),negterms(p,i))
        enddo
      enddo
      !$acc loop seq
      do i = nfracterms(p), 1, -1
        !$acc loop gang vector independent
        do idx_5 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
          spfj(idx_5,p) = spfj(idx_5,p) + ffrac(p,i) * prk(idx_5 + (LBOUND(prk, dim=1) - LBOUND(spfj, dim=1)),fracterms(p,i))
        enddo
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(f)
      print*, 'Checksum - f in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(spfj)
      print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
    #endif

    ! issue with kernel region cover whole section below - solution is to move kernel regions to inner loop
    do jc = 1, nstst, 1
      do j3 = 1, nmsjac(jc), 1
        irj = nsjac1(j3,jc)
        do jn = 1, jpmsp, 1
          ij(jn) = njcoth(irj,jn)
          lj(jn) = ij(jn) /= 0
        enddo
        do jn = 1, 2, 1
          if (lj(jn)) then
            print *, "RW: P = ", nonzero_map(ij(jn),ntro3)
            print *, "RW: ij(jn) = ", ij(jn)
            print *, "RW: ntro3 = ", ntro3

            !$acc kernels
            p = nonzero_map(ij(jn),ntro3)
            ! !$acc loop gang vector independent
            ! do idx_6 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
            !   spfj(idx_6,p) = spfj(idx_6,p) - prk(idx_6 + (LBOUND(prk, dim=1) - LBOUND(spfj, dim=1)),irj) * deriv(idx_6 + (LBOUND(deriv, dim=1) - LBOUND(spfj, dim=1)),jc,1)
            !   ! spfj(idx_6,p) = spfj(idx_6,p) - prk(idx_6, irj)!  * deriv(idx_6,jc,1)
            ! enddo
            spfj(:,p) = spfj(:,p) -                                                &
              prk(:,irj)*deriv(:,jc,1) ! *prkrat(jc)
            !$acc end kernels

            #ifdef ARRAYCHECKSUM
              checksum_result = array_xor_checksum(spfj)
              print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
            #endif

          end if
        enddo

        do jn = 3, jpmsp, 1
          if (lj(jn)) then
            !$acc kernels
            p = nonzero_map(ij(jn),ntro3)
            !$acc loop gang vector independent
            do idx_7 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
              spfj(idx_7,p) = spfj(idx_7,p) + prk(idx_7 + (LBOUND(prk, dim=1) - LBOUND(spfj, dim=1)),irj) * deriv(idx_7 + (LBOUND(deriv, dim=1) - LBOUND(spfj, dim=1)),jc,1)
              ! spfj(idx_7,p) = spfj(idx_7,p) + prk(idx_7,irj) ! * deriv(idx_7,jc,1)
            enddo
            !$acc end kernels

            #ifdef ARRAYCHECKSUM
              checksum_result = array_xor_checksum(spfj)
              print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
            #endif

          end if
        enddo
      enddo
    enddo

    if (ndepw /= 0 .OR. ndepd /= 0) then
      do js = 1, jpspec, 1
        ifamd = moffam(js)
        itrd = madvtr(js)
        ityped(:) = ctype(js)
        if (ifamd /= 0) then
          !$acc kernels
          !$acc loop gang vector independent
          do jl = 1, n_points, 1
            if (ityped(:) == jpfm(:) .OR. ityped(:) == jpif(:) .AND. linfam(jl,itrd)) then
              p = nonzero_map(ifamd,ifamd)
              spfj(jl,p) = spfj(jl,p) - nodd(js) * (dpd(jl,js) + dpw(jl,js)) * y(jl,js)
            end if
          enddo
          !$acc end kernels

          #ifdef ARRAYCHECKSUM
            checksum_result = array_xor_checksum(spfj)
            print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
          #endif

        end if
        if (itrd /= 0) then
          !$acc kernels
          p = nonzero_map(itrd,itrd)
          !$acc loop gang vector independent
          do idx_8 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
            spfj(idx_8,p) = spfj(idx_8,p) - (dpd(idx_8 + (LBOUND(dpd, dim=1) - LBOUND(spfj, dim=1)),js) + dpw(idx_8 + (LBOUND(dpw, dim=1) - LBOUND(spfj, dim=1)),js)) * y(idx_8 + (LBOUND(y, dim=1) - LBOUND(spfj, dim=1)),js)
          enddo
          !$acc end kernels

          #ifdef ARRAYCHECKSUM
            checksum_result = array_xor_checksum(spfj)
            print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
          #endif

        end if
      enddo
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do p = 1, total, 1
      do idx_9 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
        spfj(idx_9,p) = spfj(idx_9,p) / f(idx_9 + (LBOUND(f, dim=1) - LBOUND(spfj, dim=1)),base_tracer(p))
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(spfj)
      print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine spfuljac
  subroutine splinslv2(n_points, bb, xx, min_pivot, max_val, nonzero_map_unordered, modified_map, spfj)
    use asad_mod, only : jpcspf, reorder, spfjsize_max, total
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points,jpcspf), intent(inout) :: bb
    real, dimension(n_points,jpcspf), intent(out) :: xx
    real, intent(in) :: min_pivot
    real, intent(in) :: max_val
    integer, dimension(jpcspf,jpcspf), intent(in) :: nonzero_map_unordered
    integer, dimension(jpcspf,jpcspf), intent(inout) :: modified_map
    real, dimension(n_points,spfjsize_max), intent(inout) :: spfj
    integer :: total1
    integer :: kr
    integer :: i
    integer :: j
    integer :: ikr
    integer :: krj
    integer :: ij
    real, dimension(n_points,jpcspf) :: bb1
    real, dimension(n_points,jpcspf) :: xx1
    real, dimension(n_points) :: pivot
    real, dimension(n_points) :: kfact
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'SPLINSLV2'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: widx1
    integer :: widx1_1

    #ifdef ARRAYCHECKSUM
      integer :: checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(spfj)
      print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(modified_map)
      print*, 'Checksum - modified_map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(pivot)
      print*, 'Checksum - pivot in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(kfact)
      print*, 'Checksum - kfact in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(bb1)
      print*, 'Checksum - bb1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx1)
      print*, 'Checksum - xx1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx)
      print*, 'Checksum - xx in ', RoutineName, ' = ', checksum_result
    #endif

    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do j = 1, total, 1
      do i = 1, n_points, 1
        spfj(i,j) = MIN(MAX(spfj(i,j), -max_val), max_val)
      enddo
    enddo
    total1 = total
    !$acc loop gang vector independent collapse(2)
    do idx = LBOUND(modified_map, dim=2), UBOUND(modified_map, dim=2), 1
      do idx_1 = LBOUND(modified_map, dim=1), UBOUND(modified_map, dim=1), 1
        modified_map(idx_1,idx) = nonzero_map_unordered(idx_1,idx)
      enddo
    enddo
    !$acc loop seq
    do kr = 1, jpcspf, 1
      !$acc loop gang vector independent
      do idx_2 = 1, n_points, 1
        pivot(idx_2) = spfj(idx_2,modified_map(kr,kr))
      enddo
      !$acc loop gang vector independent
      do widx1 = 1, n_points, 1
        if (ABS(pivot(widx1)) > min_pivot) then
          pivot(widx1) = 1.0 / pivot(widx1)
        else
          pivot(widx1) = max_val
        end if
      enddo
      !$acc loop seq
      do i = kr + 1, jpcspf, 1
        ikr = modified_map(i,kr)
        if (ikr > 0) then
          !$acc loop gang vector independent
          do idx_3 = 1, n_points, 1
            kfact(idx_3) = spfj(idx_3,ikr) * pivot(idx_3)
          enddo
          !$acc loop gang vector independent
          do idx_4 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
            spfj(idx_4,ikr) = kfact(idx_4)
          enddo
          !$acc loop seq
          do j = kr + 1, jpcspf, 1
            krj = modified_map(kr,j)
            if (krj > 0) then
              ij = modified_map(i,j)
              if (ij > 0) then
                !$acc loop gang vector independent
                do idx_5 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
                  spfj(idx_5,ij) = spfj(idx_5,ij) - kfact(idx_5) * spfj(idx_5,krj)
                enddo
              else
                total1 = total1 + 1
                modified_map(i,j) = total1
                !$acc loop gang vector independent
                do idx_6 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
                  spfj(idx_6,total1) = -kfact(idx_6) * spfj(idx_6,krj)
                enddo
              end if
            end if
          enddo
        end if
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do j = 1, total1, 1
      do i = 1, n_points, 1
        spfj(i,j) = MIN(MAX(spfj(i,j), -max_val), max_val)
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do i = 1, jpcspf, 1
      do idx_7 = LBOUND(bb1, dim=1), UBOUND(bb1, dim=1), 1
        bb1(idx_7,i) = bb(idx_7,reorder(i))
      enddo
    enddo
    !$acc loop seq
    do kr = 1, jpcspf - 1, 1
      !$acc loop seq
      do i = kr + 1, jpcspf, 1
        ikr = modified_map(i,kr)
        if (ikr > 0) then
          !$acc loop gang vector independent
          do idx_8 = LBOUND(bb1, dim=1), UBOUND(bb1, dim=1), 1
            bb1(idx_8,i) = bb1(idx_8,i) - spfj(idx_8,ikr) * bb1(idx_8,kr)
          enddo
        end if
      enddo
    enddo
    !$acc loop seq
    do kr = jpcspf, 1, -1
      !$acc loop gang vector independent
      do idx_9 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
        xx1(idx_9,kr) = bb1(idx_9,kr)
      enddo
      !$acc loop seq
      do j = kr + 1, jpcspf, 1
        krj = modified_map(kr,j)
        if (krj > 0) then
          !$acc loop gang vector independent
          do idx_10 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
            xx1(idx_10,kr) = xx1(idx_10,kr) - spfj(idx_10,krj) * xx1(idx_10,j)
          enddo
        end if
      enddo
      !$acc loop gang vector independent
      do idx_11 = 1, n_points, 1
        pivot(idx_11) = spfj(idx_11,modified_map(kr,kr))
      enddo
      !$acc loop gang vector independent
      do widx1_1 = 1, n_points, 1
        if (ABS(pivot(widx1_1)) < min_pivot) then
          pivot(widx1_1) = min_pivot
        end if
      enddo
      !$acc loop gang vector independent
      do idx_12 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
        xx1(idx_12,kr) = MIN(MAX(xx1(idx_12,kr), -max_val), max_val) / pivot(idx_12)
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_13 = 1, jpcspf, 1
      do idx_14 = 1, n_points, 1
        xx(idx_14,idx_13) = 0.0
      enddo
    enddo
    !$acc loop seq
    do j = 1, jpcspf, 1
      !$acc loop gang vector independent
      do idx_15 = LBOUND(xx, dim=1), UBOUND(xx, dim=1), 1
        xx(idx_15,reorder(j)) = xx1(idx_15,j)
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(spfj)
      print*, 'Checksum - spfj in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(modified_map)
      print*, 'Checksum - modified_map in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(pivot)
      print*, 'Checksum - pivot in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(kfact)
      print*, 'Checksum - kfact in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(bb1)
      print*, 'Checksum - bb1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx1)
      print*, 'Checksum - xx1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx)
      print*, 'Checksum - xx in ', RoutineName, ' = ', checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine splinslv2
  subroutine spresolv2(n_points, bb, xx, min_pivot, modified_map, spfj)
    use asad_mod, only : jpcspf, reorder, spfjsize_max
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real, dimension(n_points,jpcspf), intent(in) :: bb
    real, dimension(n_points,jpcspf), intent(out) :: xx
    real, intent(in) :: min_pivot
    integer, dimension(jpcspf,jpcspf), intent(in) :: modified_map
    real, dimension(n_points,spfjsize_max), intent(in) :: spfj
    integer :: kr
    integer :: i
    integer :: j
    integer :: krj
    integer :: ikr
    real, dimension(n_points,jpcspf) :: bb1
    real, dimension(n_points,jpcspf) :: xx1
    real, dimension(n_points) :: pivot
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'SPRESOLV2'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: widx1

    #ifdef ARRAYCHECKSUM
      integer :: checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(bb1)
      print*, 'Checksum - bb1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx1)
      print*, 'Checksum - xx1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(pivot)
      print*, 'Checksum - pivot in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx)
      print*, 'Checksum - xx in ', RoutineName, ' = ', checksum_result
    #endif

    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do i = 1, jpcspf, 1
      do idx = LBOUND(bb1, dim=1), UBOUND(bb1, dim=1), 1
        bb1(idx,i) = bb(idx,reorder(i))
      enddo
    enddo
    !$acc loop seq
    do kr = 1, jpcspf - 1, 1
      !$acc loop seq
      do i = kr + 1, jpcspf, 1
        ikr = modified_map(i,kr)
        if (ikr > 0) then
          !$acc loop gang vector independent
          do idx_1 = LBOUND(bb1, dim=1), UBOUND(bb1, dim=1), 1
            bb1(idx_1,i) = bb1(idx_1,i) - spfj(idx_1,ikr) * bb1(idx_1,kr)
          enddo
        end if
      enddo
    enddo
    !$acc loop seq
    do kr = jpcspf, 1, -1
      !$acc loop gang vector independent
      do idx_2 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
        xx1(idx_2,kr) = bb1(idx_2,kr)
      enddo
      !$acc loop seq
      do j = kr + 1, jpcspf, 1
        krj = modified_map(kr,j)
        if (krj > 0) then
          !$acc loop gang vector independent
          do idx_3 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
            xx1(idx_3,kr) = xx1(idx_3,kr) - spfj(idx_3,krj) * xx1(idx_3,j)
          enddo
        end if
      enddo
      !$acc loop gang vector independent
      do idx_4 = 1, n_points, 1
        pivot(idx_4) = spfj(idx_4,modified_map(kr,kr))
      enddo
      !$acc loop gang vector independent
      do widx1 = 1, n_points, 1
        if (ABS(pivot(widx1)) < min_pivot) then
          pivot(widx1) = min_pivot
        end if
      enddo
      !$acc loop gang vector independent
      do idx_5 = LBOUND(xx1, dim=1), UBOUND(xx1, dim=1), 1
        xx1(idx_5,kr) = xx1(idx_5,kr) / pivot(idx_5)
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_6 = 1, jpcspf, 1
      do idx_7 = 1, n_points, 1
        xx(idx_7,idx_6) = 0.0
      enddo
    enddo
    !$acc loop seq
    do i = 1, jpcspf, 1
      !$acc loop gang vector independent
      do idx_8 = LBOUND(xx, dim=1), UBOUND(xx, dim=1), 1
        xx(idx_8,reorder(i)) = xx1(idx_8,i)
      enddo
    enddo
    !$acc end kernels

    #ifdef ARRAYCHECKSUM
      checksum_result = array_xor_checksum(bb1)
      print*, 'Checksum - bb1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx1)
      print*, 'Checksum - xx1 in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(pivot)
      print*, 'Checksum - pivot in ', RoutineName, ' = ', checksum_result
      checksum_result = array_xor_checksum(xx)
      print*, 'Checksum - xx in ', RoutineName, ' = ', checksum_result
    #endif

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine spresolv2

end module asad_sparse_vars
