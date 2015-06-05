!***********************************************************************

!     file  lusol8b.f90

!     lu8adc   lu8adr   lu8dlc   lu8dlr   lu8mod   lu8rpr

!     These routines call lu7asv and lu7bak in lusol7b.f90.

! 08 Jun 2004: lusol8b.f is essentially lu8b.for from VMS days.
!              integer*4 changed to integer  .
! 15 Sep 2004: Test nout. gt. 0 to protect write statements.
! 27 Mar 2015: Converted to F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module lusol8b
  use  lusol_precision, only : ip, rp
  use  lusol
  use  lusol6b
  use  lusol7b
  implicit none
  private
  public    :: lu8adc, lu8adr, lu8dlc, lu8dlr, lu8mod, lu8rpr
  intrinsic :: abs, int, max, min, real

contains

    subroutine lu8adc( mode, m, n, v, w, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform, diag, vnorm )
    
    implicit real(rp) (a-h,o-z)
    
    integer(ip),   intent(in)    :: mode, m, n, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena), v(m), w(n), vnorm, diag
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform

!     ------------------------------------------------------------------
!     lu8adc  updates the LU factorization  A = L*U  when the vector
!     a(new)  is added to become the  n-th  column of  A.
!     The dimension of  A  is assumed to increase from  n-1  to  n
!     (not from  n  to  n+1).

!     If  mode = 1,  v(*)  must contain  a(new).
!     If  mode = 2,  v(*)  must satisfy  L*v = a(new).
!     On exit,  L*v = a(new)  in both cases.

!     The array  w(*)  is not used or altered.

!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.

!     On exit:
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  7  if the update was not completed (lack of storage).

!     -- Feb 1985: Last  F66 version.
!     12 May 1988: First F77 version.  Now uses lu8rpc.
!     ------------------------------------------------------------------    
    integer(ip)            :: nout, lprint, nrank, mode2
    integer(ip), parameter :: izero = 0
    
    nout   = luparm(1)
    lprint = luparm(2)
    nrank  = luparm(16)

! Set  locc(n)  in case the user forgot.
! Then let lu8rpc do the job (but stop it from printing messages).

    locc(n)   =  0
    q(n)     =  n
    luparm(2) = -1

    call lu8rpc( izero, mode, m, n, n, v, w, &
                 lena, luparm, parmlu, &
                 a, indc, indr, p, q, &
                 lenc, lenr, locc, locr, &
                 inform, diag, vnorm )

    if (inform <= 0) then
        if (nrank == n - 1) then
            if (nout > 0  .AND.  lprint >= 0) &
            write(nout, 1100) m, n, diag, vnorm
        end if
    else if (inform == 7) then
        if (nout > 0  .AND.  lprint >= 0) &
        write(nout, 1700) lena
    end if

    luparm(2) = lprint
    return

    1100 format(/ ' lu8adc  warning.  Rank did not increase', &
    ' after adding a column.' &
    / ' m =', I8, '    n =', i8, &
    '    diag =', 1p, e12.2, '    vnorm =', e12.2)
    1700 format(/ ' lu8adc  error...  Insufficient storage.', &
    '    lena =', I8)

    end ! subroutine lu8adc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine lu8adr( m, n, w, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform, diag )
                       
    implicit real(rp) (a-h,o-z)
    
    integer(ip),   intent(in)    :: m, n, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena), diag
    real(rp),      intent(in)    :: w(n)
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform

!     ------------------------------------------------------------------
!     lu8adr  updates the LU factorization  A = L*U  when the vector
!     w(*)  is added to become the  m-th  row of  A.
!     The dimension of  A  is assumed to increase from  m-1  to  m
!     (not from  m  to  m+1).

!     w(*)  is not altered.

!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.

!     On exit:
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  7  if the update was not completed (lack of storage).

!     -- Feb 1985: Last  F66 version.
!     17 May 1988: First F77 version.
!     ------------------------------------------------------------------    
    integer(ip)            :: nout, lprint, nrank, lenl, lenu, lrow, &
                              nrank1, minfre, nfree, j, k, kfirst, lenw
    integer(rp), parameter :: izero = 0
    
    real(rp)               :: small
    
    real(rp), parameter    :: zero = 0.0
          
    nout   = luparm(1)
    lprint = luparm(2)
    nrank  = luparm(16)
    lenl   = luparm(23)
    lenu   = luparm(24)
    lrow   = luparm(25)
    small  = parmlu(3)
    diag   = zero
    nrank1 = nrank + 1

! Compress row file if necessary.

    minfre = n
    nfree  = lena - lenl - lrow
    if (nfree >= minfre) go to 100
    call lu1rec( m, .TRUE. , luparm, lrow, lena, a, indr, lenr, locr )
    nfree  = lena - lenl - lrow
    if (nfree < minfre) go to 970

! Pack the nonzeros of  w  at the end of the row file.
! Go backwards in order to set  kfirst.

    100 locr(m) = lrow + 1
    p(m)   = m
    lenr(m) = 0

    do 120 k = n, 1, -1
        j          = q(k)
        if (abs( w(j) ) <= small) go to 120
        kfirst     = k
        lrow       = lrow + 1
        a(lrow)    = w(j)
        indr(lrow) = j
    120 END DO

    lenw    = lrow + 1 - locr(m)
    lenu    = lenu + lenw
    lenr(m) = lenw

!     ------------------------------------------------------------------
!     Triangularize the new row.
!     ------------------------------------------------------------------

    if (lenw > 0) then

    ! Swap the new row into position nrank1.

        p(m)      = p(nrank1)
        p(nrank1) = m

    ! Perform a forward sweep to eliminate subdiagonal elements.

        if (kfirst <= nrank1) then
            call lu7for( m, n, kfirst, nrank1, &
                         lena, luparm, parmlu, &
                         lenl, lenu, lrow, &
                         a, indc, indr, p, q, lenr, locc, locr, &
                         inform, diag )
            if (inform == 7) go to 970
        end if

    ! See if the rank increases.

        if (nrank < n) then
            nrank  = nrank1
            call lu7rnk( m, n, izero, lena, parmlu, &
                         lenl, lenu, lrow, nrank, &
                         a, indc, indr, p, q, lenr, locc, locr, &
                         inform, diag )
        end if
    end if

!     ------------------------------------------------------------------
!     Set inform for normal exit.
!     ------------------------------------------------------------------
    if (nrank /= nrank1) then
        inform = 0
        if (nrank == m - 1) then
            if (nout > 0  .AND.  lprint >= 0) &
            write(nout, 1100) m, n, diag
        end if
    else
        inform = 1
    end if
    go to 990

! Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena

! Exit.

    990 luparm(10) = inform
    luparm(15) = luparm(15) + 1
    luparm(16) = nrank
    luparm(23) = lenl
    luparm(24) = lenu
    luparm(25) = lrow
    return

    1100 format(/ ' lu8adr  warning.  Rank did not increase', &
    ' after adding a row.' &
    / ' m =', i8, '    n =', i8, &
    '    diag =', 1p, e12.2)
    1700 format(/ ' lu8adr  error...  Insufficient storage.', &
    '    lena =', i8)

    end ! subroutine lu8adr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu8dlc( m, n, jdel, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform )
    
    implicit real(rp) (a-h,o-z)
        
    integer(ip),   intent(in)    :: m, n, jdel, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena)
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform

!     ------------------------------------------------------------------
!     lu8dlc  updates the LU factorization  A = L*U  when column  jdel
!     is deleted.  The dimension of  A  decreases from  n  to  n-1.
!     The column permutation  q(*)  and the column indices in  U
!     are altered accordingly.

!     NOTE:  The calling program must change  n  to  n-1.

!     In some cases it may be more efficient to use  lu8rpc  to
!     replace column  jdel  by zero, leaving it in the current position.

!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.

!     On exit:
!     inform = -1  if the rank of U decreased by 1.
!     inform =  0  if the rank of U stayed the same.
!     inform =  7  if the update was not completed (lack of storage).

!     -- Feb 1985: Last  F66 version.
!     17 May 1988: First F77 version.
!     ------------------------------------------------------------------   
    integer(ip)            :: nout, lprint, nrank, lenl, lenu, lrow, &
                              n1, k, kdel, lr1, lr2, &
                              i, j, l 
    integer(ip), parameter :: izero = 0
                              
    nout   = luparm(1)
    lprint = luparm(2)
    nrank  = luparm(16)
    lenl   = luparm(23)
    lenu   = luparm(24)
    lrow   = luparm(25)
    n1     = n - 1
    if (jdel < 1) go to 980
    if (jdel > n) go to 980

! Remove column jdel from U, and set kdel so that q(kdel) = jdel.

    call lu7zap( m, n, jdel, kdel, &
                 lena, lenu, lrow, nrank, &
                 a, indr, p, q, lenr, locr )

! Renumber columns of  U  that are to the right of column  jdel.
! In effect, those columns are shifted one place to the left.

    if (jdel < n) then
        do k = 1, nrank
            i      = p(k)
            lr1    = locr(i)
            lr2    = lr1 + lenr(i) - 1

            do l = lr1, lr2
                j     = indr(l)
                if (j > jdel) indr(l) = j - 1
            end do
        end do

        do k = 1, n
            j     = q(k)
            if (j > jdel) q(k) = j - 1
        end do
    end if

! Perform cyclic permutations to move column kdel to the end
! and the corresponding row to position nrank.
! Then eliminate the resulting row spike.

    call lu7cyc( kdel, nrank, p )
    call lu7cyc( kdel, n    , q )

    call lu7for( m, n1, kdel, nrank, &
                 lena, luparm, parmlu, &
                 lenl, lenu, lrow, &
                 a, indc, indr, p, q, lenr, locc, locr, &
                 inform, diag )
    if (inform == 7) go to 970

! See if the rank decreased.

    call lu7rnk( m, n1, izero, lena, parmlu, &
                 lenl, lenu, lrow, nrank, &
                 a, indc, indr, p, q, lenr, locc, locr, &
                 inform, diag )
    go to 990

! Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena
    go to 990

! jdel  is out of range.

    980 inform = 8
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1800) m, n, jdel

! Exit.

    990 luparm(10) = inform
    luparm(15) = luparm(15) + 1
    luparm(16) = nrank
    luparm(23) = lenl
    luparm(24) = lenu
    luparm(25) = lrow
    return

    1700 format(/ ' lu8dlc  error...  Insufficient storage.', &
    '    lena =', i8)
    1800 format(/ ' lu8dlc  error...  jdel  is out of range.', &
    '    m =', i8, '    n =', i8, '    jdel =', i8)

    end ! subroutine lu8dlc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu8dlr( mode, m, n, idel, v, w, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform )

    implicit real(rp) (a-h,o-z)

    integer(ip),   intent(in)    :: mode, m, n, idel, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena), v(m), w(n)
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform
    
!     ------------------------------------------------------------------
!     lu8dlr  updates the LU factorization  A = L*U  when row  idel
!     (the vector  w) is deleted from  A. The update is implemented as
!     the rank-one modification

!           A(new)  =  A  -  e(idel) * w',

!     followed by a renumbering that makes row  idel  the last row of  A
!     and shifts rows  idel + 1,  idel + 2,  ...,  m  one place up.
!     Thus, row  idel  is replaced by the zero vector (rather than being
!     deleted), and is then cyclically permuted to the bottom of  A.
!     The dimensions of  A  do not alter, but  A  and  U  may become
!     singular.

!     If  mode = 1,  the old row is assumed to be unknown.  It will be
!                    computed from the LU factors of  A.
!     If  mode = 2,  w(*)  must contain the old row.

!     v(*)  is a work array of length  m.

!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform = 0), this will again be true.

!     Note --- significant overhead is involved in permuting row  idel
!     to the bottom.  In some cases it may be better to use  lu8rpr  to
!     replace row  idel  by zero, leaving it in the current position.
!     The growth of nonzeros in  L  and  U  is identical, but less
!     housekeeping is required than with  lu8dlr.


!     On exit:
!     inform = -1  if the rank of U decreased by 1.
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  7  if the update was not completed (lack of storage).

!     -- Feb 1985: Last  F66 version.
!     18 May 1988: First F77 version.
!     ------------------------------------------------------------------

    real(rp),    parameter :: zero = 0.0,  one = 1.0
    integer(ip)            :: nout, lprint, lenl, mode2 ,&
                              i, k, l, l1
                      
    nout   = luparm(1)
    lprint = luparm(2)
    if (idel < 1) go to 980
    if (idel > m) go to 980

    if (mode == 1) then

    ! Compute row idel as the vector w = A(transpose) * e(idel).

        do i = 1, m
            v(i)  = zero
        end do
        v(idel) = one
        mode2 = 6
        call lu6mul( mode2, m, n, v, w, lena, luparm, parmlu, &
                     a, indc, indr, p, q, lenc, lenr, locc, locr )
    end if

! Set up the required vectors and do the rank-one mod
! (but don't let lu8mod print anything).

    do i = 1, m
        v(i)  = zero
    end do

    v(idel)   =   one
    beta      = - one
    luparm(2) = - 1

    mode2 = 1
    call lu8mod( mode2, m, n, beta, v, w, &
                 lena, luparm, parmlu, &
                 a, indc, indr, p, q, lenc, lenr, locc, locr, &
                 inform )
    if (inform == 7) go to 970

!     ------------------------------------------------------------------
!     Permute the deleted row to the bottom.
!     ------------------------------------------------------------------
    if (idel < m) then
        call lu7cyc( idel, m, lenr )
        call lu7cyc( idel, m, locr )

        do 320 k = 1, m
            i     = p(k)
            if (i < idel) go to 320
            p(k) = i - 1
            if (i == idel) p(k) = m
        320 END DO

        lenl   = luparm(23)
        l1     = lena + 1 - lenl

        do 400 l = l1, lena
            i       = indc(l)
            if (i < idel) go to 350
            indc(l) = i - 1
            if (i == idel) indc(l) = m

            350 i       = indr(l)
            if (i < idel) go to 400
            indr(l) = i - 1
            if (i == idel) indr(l) = m
        400 END DO
    end if

    go to 990

! Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena
    go to 990

! idel  is out of range.

    980 inform = 8
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1800) m, n, idel

! Exit.

    990 luparm(2)  = lprint
    luparm(10) = inform
    return

    1700 format(/ ' lu8dlr  error...  Insufficient storage.', &
    '    lena =', i8)
    1800 format(/ ' lu8dlr  error...  idel  is out of range.', &
    '    m =', i8, '    n =', i8, '    idel =', i8)

    end ! subroutine lu8dlr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu8mod( mode, m, n, beta, v, w, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform )

    implicit real(rp) (a-h,o-z)

    integer(ip),   intent(in)    :: mode, m, n, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena), v(m), w(n)
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform
    
!     ------------------------------------------------------------------
!     lu8mod  updates the LU factorization  A = L*U  when the
!     m by n matrix  A  is subjected to a rank-one modification
!     to become
!                    A  +  beta * v * w(transpose)
!     for the given scalar  beta  and vectors  v, w.

!     If mode = 1, v(*) must contain the vector  v.
!                  On exit, v(*) will contain  y  satisfying  L*y = v.

!     If MODE = 2, v(*) must already contain   y  satisfying  L*y = v.

!     In both cases, w(*) is unaltered.

!     On entry, the elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.

!     On exit:
!     inform = -1  if the rank of U decreased by 1.
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  7  if the update was not completed (lack of storage).

!     -- Mar 1985: Last  F66 version.
!     18 May 1988: First F77 version.
!     ------------------------------------------------------------------
    integer(ip)    :: nout, lprint, nrank, lenl, lenu, &
                      lrow, utol1, nrank0, i, j, k, iv, iw, j1, jelm, &
                      jrep, &! jrep is never defined but it is used!!!
                      jsing, jw, kfirst, klast, l1, lenw, lfree, &
                      lw, lw1, lw2, marker, minfre, nfree
    real(rp)       :: small, cv
    logical        :: singlr

    nout   = luparm(1)
    lprint = luparm(2)
    nrank  = luparm(16)
    lenl   = luparm(23)
    lenu   = luparm(24)
    lrow   = luparm(25)
    small  = parmlu(3)
    utol1  = parmlu(4)
    nrank0 = nrank

! If necessary, solve  L*v(new) = v.

    if (mode == 1) then
        call lu6sol( mode, m, n, v, w, lena, luparm, parmlu, &
                     a, indc, indr, p, q, lenc, lenr, locc, locr, &
                     inform )
    end if

!     ------------------------------------------------------------------
!     Find the first nonzero in  w  (in pivotal column order).
!     ------------------------------------------------------------------
    do k = 1, n
        kfirst = k
        j      = q(k)
        if (abs( w(j) ) > small) go to 120
    end do
    go to 900

!     ------------------------------------------------------------------
!     Eliminate any nonzeros in  v  below the trapezoid.
!     ------------------------------------------------------------------
    120 if (nrank < m) then
        nrank  = nrank + 1
        jelm   = 0

        call lu7elm( m, n, jelm, v, &
                     lena, luparm, parmlu, &
                     lenl, lenu, lrow, nrank, &
                     a, indc, indr, p, q, lenr, locc, locr, &
                     inform, diag )
        if (inform == 7) go to 970

        if (inform == 0) nrank = nrank - 1
    end if

!     ------------------------------------------------------------------
!     Find the last nonzero in  v  (in pivotal row order).
!     ------------------------------------------------------------------
    do k = nrank, 1, -1
        klast  = k
        i      = p(k)
        if (abs( v(i) ) > small) go to 220
    end do
    go to 900

!     ------------------------------------------------------------------
!     Perform a backward sweep of eliminations to reduce part of  v
!     to a multiple of the unit vector  e(iw),  where  iw = p(klast).
!     Elements  p(kfirst+1),  p(kfirst+2),  ...,  p(klast)  of  v
!     are involved.
!     L, U  and  p  are updated accordingly.
!     U  will then be trapezoidal except for row  iw = p(klast).
!     ------------------------------------------------------------------
    220 if (kfirst + 1  <  klast) then
        call lu7bak( m, n, kfirst, klast, v, &
                    lena, luparm, parmlu, &
                    lenl, lenu, lrow, &
                    a, indc, indr, p, q, lenr, locc, locr, &
                    inform )
        if (inform /= 0) go to 970
    end if

!     ------------------------------------------------------------------
!     Pack the nonzeros of  w  in pivotal order in front of  L.
!     (We will treat the packed  w  much like a normal row of  U.)
!     Set markers on  w  and initialize the corresponding
!     elements of  indc(*),  which are later used by  lu7asv.
!     ------------------------------------------------------------------

    lfree  = lena - lenl
    minfre = n + 1 - kfirst
    nfree  = lfree - lrow
    if (nfree >= minfre) go to 310
    call lu1rec( m, .TRUE. , luparm, lrow, lena, a, indr, lenr, locr )
                 nfree  = lfree - lrow
    if (nfree < minfre) go to 970

    310 lw     = lfree + 1

    do 320   k  = n, kfirst, -1
        j        = q(k)
        if (abs( w(j) ) <= small) go to 320
        lw       = lw - 1
        a(lw)    = w(j)
        indr(lw) = j
        indc(lw) = 0
        locc(j)  = lw
    320 END DO

    lw1    = lw
    lw2    = lfree
    lenw   = lw2 + 1 - lw1
    lfree  = lfree - lenw

!     ------------------------------------------------------------------
!     Add multiples of  w  to the first  kfirst  rows of  U.
!     (This does not alter the trapezoidal form of  U.)
!     ------------------------------------------------------------------

    do 450 k  = 1, kfirst
        iv     = p(k)
        cv     = v(iv)
        if (abs( cv ) <= small) go to 450

    !        ===============================================================
    !        Compress storage if necessary, so there will be room if
    !        row  iv  has to be moved to the end.
    !        ===============================================================
        minfre = n
        nfree  = lfree - lrow
        if (nfree >= minfre) go to 420
        call lu1rec( m, .TRUE. , luparm, lrow, lena, a,indr,lenr,locr )
                     nfree  = lfree - lrow
        if (nfree < minfre) go to 970

    !        ===============================================================
    !        Set  v  =  v  +  wmult * w.
    !        ===============================================================
        420 wmult  = beta * cv
        call lu7asv( m, n, iv, lenw, lw1, lw2, k, wmult, &
                     lena, luparm, parmlu, &
                     lenu, lrow, &
                     a, indc, indr, lenr, locc, locr )
    450 END DO

!     ------------------------------------------------------------------
!     Add a multiple of  w  to row  iw  of  U.
!     ------------------------------------------------------------------

    if (kfirst < klast) then
        minfre = n
        nfree  = lfree - lrow
        if (nfree >= minfre) go to 500
        call lu1rec( m, .TRUE. , luparm, lrow, lena, a,indr,lenr,locr )
                     nfree  = lfree - lrow
        if (nfree < minfre) go to 970

        500 iw     = p(klast)
        marker = m + 1
        wmult  = beta * v(iw)
        call lu7asv( m, n, iw, lenw, lw1, lw2, marker, wmult, &
                     lena, luparm, parmlu, &
                     lenu, lrow, &
                     a, indc, indr, lenr, locc, locr )
    end if

!     ------------------------------------------------------------------
!     Cancel the markers on  w.
!     ------------------------------------------------------------------
    do lw  = lw1, lw2
        jw       = indr(lw)
        locc(jw) = 0
    end do

!     ------------------------------------------------------------------
!     Apply a forward sweep to eliminate the nonzeros in row  iw.
!     ------------------------------------------------------------------

    if (kfirst > klast) then
        if (klast < nrank) go to 900
    else
        call lu7for( m, n, kfirst, klast, &
                     lena, luparm, parmlu, &
                     lenl, lenu, lrow, &
                     a, indc, indr, p, q, lenr, locc, locr, &
                     inform, diag )
        if (inform == 7) go to 970
    end if

!     ------------------------------------------------------------------
!     Test for singularity in column klast (if klast .le. nrank).
!     The code is similar to part of lu8rpc with klast in place of krep.
!     ------------------------------------------------------------------

    if (klast <= nrank) then
        diag   = zero
        iw     = p(klast)
        singlr = lenr(iw) == 0

        if ( .NOT. singlr) then
            l1     = locr(iw)
            j1     = indr(l1)
            singlr = j1 /= jrep        ! BUG  HERE: jrep never defined

            if ( .NOT. singlr) then
                diag   = a(l1)
                singlr = abs( diag ) <= utol1
            end if
        end if

        if ( singlr  .AND.  klast < nrank ) then

        ! Perform cyclic permutations to move column klast
        ! to the end and the corresponding row to position nrank.
        ! Then eliminate the resulting row spike.

            call lu7cyc( klast, nrank, p )
            call lu7cyc( klast, n    , q )

            call lu7for( m, n, klast, nrank, &
                         lena, luparm, parmlu, &
                         lenl, lenu, lrow, &
                         a, indc, indr, p, q, lenr, locc, locr, &
                         inform, diag )
            if (inform == 7) go to 970
        end if

    ! Find the best column to be in position nrank.
    ! If nothing satisfactory exists, nrank will be decreased.

        jsing  = 0
        call lu7rnk( m, n, jsing, lena, parmlu, &
                     lenl, lenu, lrow, nrank, &
                     a, indc, indr, p, q, lenr, locc, locr, &
                     inform, diag )
    end if

!     ------------------------------------------------------------------
!     Set inform for exit.
!     ------------------------------------------------------------------

    900 if (nrank == nrank0) then
        inform =  0
    else if (nrank < nrank0) then
        inform = -1
    else
        inform =  1
    end if
    go to 990

! Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena

! Exit.

    990 luparm(10) = inform
    luparm(15) = luparm(15) + 1
    luparm(16) = nrank
    luparm(23) = lenl
    luparm(24) = lenu
    luparm(25) = lrow
    return

    1700 format(/ ' lu8mod  error...  Insufficient storage.', &
    '    lena =', i8)

    end ! subroutine lu8mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu8rpr( mode1, mode2, m, n, irep, v, w, wnew, &
                       lena, luparm, parmlu, &
                       a, indc, indr, p, q, &
                       lenc, lenr, locc, locr, &
                       inform )

    implicit real(rp) (a-h,o-z)
    
    integer(ip),   intent(in)    :: mode1, mode2, m, n, irep, lena
    integer(ip),   intent(inout) :: luparm(30)
    real(rp),      intent(inout) :: parmlu(30), a(lena), v(m), &
                                    w(n), wnew(n)
    integer(ip),   intent(inout) :: indc(lena), indr(lena), &
                                    p(m)     , q(n)     , &
                                    lenc(n)   , lenr(m)   , &
                                    locc(n)   , locr(m)
    integer(ip),   intent(out)   :: inform
!     ------------------------------------------------------------------
!     lu8rpr  updates the LU factorization  A = L*U  when row  irep
!     (a vector  w(old)) is replaced by some vector  w(new).
!     The update is implemented as the rank-one modification

!           A(new)  =  A  -  e(irep) * ( w(old) - w(new) )'

!     with variations determined by  mode1  and  mode2  as follows.

!     mode1
!     -----
!       0     w(old)  is assumed to be zero.
!             w(*)    need not be set on entry, but will be altered.

!       1     w(old)  is assumed to be unknown.  it will be computed
!                     from the LU factors of  A.
!             w(*)    need not be set on entry, but will be altered.

!       2     w(*)    must contain  w(old).
!                     On exit, it will contain  w(old) - w(new).

!       3     w(*)    must contain  w(old) - w(new).  It is not altered.
!                     wnew(*) is not used.

!       4     w(*)    must contain  w(new) - w(old).  It is not altered.
!                     wnew(*) is not used.

!     If  mode1 = 3 or 4,  mode2  is not used.  It may be set to 0 or 1.
!     Otherwise,  mode2  is used as follows.

!     mode2
!     -----
!       0     w(new)  is assumed to be zero.   wnew(*)  is not used.

!       1     wnew(*) must contain the new row.

!     v(*)    is a work array of length  m.
!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.


!     On exit:
!     inform = -1  if the rank of U decreased by 1.
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  7  if the update was not completed (lack of storage).
!     inform =  8  if jrep is not between 1 and n.

!     -- Mar 1985: Last  F66 version.
!     18 May 1988: First F77 version.
!     ------------------------------------------------------------------
    integer(ip)                 :: nout, lprint, i, j, mode
    real(rp),   parameter       :: zero = 0.0,  one = 1.0

    nout   = luparm(1)
    lprint = luparm(2)
    if (irep < 1) go to 980
    if (irep > m) go to 980

    if (mode1 == 0) then

    ! The old row  irep  is zero.

        do j = 1, n
            w(j) = zero
        end do
    else if (mode1 == 1) then

    ! Compute row irep as the vector  w  =  A(transpose) * e(irep).

        do i = 1, m
            v(i)  = zero
        end do
        v(irep) = one
        mode = 6
        call lu6mul( mode, m, n, v, w, lena, luparm, parmlu, &
                     a, indc, indr, p, q, lenc, lenr, locc, locr )
    end if

    if (mode1 <= 2) then

    ! Compute the difference except if  wnew = 0.

        if (mode2 > 0) then
            do j = 1, n
                w(j)  = w(j) - wnew(j)
            end do
        end if
    end if

!     ------------------------------------------------------------------
!     Set  v = a unit vector  and do the rank-one modification
!     (but don't let lu8mod print anything).
!     ------------------------------------------------------------------

    do i = 1, m
        v(i)  = zero
    end do

    v(irep)   =   one
    beta      = - one
    if (mode1 == 4) beta = one
    luparm(2) = - 1

    mode = 1
    call lu8mod( mode, m, n, beta, v, w, &
                 lena, luparm, parmlu, &
                 a, indc, indr, p, q, lenc, lenr, locc, locr, &
                 inform )
    if (inform == 7) go to 970
    go to 990

! Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena
    go to 990

! irep  is out of range.

    980 inform = 8
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1800) m, n, irep

! Exit.

    990 luparm(2)  = lprint
    luparm(10) = inform
    return

    1700 format(/ ' lu8rpr  error...  Insufficient storage.', &
    '    lena =', i8)
    1800 format(/ ' lu8rpr  error...  irep  is out of range.', &
    '    m =', i8, '    n =', i8, '    irep =', i8)

    end ! subroutine lu8rpr

end module lusol8b