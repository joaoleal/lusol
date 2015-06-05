!***********************************************************************
!     File lusol_util.f90.

!     This file contains most of the sparse LU package LUSOL
!     (the parts needed by MINOS, SQOPT and SNOPT).
!     The parts included are

!        27HEAD.f    These comments.
!        lusol1.f    Factor a given matrix A from scratch (lu1fac).
!        lusol2.f    Heap-management routines for lu1fac.
!        lusol6a.f   Solve with the current LU factors.
!        lusol7a.f   Utilities for all update routines.
!        lusol8a.f   Replace a column (Bartels-Golub update).

!***********************************************************************

!     File lusol1.f

!     lu1fac   lu1fad   lu1gau   lu1mar   lu1mRP   lu1mCP   lu1mSP
!     lu1pen   lu1mxc   lu1mxr   lu1or1   lu1or2   lu1or3   lu1or4
!     lu1pq1   lu1pq2   lu1pq3   lu1rec   lu1slk
!     lu1ful   lu1DPP   lu1DCP

! 26 Apr 2002: TCP implemented using heap data structure.
! 01 May 2002: lu1DCP implemented.
! 07 May 2002: lu1mxc must put 0.0 at top of empty columns.
! 09 May 2002: lu1mCP implements Markowitz with cols searched
!              in heap order.
!              Often faster (searching 20 or 40 cols) but more dense.
! 11 Jun 2002: TRP implemented.
!              lu1mRP implements Markowitz with Threshold Rook Pivoting.
!              lu1mxc maintains max col elements.  (Previously lu1max.)
!              lu1mxr maintains max row elements.
! 12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus).
!              Disabled it for the moment.  (Use lu1mar + TCP.)
! 14 Dec 2002: TSP implemented.
!              lu1mSP implements Markowitz with
!              Threshold Symmetric Pivoting.
! 07 Mar 2003: character*1, character*2 changed to f90 form.
!              Comments changed from * in column to ! in column 1.
!              Comments kept within column 72 to avoid compiler warning.
! 19 Dec 2004: Hdelete(...) has new input argument Hlenin.
! 21 Dec 2004: Print Ltol and Lmax with e10.2 instead of e10.1.
! 26 Mar 2006: lu1fad: Ignore nsing from lu1ful.
!              lu1DPP: nsing redefined (but not used by lu1fad).
!              lu1DCP: nsing redefined (but not used by lu1fad).
! 30 Mar 2011: In lu1pen, loop 620 revised following advice from Ralf �stermark,
!              School of Business and Economics at �bo Akademi University.
!              The previous version produced a core dump on the Cray XT.
! 06 Jun 2013: sn28lusol.f90 of 03 Apr 2013 contains updated f90 lu1mxr
!              to improve the efficiency of TRP.  Adapt this for lusol.f
!              = mi27lu.f and sn27lu.f.  The new lu1fad and lu1mxr need
!              3 extra work vectors that are implemented as automatic arrays
!              (not a feature of f77) but this should be ok as long as we
!              use an f90 compiler.
!              Ding Ma and Michael Saunders, Stanford University.
! 07 Jul 2013: sn27lusol.f derived from sn28lu90.f.
!              This version assumes the 3 work vectors for lu1mxr
!              are passed to lu1fac (so lu1fac has 3 extra parameters).
!              The calls to lu1fac in subroutine s2BLU, file sn25bfac
!              must borrow storage from indr, indc (and hence a,
!              because a, indc, indr must seem to have the same length).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu1rec( n, reals, luparm, ltop, &
    lena, a, ind, len, loc )

    logical ::            reals
    integer ::            luparm(30), ltop
    double precision ::   a(lena)
    integer ::            ind(lena), len(n)
    integer ::            loc(n)

!     ------------------------------------------------------------------
!     00 Jun 1983: Original version of lu1rec followed John Reid's
!                  compression routine in LA05.  It recovered
!                  space in ind(*) and optionally a(*)
!                  by eliminating entries with ind(l) = 0.
!                  The elements of ind(*) could not be negative.
!                  If len(i) was positive, entry i contained
!                  that many elements, starting at  loc(i).
!                  Otherwise, entry i was eliminated.

!     23 Mar 2001: Realised we could have len(i) = 0 in rare cases!
!                  (Mostly during TCP when the pivot row contains
!                  a column of length 1 that couldn't be a pivot.)
!                  Revised storage scheme to
!                     keep        entries with       ind(l) >  0,
!                     squeeze out entries with -n <= ind(l) <= 0,
!                  and to allow len(i) = 0.
!                  Empty items are moved to the end of the compressed
!                  ind(*) and/or a(*) arrays are given one empty space.
!                  Items with len(i) < 0 are still eliminated.

!     27 Mar 2001: Decided to use only ind(l) > 0 and = 0 in lu1fad.
!                  Still have to keep entries with len(i) = 0.

!     On exit:
!     ltop         is the length of useful entries in ind(*), a(*).
!     ind(ltop+1)  is "i" such that len(i), loc(i) belong to the last
!                  item in ind(*), a(*).
!     ------------------------------------------------------------------

    nempty = 0

    do 10 i = 1, n
        leni = len(i)
        if (leni > 0) then
            l      = loc(i) + leni - 1
            len(i) = ind(l)
            ind(l) = - (n + i)
        else if (leni == 0) then
            nempty = nempty + 1
        end if
    10 END DO

    k      = 0
    klast  = 0    ! Previous k
    ilast  = 0    ! Last entry moved.

    do 20 l = 1, ltop
        i    = ind(l)
        if (i > 0) then
            k      = k + 1
            ind(k) = i
            if (reals) a(k) = a(l)

        else if (i < -n) then

        !           This is the end of entry  i.

            i      = - (i + n)
            ilast  = i
            k      = k + 1
            ind(k) = len(i)
            if (reals) a(k) = a(l)
            loc(i) = klast + 1
            len(i) = k     - klast
            klast  = k
        end if
    20 END DO

!     Move any empty items to the end, adding 1 free entry for each.

    if (nempty > 0) then
        do i = 1, n
            if (len(i) == 0) then
                k      = k + 1
                loc(i) = k
                ind(k) = 0
                ilast  = i
            end if
        end do
    end if

    nout   = luparm(1)
    lprint = luparm(2)
    if (lprint >= 50) write(nout, 1000) ltop, k, reals, nempty
    luparm(26) = luparm(26) + 1  ! ncp

!     Return ilast in ind(ltop + 1).

    ltop        = k
    ind(ltop+1) = ilast
    return

    1000 format(' lu1rec.  File compressed from', i10, '   to', i10, l3, &
    '  nempty =', i8)

    end ! subroutine lu1rec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu1slk( m, n, lena, iq, iqloc, a, locc, w )

    implicit           none
    integer ::            m, n, lena
    integer ::            iq(n), iqloc(m), locc(n)
    double precision ::   a(lena), w(n)

!     ------------------------------------------------------------------
!     lu1slk  sets w(j) > 0 if column j is a unit vector.

!     21 Nov 2000: First version.  lu1fad needs it for TCP.
!                  Note that w(*) is nominally an integer array,
!                  but the only spare space is the double array w(*).
!     ------------------------------------------------------------------

    integer ::            j, lc1, lq, lq1, lq2

    do j = 1, n
        w(j) = 0.0d+0
    end do

    lq1    = iqloc(1)
    lq2    = n
    if (m > 1) lq2 = iqloc(2) - 1

    do lq = lq1, lq2
        j      = iq(lq)
        lc1    = locc(j)
        if (abs( a(lc1) ) == 1.0d+0) then
            w(j) = 1.0d+0
        end if
    end do

    end ! subroutine lu1slk

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu1ful( m     , n    , lena , lenD , lu1 , TPP, &
    mleft , nleft, nrank, nrowu, &
    lenL  , lenU , nsing, &
    keepLU, small, &
    a     , d    , indc , indr , ip  , iq, &
    lenc  , lenr , locc , ipinv, ipvt )

    implicit           double precision (a-h,o-z)
    logical ::            TPP       , keepLU
    double precision ::   a(lena)   , d(lenD)
    integer ::            indc(lena), indr(lena), ip(m)   , iq(n), &
    lenc(n)   , lenr(m)   , ipinv(m)
    integer ::            locc(n)   , ipvt(m)

!     ------------------------------------------------------------------
!     lu1ful computes a dense (full) LU factorization of the
!     mleft by nleft matrix that remains to be factored at the
!     beginning of the nrowu-th pass through the main loop of lu1fad.

!     02 May 1989: First version.
!     05 Feb 1994: Column interchanges added to lu1DPP.
!     08 Feb 1994: ipinv reconstructed, since lu1pq3 may alter ip.
!     ------------------------------------------------------------------

    parameter        ( zero = 0.0d+0 )


!     ------------------------------------------------------------------
!     If lu1pq3 moved any empty rows, reset ipinv = inverse of ip.
!     ------------------------------------------------------------------
    if (nrank < m) then
        do 100 l    = 1, m
            i        = ip(l)
            ipinv(i) = l
        100 END DO
    end if

!     ------------------------------------------------------------------
!     Copy the remaining matrix into the dense matrix D.
!     ------------------------------------------------------------------
!     call dload ( lenD, zero, d, 1 )
    do j = 1, lenD
        d(j) = zero
    end do
    ipbase = nrowu - 1
    ldbase = 1 - nrowu

    do 200 lq = nrowu, n
        j      = iq(lq)
        lc1    = locc(j)
        lc2    = lc1 + lenc(j) - 1

        do 150 lc = lc1, lc2
            i      = indc(lc)
            ld     = ldbase + ipinv(i)
            d(ld)  = a(lc)
        150 END DO

        ldbase = ldbase + mleft
    200 END DO

!     ------------------------------------------------------------------
!     Call our favorite dense LU factorizer.
!     ------------------------------------------------------------------
    if ( TPP ) then
        call lu1DPP( d, mleft, mleft, nleft, small, nsing, &
        ipvt, iq(nrowu) )
    else
        call lu1DCP( d, mleft, mleft, nleft, small, nsing, &
        ipvt, iq(nrowu) )
    end if

!     ------------------------------------------------------------------
!     Move D to the beginning of A,
!     and pack L and U at the top of a, indc, indr.
!     In the process, apply the row permutation to ip.
!     lkk points to the diagonal of U.
!     ------------------------------------------------------------------
    call dcopy ( lenD, d, 1, a, 1 )

    ldiagU = lena   - n
    lkk    = 1
    lkn    = lenD  - mleft + 1
    lu     = lu1

    do 450  k = 1, min( mleft, nleft )
        l1     = ipbase + k
        l2     = ipbase + ipvt(k)
        if (l1 /= l2) then
            i      = ip(l1)
            ip(l1) = ip(l2)
            ip(l2) = i
        end if
        ibest  = ip( l1 )
        jbest  = iq( l1 )

        if ( keepLU ) then
        !           ===========================================================
        !           Pack the next column of L.
        !           ===========================================================
            la     = lkk
            ll     = lu
            nrowd  = 1

            do 410  i = k + 1, mleft
                la     = la + 1
                ai     = a(la)
                if (abs( ai ) > small) then
                    nrowd    = nrowd + 1
                    ll       = ll    - 1
                    a(ll)    = ai
                    indc(ll) = ip( ipbase + i )
                    indr(ll) = ibest
                end if
            410 END DO

        !           ===========================================================
        !           Pack the next row of U.
        !           We go backwards through the row of D
        !           so the diagonal ends up at the front of the row of  U.
        !           Beware -- the diagonal may be zero.
        !           ===========================================================
            la     = lkn + mleft
            lu     = ll
            ncold  = 0

            do 420  j = nleft, k, -1
                la     = la - mleft
                aj     = a(la)
                if (abs( aj ) > small  .OR.  j == k) then
                    ncold    = ncold + 1
                    lu       = lu    - 1
                    a(lu)    = aj
                    indr(lu) = iq( ipbase + j )
                end if
            420 END DO

            lenr(ibest) = - ncold
            lenc(jbest) = - nrowd
            lenL        =   lenL + nrowd - 1
            lenU        =   lenU + ncold
            lkn         =   lkn  + 1

        else
        !           ===========================================================
        !           Store just the diagonal of U, in natural order.
        !           ===========================================================
            a(ldiagU + jbest) = a(lkk)
        end if

        lkk    = lkk  + mleft + 1
    450 END DO

    end ! subroutine lu1ful

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu1DPP( a, lda, m, n, small, nsing, &
    ipvt, iq )

    implicit           none
    integer ::            lda, m, n, nsing
    integer ::            ipvt(m), iq(n)
    double precision ::   a(lda,n), small

!     ------------------------------------------------------------------
!     lu1DPP factors a dense m x n matrix A by Gaussian elimination,
!     using row interchanges for stability, as in dgefa from LINPACK.
!     This version also uses column interchanges if all elements in a
!     pivot column are smaller than (or equal to) "small".  Such columns
!     are changed to zero and permuted to the right-hand end.

!     As in LINPACK, ipvt(*) keeps track of pivot rows.
!     Rows of U are interchanged, but we don't have to physically
!     permute rows of L.  In contrast, column interchanges are applied
!     directly to the columns of both L and U, and to the column
!     permutation vector iq(*).

!     02 May 1989: First version derived from dgefa
!                  in LINPACK (version dated 08/14/78).
!     05 Feb 1994: Generalized to treat rectangular matrices
!                  and use column interchanges when necessary.
!                  ipvt is retained, but column permutations are applied
!                  directly to iq(*).
!     21 Dec 1994: Bug found via example from Steve Dirkse.
!                  Loop 100 added to set ipvt(*) for singular rows.
!     26 Mar 2006: nsing redefined (see below).
!                  Changed to implicit none.
!     ------------------------------------------------------------------

!     On entry:

!        a       Array holding the matrix A to be factored.

!        lda     The leading dimension of the array  a.

!        m       The number of rows    in  A.

!        n       The number of columns in  A.

!        small   A drop tolerance.  Must be zero or positive.

!     On exit:

!        a       An upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        nsing   Number of singularities detected.
!                26 Mar 2006: nsing redefined to be more meaningful.
!                Users may define rankU = n - nsing and regard
!                U as upper-trapezoidal, with the first rankU columns
!                being triangular and the rest trapezoidal.
!                It would be better to return rankU, but we still
!                return nsing for compatibility (even though lu1fad
!                no longer uses it).

!        ipvt    Records the pivot rows.

!        iq      A vector to which column interchanges are applied.
!     ------------------------------------------------------------------

    double precision ::    t
    integer ::             idamax, i, j, k, kp1, l, last, lencol, rankU
    double precision ::    zero         ,  one
    parameter         ( zero = 0.0d+0,  one = 1.0d+0 )


    rankU  = 0
    k      = 1
    last   = n

!     ------------------------------------------------------------------
!     Start of elimination loop.
!     ------------------------------------------------------------------
    10 kp1    = k + 1
    lencol = m - k + 1

! Find l, the pivot row.

    l       = idamax( lencol, a(k,k), 1 ) + k - 1
    ipvt(k) = l

    if (abs( a(l,k) ) <= small) then
    !==============================================================
    ! Do column interchange, changing old pivot column to zero.
    ! Reduce "last" and try again with same k.
    !==============================================================
        j        = iq(last)
        iq(last) = iq(k)
        iq(k)    = j

        do i = 1, k - 1
            t         = a(i,last)
            a(i,last) = a(i,k)
            a(i,k)    = t
        end do

        do i = k, m
            t         = a(i,last)
            a(i,last) = zero
            a(i,k)    = t
        end do

        last     = last - 1
        if (k <= last) go to 10

    else
        rankU  = rankU + 1
        if (k < m) then
        !===========================================================
        ! o row interchange if necessary.
        !===========================================================
            if (l /= k) then
                t      = a(l,k)
                a(l,k) = a(k,k)
                a(k,k) = t
            end if

        !===========================================================
        ! Compute multipliers.
        ! Do row elimination with column indexing.
        !===========================================================
            t = - one / a(k,k)
            call dscal ( m-k, t, a(kp1,k), 1 )

            do j = kp1, last
                t    = a(l,j)
                if (l /= k) then
                    a(l,j) = a(k,j)
                    a(k,j) = t
                end if
                call daxpy ( m-k, t, a(kp1,k), 1, a(kp1,j), 1 )
            end do

            k = k + 1
            if (k <= last) go to 10
        end if
    end if

! Set ipvt(*) for singular rows.

    do 100 k = last + 1, m
        ipvt(k) = k
    100 END DO

    nsing  = n - rankU

    end ! subroutine lu1DPP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu1DCP( a, lda, m, n, small, nsing, &
    ipvt, iq )

    implicit           none
    integer ::            lda, m, n, nsing
    integer ::            ipvt(m), iq(n)
    double precision ::   a(lda,n), small

!     ------------------------------------------------------------------
!     lu1DCP factors a dense m x n matrix A by Gaussian elimination,
!     using Complete Pivoting (row and column interchanges) for
!     stability.
!     This version also uses column interchanges if all elements in a
!     pivot column are smaller than (or equal to) "small".  Such columns
!     are changed to zero and permuted to the right-hand end.

!     As in LINPACK's dgefa, ipvt(*) keeps track of pivot rows.
!     Rows of U are interchanged, but we don't have to physically
!     permute rows of L.  In contrast, column interchanges are applied
!     directly to the columns of both L and U, and to the column
!     permutation vector iq(*).

!     01 May 2002: First dense Complete Pivoting, derived from lu1DPP.
!     07 May 2002: Another break needed at end of first loop.
!     26 Mar 2006: Cosmetic mods while looking for "nsing" bug when m<n.
!                  nsing redefined (see below).
!                  Changed to implicit none.
!     ------------------------------------------------------------------

!     On entry:

!        a       Array holding the matrix A to be factored.

!        lda     The leading dimension of the array  a.

!        m       The number of rows    in  A.

!        n       The number of columns in  A.

!        small   A drop tolerance.  Must be zero or positive.

!     On exit:

!        a       An upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        nsing   Number of singularities detected.
!                26 Mar 2006: nsing redefined to be more meaningful.
!                Users may define rankU = n - nsing and regard
!                U as upper-trapezoidal, with the first rankU columns
!                being triangular and the rest trapezoidal.
!                It would be better to return rankU, but we still
!                return nsing for compatibility (even though lu1fad
!                no longer uses it).

!        ipvt    Records the pivot rows.

!        iq      A vector to which column interchanges are applied.
!     ------------------------------------------------------------------

    double precision ::    aijmax, ajmax, t
    integer ::             idamax, i, imax, j, jlast, jmax, jnew, &
    k, kp1, l, last, lencol, rankU
    double precision ::    zero         ,  one
    parameter         ( zero = 0.0d+0,  one = 1.0d+0 )


    rankU  = 0
    lencol = m + 1
    last   = n

!-----------------------------------------------------------------
! Start of elimination loop.
!-----------------------------------------------------------------
    do k = 1, n
        kp1    = k + 1
        lencol = lencol - 1

    ! Find the biggest aij in row imax and column jmax.

        aijmax = zero
        imax   = k
        jmax   = k
        jlast  = last

        do j = k, jlast
            10 l      = idamax( lencol, a(k,j), 1 ) + k - 1
            ajmax  = abs( a(l,j) )

            if (ajmax <= small) then
            !========================================================
            ! Do column interchange, changing old column to zero.
            ! Reduce  "last"  and try again with same j.
            !========================================================
                jnew     = iq(last)
                iq(last) = iq(j)
                iq(j)    = jnew

                do i = 1, k - 1
                    t         = a(i,last)
                    a(i,last) = a(i,j)
                    a(i,j)    = t
                end do

                do i = k, m
                    t         = a(i,last)
                    a(i,last) = zero
                    a(i,j)    = t
                end do

                last   = last - 1
                if (j <= last) go to 10 ! repeat
                go to 200                 ! break
            end if

        ! Check if this column has biggest aij so far.

            if (aijmax < ajmax) then
                aijmax  =   ajmax
                imax    =   l
                jmax    =   j
            end if

            if (j >= last) go to 200   ! break
        end do

        200 ipvt(k) = imax

        if (jmax /= k) then
        !==========================================================
        ! Do column interchange (k and jmax).
        !==========================================================
            jnew     = iq(jmax)
            iq(jmax) = iq(k)
            iq(k)    = jnew

            do i = 1, m
                t         = a(i,jmax)
                a(i,jmax) = a(i,k)
                a(i,k)    = t
            end do
        end if

        if (k < m) then
        !===========================================================
        ! Do row interchange if necessary.
        !===========================================================
            t         = a(imax,k)
            if (imax /= k) then
                a(imax,k) = a(k,k)
                a(k,k)    = t
            end if

        !===========================================================
        ! Compute multipliers.
        ! Do row elimination with column indexing.
        !===========================================================
            t      = - one / t
            call dscal ( m-k, t, a(kp1,k), 1 )

            do j = kp1, last
                t         = a(imax,j)
                if (imax /= k) then
                    a(imax,j) = a(k,j)
                    a(k,j)    = t
                end if
                call daxpy ( m-k, t, a(kp1,k), 1, a(kp1,j), 1 )
            end do

        else
            go to 500               ! break
        end if

        if (k >= last) go to 500 ! break
    end do

! Set ipvt(*) for singular rows.

    500 do k = last + 1, m
        ipvt(k) = k
    end do

    nsing  = n - rankU

    end ! subroutine lu1DCP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     File  lusol2.f

!     Hbuild   Hchange  Hdelete  Hdown    Hinsert  Hup

!     Heap-management routines for LUSOL's lu1fac.
!     May be useful for other applications.

! 11 Feb 2002: MATLAB version derived from "Algorithms" by R. Sedgewick.
! 03 Mar 2002: F77    version derived from MATLAB version.
! 07 May 2002: Safeguard input parameters k, N, Nk.
!              We don't want them to be output!
! 19 Dec 2004: Hdelete: Nin is new input parameter for length of Hj, Ha.
! 19 Dec 2004: Current version of lusol2.f.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     For LUSOL, the heap structure involves three arrays of length N.
!     N        is the current number of entries in the heap.
!     Ha(1:N)  contains the values that the heap is partially sorting.
!              For LUSOL they are double precision values -- the largest
!              element in each remaining column of the updated matrix.
!              The biggest entry is in Ha(1), the top of the heap.
!     Hj(1:N)  contains column numbers j.
!              Ha(k) is the biggest entry in column j = Hj(k).
!     Hk(1:N)  contains indices within the heap.  It is the
!              inverse of Hj(1:N), so  k = Hk(j)  <=>  j = Hj(k).
!              Column j is entry k in the heap.
!     hops     is the number of heap operations,
!              i.e., the number of times an entry is moved
!              (the number of "hops" up or down the heap).
!     Together, Hj and Hk let us find values inside the heap
!     whenever we want to change one of the values in Ha.
!     For other applications, Ha may need to be some other data type,
!     like the keys that sort routines operate on.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hbuild( Ha, Hj, Hk, N, Nk, hops )

    implicit &
    none
    integer :: &
    N, Nk, hops, Hj(N), Hk(Nk)
    double precision :: &
    Ha(N)

!     ==================================================================
!     Hbuild initializes the heap by inserting each element of Ha.
!     Input:  Ha, Hj.
!     Output: Ha, Hj, Hk, hops.

!     01 May 2002: Use k for new length of heap, not k-1 for old length.
!     05 May 2002: Use kk in call to stop loop variable k being altered.
!                  (Actually Hinsert no longer alters that parameter.)
!     07 May 2002: ftnchek wants us to protect Nk, Ha(k), Hj(k) too.
!     07 May 2002: Current version of Hbuild.
!     ==================================================================

    integer :: &
    h, jv, k, kk, Nkk
    double precision :: &
    v

    Nkk  = Nk
    hops = 0
    do k = 1, N
        kk    = k
        v     = Ha(k)
        jv    = Hj(k)
        call Hinsert( Ha, Hj, Hk, kk, Nkk, v, jv, h )
        hops  = hops + h
    end do

    end ! subroutine Hbuild

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hchange( Ha, Hj, Hk, N, Nk, k, v, jv, hops )

    implicit &
    none
    integer :: &
    N, Nk, k, jv, hops, Hj(N), Hk(Nk)
    double precision :: &
    v, Ha(N)

!     ==================================================================
!     Hchange changes Ha(k) to v in heap of length N.

!     01 May 2002: Need Nk for length of Hk.
!     07 May 2002: Protect input parameters N, Nk, k.
!     07 May 2002: Current version of Hchange.
!     ==================================================================

    integer :: &
    kx, Nx, Nkx
    double precision :: &
    v1

    Nx     = N
    Nkx    = Nk
    kx     = k
    v1     = Ha(k)
    Ha(k)  = v
    Hj(k)  = jv
    Hk(jv) = k
    if (v1 < v) then
        call Hup   ( Ha, Hj, Hk, Nx, Nkx, kx, hops )
    else
        call Hdown ( Ha, Hj, Hk, Nx, Nkx, kx, hops )
    end if

    end ! subroutine Hchange

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hdelete( Ha, Hj, Hk, Nin, N, Nk, k, hops )

    implicit &
    none
    integer :: &
    N, Nin, Nk, k, hops, Hj(Nin), Hk(Nk)
    double precision :: &
    Ha(Nin)

!     ==================================================================
!     Hdelete deletes Ha(k) from heap of length N.

!     03 Apr 2002: Current version of Hdelete.
!     01 May 2002: Need Nk for length of Hk.
!     07 May 2002: Protect input parameters N, Nk, k.
!     19 Dec 2004: Nin is new input parameter for length of Hj, Ha.
!     19 Dec 2004: Current version of Hdelete.
!     ==================================================================

    integer :: &
    jv, kx, Nkx, Nx
    double precision :: &
    v

    kx    = k
    Nkx   = Nk
    Nx    = N
    v     = Ha(N)
    jv    = Hj(N)
    N     = N - 1
    hops  = 0
    if (k <= N) then
        call Hchange( Ha, Hj, Hk, Nx, Nkx, kx, v, jv, hops )
    end if

    end ! subroutine Hdelete

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hdown ( Ha, Hj, Hk, N, Nk, kk, hops )

    implicit &
    none
    integer :: &
    N, Nk, kk, hops, Hj(N), Hk(Nk)
    double precision :: &
    Ha(N)

!     ==================================================================
!     Hdown  updates heap by moving down tree from node k.

!     01 May 2002: Need Nk for length of Hk.
!     05 May 2002: Change input paramter k to kk to stop k being output.
!     05 May 2002: Current version of Hdown.
!     ==================================================================

    integer :: &
    j, jj, jv, k, N2
    double precision :: &
    v

    k     = kk
    hops  = 0
    v     = Ha(k)
    jv    = Hj(k)
    N2    = N/2

!     while 1
    100 if (k > N2   ) go to 200   ! break
    hops   = hops + 1
    j      = k+k
    if (j < N) then
        if (Ha(j) < Ha(j+1)) j = j+1
    end if
    if (v >= Ha(j)) go to 200   ! break
    Ha(k)  = Ha(j)
    jj     = Hj(j)
    Hj(k)  = jj
    Hk(jj) =  k
    k      =  j
    go to 100
!     end while

    200 Ha(k)  =  v
    Hj(k)  = jv
    Hk(jv) =  k

    end ! subroutine Hdown

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hinsert( Ha, Hj, Hk, N, Nk, v, jv, hops )

    implicit &
    none
    integer :: &
    N, Nk, jv, hops, Hj(N), Hk(Nk)
    double precision :: &
    v, Ha(N)

!     ==================================================================
!     Hinsert inserts (v,jv) into heap of length N-1
!     to make heap of length N.

!     03 Apr 2002: First version of Hinsert.
!     01 May 2002: Require N to be final length, not old length.
!                  Need Nk for length of Hk.
!     07 May 2002: Protect input parameters N, Nk.
!     07 May 2002: Current version of Hinsert.
!     ==================================================================

    integer :: &
    kk, Nkk, Nnew

    Nnew     = N
    Nkk      = Nk
    kk       = Nnew
    Ha(Nnew) =  v
    Hj(Nnew) = jv
    Hk(jv)   = Nnew
    call Hup   ( Ha, Hj, Hk, Nnew, Nkk, kk, hops )

    end ! subroutine Hinsert

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Hup   ( Ha, Hj, Hk, N, Nk, kk, hops )

    implicit &
    none
    integer :: &
    N, Nk, kk, hops, Hj(N), Hk(Nk)
    double precision :: &
    Ha(N)

!     ==================================================================
!     Hup updates heap by moving up tree from node k.

!     01 May 2002: Need Nk for length of Hk.
!     05 May 2002: Change input paramter k to kk to stop k being output.
!     05 May 2002: Current version of Hup.
!     ==================================================================

    integer :: &
    j, jv, k, k2
    double precision :: &
    v

    k     = kk
    hops  = 0
    v     = Ha(k)
    jv    = Hj(k)
!     while 1
    100 if (k <  2    ) go to 200   ! break
    k2    = k/2
    if (v < Ha(k2)) go to 200   ! break
    hops  = hops + 1
    Ha(k) = Ha(k2)
    j     = Hj(k2)
    Hj(k) =  j
    Hk(j) =  k
    k     = k2
    go to 100
!     end while

    200 Ha(k)  =  v
    Hj(k)  = jv
    Hk(jv) =  k

    end ! subroutine Hup
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     File  lusol6a.f

!     lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD
!     lu6chk

! 26 Apr 2002: lu6 routines put into a separate file.
! 15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut.
!              lu6LD implemented to allow solves with LDL' or L|D|L'.
! 23 Apr 2004: lu6chk modified.  TRP can judge singularity better
!              by comparing all diagonals to DUmax.
! 27 Jun 2004: lu6chk.  Allow write only if nout .gt. 0 .
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6sol( mode, m, n, v, w, &
    lena, luparm, parmlu, &
    a, indc, indr, ip, iq, &
    lenc, lenr, locc, locr, &
    inform )

    implicit &
    none
    integer :: &
    luparm(30), mode, m, n, lena, inform, &
    indc(lena), indr(lena), ip(m), iq(n), &
    lenc(n), lenr(m), locc(n), locr(m)
    double precision :: &
    parmlu(30), a(lena), v(m), w(n)

!-----------------------------------------------------------------------
!     lu6sol  uses the factorization  A = L U  as follows:

!     mode
!      1    v  solves   L v = v(input).   w  is not touched.
!      2    v  solves   L'v = v(input).   w  is not touched.
!      3    w  solves   U w = v.          v  is not altered.
!      4    v  solves   U'v = w.          w  is destroyed.
!      5    w  solves   A w = v.          v  is altered as in 1.
!      6    v  solves   A'v = w.          w  is destroyed.

!     If mode = 3,4,5,6, v and w must not be the same arrays.

!     If lu1fac has just been used to factorize a symmetric matrix A
!     (which must be definite or quasi-definite), the factors A = L U
!     may be regarded as A = LDL', where D = diag(U).  In such cases,

!     mode
!      7    v  solves   A v = L D L'v = v(input).   w  is not touched.
!      8    v  solves       L |D| L'v = v(input).   w  is not touched.

!     ip(*), iq(*)      hold row and column numbers in pivotal order.
!     lenc(k)           is the length of the k-th column of initial L.
!     lenr(i)           is the length of the i-th row of U.
!     locc(*)           is not used.
!     locr(i)           is the start  of the i-th row of U.

!     U is assumed to be in upper-trapezoidal form (nrank by n).
!     The first entry for each row is the diagonal element
!     (according to the permutations  ip, iq).  It is stored at
!     location locr(i) in a(*), indr(*).

!     On exit, inform = 0 except as follows.
!     If mode = 3,4,5,6 and if U (and hence A) is singular, then
!     inform = 1 if there is a nonzero residual in solving the system
!     involving U.  parmlu(20) returns the norm of the residual.

!       July 1987: Early version.
!     09 May 1988: f77 version.
!     27 Apr 2000: Abolished the dreaded "computed go to".
!                  But hard to change other "go to"s to "if then else".
!     15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol.
!-----------------------------------------------------------------------

    if      (mode == 1) then             ! Solve  L v(new) = v.
        call lu6L  ( &
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )

    else if (mode == 2) then             ! Solve  L'v(new) = v.
        call lu6Lt ( &
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )

    else if (mode == 3) then             ! Solve  U w = v.
        call lu6U  ( &
        inform, m, n, v, w, &
        lena, luparm, parmlu, &
        a, indr, ip, iq, lenr, locr )

    else if (mode == 4) then             ! Solve  U'v = w.
        call lu6Ut ( &
        inform, m, n, v, w, &
        lena, luparm, parmlu, &
        a, indr, ip, iq, lenr, locr )

    else if (mode == 5) then             ! Solve  A w      = v
        call lu6L  (                         & ! via    L v(new) = v
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )
        call lu6U  (                         & ! and    U w = v(new).
        inform, m, n, v, w, &
        lena, luparm, parmlu, &
        a, indr, ip, iq, lenr, locr )

    else if (mode == 6) then             ! Solve  A'v = w
        call lu6Ut (                         & ! via    U'v = w
        inform, m, n, v, w, &
        lena, luparm, parmlu, &
        a, indr, ip, iq, lenr, locr )
        call lu6Lt (                         & ! and    L'v(new) = v.
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )

    else if (mode == 7) then
        call lu6LD (                         & ! Solve  LDv(bar) = v
        inform, 1, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc, locr )
        call lu6Lt (                         & ! and    L'v(new) = v(bar).
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )

    else if (mode == 8) then
        call lu6LD (                         & ! Solve  L|D|v(bar) = v
        inform, 2, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc, locr )
        call lu6Lt (                         & ! and    L'v(new) = v(bar).
        inform, m, n, v, &
        lena, luparm, parmlu, &
        a, indc, indr, lenc )
    end if

    end ! subroutine lu6sol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6L  ( &
    inform, m, n, v, &
    lena, luparm, parmlu, &
    a, indc, indr, lenc )

    implicit &
    none
    integer :: &
    inform, m, n, lena, luparm(30), &
    indc(lena), indr(lena), lenc(n)
    double precision :: &
    parmlu(30), a(lena), v(m)

!     ------------------------------------------------------------------
!     lu6L   solves   L v = v(input).

!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

    integer :: &
    i, ipiv, j, k, l, l1, ldummy, len, lenL, lenL0, numL, numL0
    double precision :: &
    small, vpiv

    numL0  = luparm(20)
    lenL0  = luparm(21)
    lenL   = luparm(23)
    small  = parmlu(3)
    inform = 0
    l1     = lena + 1

    do k = 1, numL0
        len   = lenc(k)
        l     = l1
        l1    = l1 - len
        ipiv  = indr(l1)
        vpiv  = v(ipiv)

        if (abs( vpiv ) > small) then
        !***** This loop could be coded specially.
            do ldummy = 1, len
                l    = l - 1
                j    = indc(l)
                v(j) = v(j)  +  a(l) * vpiv
            end do
        end if
    end do

    l      = lena - lenL0 + 1
    numL   = lenL - lenL0

!***** This loop could be coded specially.

    do ldummy = 1, numL
        l      = l - 1
        i      = indr(l)
        if (abs( v(i) ) > small) then
            j    = indc(l)
            v(j) = v(j)  +  a(l) * v(i)
        end if
    end do

!     Exit.

    luparm(10) = inform

    end ! subroutine lu6L

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6Lt ( &
    inform, m, n, v, &
    lena, luparm, parmlu, &
    a, indc, indr, lenc )

    implicit &
    none
    integer :: &
    inform, m, n, lena, luparm(30), &
    indc(lena), indr(lena), lenc(n)
    double precision :: &
    parmlu(30), a(lena), v(m)

!     ------------------------------------------------------------------
!     lu6Lt  solves   L'v = v(input).

!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

    integer :: &
    i, ipiv, j, k, l, l1, l2, len, lenL, lenL0, numL0
    double precision :: &
    small, sum

!     ------------------------------------------------------------------
    double precision ::   zero
    parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

    numL0  = luparm(20)
    lenL0  = luparm(21)
    lenL   = luparm(23)
    small  = parmlu(3)
    inform = 0
    l1     = lena - lenL + 1
    l2     = lena - lenL0

!***** This loop could be coded specially.
    do l = l1, l2
        j     = indc(l)
        if (abs( v(j) ) > small) then
            i     = indr(l)
            v(i)  = v(i)  +  a(l) * v(j)
        end if
    end do

    do k = numL0, 1, -1
        len   = lenc(k)
        sum   = zero
        l1    = l2 + 1
        l2    = l2 + len

    !***** This loop could be coded specially.
        do l = l1, l2
            j     = indc(l)
            sum   = sum  +  a(l) * v(j)
        end do

        ipiv    = indr(l1)
        v(ipiv) = v(ipiv) + sum
    end do

!     Exit.

    luparm(10) = inform

    end ! subroutine lu6Lt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6U  ( &
    inform, m, n, v, w, &
    lena, luparm, parmlu, &
    a, indr, ip, iq, lenr, locr )

    implicit &
    none
    integer :: &
    inform, m, n, lena, luparm(30), &
    indr(lena), ip(m), iq(n), lenr(m), locr(m)
    double precision :: &
    parmlu(30), a(lena), v(m), w(n)

!     ------------------------------------------------------------------
!     lu6U   solves   U w = v.          v  is not altered.

!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

    integer :: &
    i, j, k, klast, l, l1, l2, l3, nrank, nrank1
    double precision :: &
    resid, small, t

!     ------------------------------------------------------------------
    double precision ::   zero
    parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

    nrank  = luparm(16)
    small  = parmlu(3)
    inform = 0
    nrank1 = nrank + 1
    resid  = zero

!     Find the first nonzero in v(1:nrank), counting backwards.

    do klast = nrank, 1, -1
        i      = ip(klast)
        if (abs( v(i) ) > small) go to 320
    end do

    320 do k = klast + 1, n
        j     = iq(k)
        w(j)  = zero
    end do

!     Do the back-substitution, using rows 1:klast of U.

    do k  = klast, 1, -1
        i      = ip(k)
        t      = v(i)
        l1     = locr(i)
        l2     = l1 + 1
        l3     = l1 + lenr(i) - 1

    !***** This loop could be coded specially.
        do l = l2, l3
            j     = indr(l)
            t     = t  -  a(l) * w(j)
        end do

        j      = iq(k)
        if (abs( t ) <= small) then
            w(j)  = zero
        else
            w(j)  = t / a(l1)
        end if
    end do

!     Compute residual for overdetermined systems.

    do k = nrank1, m
        i     = ip(k)
        resid = resid  +  abs( v(i) )
    end do

!     Exit.

    if (resid > zero) inform = 1
    luparm(10) = inform
    parmlu(20) = resid

    end ! subroutine lu6U

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6Ut ( &
    inform, m, n, v, w, &
    lena, luparm, parmlu, &
    a, indr, ip, iq, lenr, locr )

    implicit &
    none
    integer :: &
    inform, m, n, lena, luparm(30), &
    indr(lena), ip(m), iq(n), lenr(m), locr(m)
    double precision :: &
    parmlu(30), a(lena), v(m), w(n)

!     ------------------------------------------------------------------
!     lu6Ut  solves   U'v = w.          w  is destroyed.

!     15 Dec 2002: First version derived from lu6sol.
!     15 Dec 2002: Current version.
!     ------------------------------------------------------------------

    integer :: &
    i, j, k, l, l1, l2, nrank, nrank1
    double precision :: &
    resid, small, t

!     ------------------------------------------------------------------
    double precision ::   zero
    parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------

    nrank  = luparm(16)
    small  = parmlu(3)
    inform = 0
    nrank1 = nrank + 1
    resid  = zero

    do k = nrank1, m
        i     = ip(k)
        v(i)  = zero
    end do

!     Do the forward-substitution, skipping columns of U(transpose)
!     when the associated element of w(*) is negligible.

    do 480 k = 1, nrank
        i      = ip(k)
        j      = iq(k)
        t      = w(j)
        if (abs( t ) <= small) then
            v(i) = zero
            go to 480
        end if

        l1     = locr(i)
        t      = t / a(l1)
        v(i)   = t
        l2     = l1 + lenr(i) - 1
        l1     = l1 + 1

    !***** This loop could be coded specially.
        do l = l1, l2
            j     = indr(l)
            w(j)  = w(j)  -  t * a(l)
        end do
    480 END DO

!     Compute residual for overdetermined systems.

    do k = nrank1, n
        j     = iq(k)
        resid = resid  +  abs( w(j) )
    end do

!     Exit.

    if (resid > zero) inform = 1
    luparm(10) = inform
    parmlu(20) = resid

    end ! subroutine lu6Ut

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6LD ( &
    inform, mode, m, n, v, &
    lena, luparm, parmlu, &
    a, indc, indr, lenc, locr )

    implicit &
    none
    integer :: &
    luparm(30), inform, mode, m, n, lena, &
    indc(lena), indr(lena), lenc(n), locr(m)
    double precision :: &
    parmlu(30), a(lena), v(m)

!-----------------------------------------------------------------------
!     lu6LD  assumes lu1fac has computed factors A = LU of a
!     symmetric definite or quasi-definite matrix A,
!     using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3,
!     or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4.
!     It also assumes that no updates have been performed.
!     In such cases,  U = D L', where D = diag(U).
!     lu6LDL returns v as follows:

!     mode
!      1    v  solves   L D v = v(input).
!      2    v  solves   L|D|v = v(input).

!     15 Dec 2002: First version of lu6LD.
!     15 Dec 2002: Current version.
!-----------------------------------------------------------------------

! Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode.
! The code for L is the same as in lu6L,
! but when a nonzero entry of v arises, we divide by
! the corresponding entry of D or |D|.

    integer :: &
    ipiv, j, k, l, l1, ldummy, len, numL0
    double precision :: &
    diag, small, vpiv

    numL0  = luparm(20)
    small  = parmlu(3)
    inform = 0
    l1     = lena + 1

    do k = 1, numL0
        len   = lenc(k)
        l     = l1
        l1    = l1 - len
        ipiv  = indr(l1)
        vpiv  = v(ipiv)

        if (abs( vpiv ) > small) then
        !***** This loop could be coded specially.
            do ldummy = 1, len
                l    = l - 1
                j    = indc(l)
                v(j) = v(j)  +  a(l) * vpiv
            end do

        ! Find diag = U(ipiv,ipiv) and divide by diag or |diag|.

            l    = locr(ipiv)
            diag = A(l)
            if (mode == 2) diag = abs( diag )
            v(ipiv) = vpiv / diag
        end if
    end do

    end ! subroutine lu6LD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu6chk( mode, m, n, w, &
    lena, luparm, parmlu, &
    a, indc, indr, ip, iq, &
    lenc, lenr, locc, locr, &
    inform )

    implicit &
    none
    integer :: &
    mode, m, n, lena, inform, &
    luparm(30), indc(lena), indr(lena), ip(m), iq(n), &
    lenc(n), lenr(m), locc(n), locr(m)
    double precision :: &
    parmlu(30), a(lena), w(n)

!     ------------------------------------------------------------------
!     lu6chk  looks at the LU factorization  A = L*U.

!     If mode = 1, lu6chk is being called by lu1fac.
!     (Other modes not yet implemented.)
!     The important input parameters are

!                    lprint = luparm(2)
!                             luparm(6) = 1 if TRP
!                    keepLU = luparm(8)
!                    Utol1  = parmlu(4)
!                    Utol2  = parmlu(5)

!     and the significant output parameters are

!                    inform = luparm(10)
!                    nsing  = luparm(11)
!                    jsing  = luparm(12)
!                    jumin  = luparm(19)
!                    Lmax   = parmlu(11)
!                    Umax   = parmlu(12)
!                    DUmax  = parmlu(13)
!                    DUmin  = parmlu(14)
!                    and      w(*).

!     Lmax  and Umax  return the largest elements in L and U.
!     DUmax and DUmin return the largest and smallest diagonals of U
!                     (excluding diagonals that are exactly zero).

!     In general, w(j) is set to the maximum absolute element in
!     the j-th column of U.  However, if the corresponding diagonal
!     of U is small in absolute terms or relative to w(j)
!     (as judged by the parameters Utol1, Utol2 respectively),
!     then w(j) is changed to - w(j).

!     Thus, if w(j) is not positive, the j-th column of A
!     appears to be dependent on the other columns of A.
!     The number of such columns, and the position of the last one,
!     are returned as nsing and jsing.

!     Note that nrank is assumed to be set already, and is not altered.
!     Typically, nsing will satisfy      nrank + nsing = n,  but if
!     Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur.

!     If keepLU = 0,
!     Lmax  and Umax  are already set by lu1fac.
!     The diagonals of U are in the top of A.
!     Only Utol1 is used in the singularity test to set w(*).

!     inform = 0  if  A  appears to have full column rank  (nsing = 0).
!     inform = 1  otherwise  (nsing .gt. 0).

!     00 Jul 1987: Early version.
!     09 May 1988: f77 version.
!     11 Mar 2001: Allow for keepLU = 0.
!     17 Nov 2001: Briefer output for singular factors.
!     05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom).
!     06 May 2002: With keepLU = 0, diags of U are in natural order.
!                  They were not being extracted correctly.
!     23 Apr 2004: TRP can judge singularity better by comparing
!                  all diagonals to DUmax.
!     27 Jun 2004: (PEG) Allow write only if nout .gt. 0 .
!     ------------------------------------------------------------------

    character &
    mnkey
    logical :: &
    keepLU, TRP
    integer :: &
    i, j, jsing, jumin, k, l, l1, l2, ldiagU, lenL, lprint, &
    ndefic, nout, nrank, nsing
    double precision :: &
    aij, diag, DUmax, DUmin, Lmax, Umax, Utol1, Utol2

    double precision ::   zero
    parameter        ( zero = 0.0d+0 )

    nout   = luparm(1)
    lprint = luparm(2)
    TRP    = luparm(6) == 1  ! Threshold Rook Pivoting
    keepLU = luparm(8) /= 0
    nrank  = luparm(16)
    lenL   = luparm(23)
    Utol1  = parmlu(4)
    Utol2  = parmlu(5)

    inform = 0
    Lmax   = zero
    Umax   = zero
    nsing  = 0
    jsing  = 0
    jumin  = 0
    DUmax  = zero
    DUmin  = 1.0d+30

    do j = 1, n
        w(j) = zero
    end do


    if (keepLU) then
    !--------------------------------------------------------------
    ! Find  Lmax.
    !--------------------------------------------------------------
        do l = lena + 1 - lenL, lena
            Lmax  = max( Lmax, abs(a(l)) )
        end do

    !--------------------------------------------------------------
    ! Find Umax and set w(j) = maximum element in j-th column of U.
    !--------------------------------------------------------------
        do k = 1, nrank
            i     = ip(k)
            l1    = locr(i)
            l2    = l1 + lenr(i) - 1

            do l = l1, l2
                j     = indr(l)
                aij   = abs( a(l) )
                w(j)  = max( w(j), aij )
                Umax  = max( Umax, aij )
            end do
        end do

        parmlu(11) = Lmax
        parmlu(12) = Umax

    !--------------------------------------------------------------
    ! Find DUmax and DUmin, the extreme diagonals of U.
    !--------------------------------------------------------------
        do k = 1, nrank
            j      = iq(k)
            i      = ip(k)
            l1     = locr(i)
            diag   = abs( a(l1) )
            DUmax  = max( DUmax, diag )
            if (DUmin > diag) then
                DUmin  =   diag
                jumin  =   j
            end if
        end do

    else
    !--------------------------------------------------------------
    ! keepLU = 0.
    ! Only diag(U) is stored.  Set w(*) accordingly.
    ! Find DUmax and DUmin, the extreme diagonals of U.
    !--------------------------------------------------------------
        ldiagU = lena - n

        do k = 1, nrank
            j      = iq(k)
        ! diag   = abs( a(ldiagU + k) ) ! 06 May 2002: Diags
            diag   = abs( a(ldiagU + j) ) ! are in natural order
            w(j)   = diag
            DUmax  = max( DUmax, diag )
            if (DUmin > diag) then
                DUmin  =   diag
                jumin  =   j
            end if
        end do
    end if


!--------------------------------------------------------------
! Negate w(j) if the corresponding diagonal of U is
! too small in absolute terms or relative to the other elements
! in the same column of  U.

! 23 Apr 2004: TRP ensures that diags are NOT small relative to
!              other elements in their own column.
!              Much better, we can compare all diags to DUmax.
!--------------------------------------------------------------
    if (mode == 1  .AND.  TRP) then
        Utol1 = max( Utol1, Utol2*DUmax )
    end if

    if (keepLU) then
        do k = 1, n
            j     = iq(k)
            if (k > nrank) then
                diag   = zero
            else
                i      = ip(k)
                l1     = locr(i)
                diag   = abs( a(l1) )
            end if

            if (diag <= Utol1  .OR.  diag <= Utol2*w(j)) then
                nsing  =   nsing + 1
                jsing  =   j
                w(j)   = - w(j)
            end if
        end do

    else ! keepLU = 0

        do k = 1, n
            j      = iq(k)
            diag   = w(j)

            if (diag <= Utol1) then
                nsing  =   nsing + 1
                jsing  =   j
                w(j)   = - w(j)
            end if
        end do
    end if


!-----------------------------------------------------------------
! Set output parameters.
!-----------------------------------------------------------------
    if (jumin == 0) DUmin = zero
    luparm(11) = nsing
    luparm(12) = jsing
    luparm(19) = jumin
    parmlu(13) = DUmax
    parmlu(14) = DUmin

    if (nsing > 0) then  ! The matrix has been judged singular.
        inform = 1
        ndefic = n - nrank
        if (nout > 0  .AND.  lprint >= 0) then
            if (m > n) then
                mnkey  = '>'
            else if (m == n) then
                mnkey  = '='
            else
                mnkey  = '<'
            end if
            write(nout, 1100) mnkey, nrank, ndefic, nsing
        end if
    end if

!     Exit.

    luparm(10) = inform
    return

    1100 format(' Singular(m', a, 'n)', &
    '  rank', i9, '  n-rank', i8, '  nsing', i9)

    end ! subroutine lu6chk

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!***********************************************************************

!     File  lusol7a.f

!     lu7add   lu7cyc   lu7elm   lu7for   lu7rnk   lu7zap

!     Utilities for LUSOL's update routines.
!     lu7for is the most important -- the forward sweep.

! 01 May 2002: Derived from LUSOL's original lu7a.f file.
! 01 May 2002: Current version of lusol7a.f.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7add( m, n, jadd, v, &
    lena, luparm, parmlu, &
    lenL, lenU, lrow, nrank, &
    a, indr, ip, lenr, locr, &
    inform, klast, vnorm )

    implicit           double precision (a-h,o-z)
    integer ::            luparm(30)
    double precision ::   parmlu(30), a(lena), v(m)
    integer ::            indr(lena), ip(m), lenr(m)
    integer ::            locr(m)

!     ------------------------------------------------------------------
!     lu7add  inserts the first nrank elements of the vector v(*)
!     as column  jadd  of  U.  We assume that  U  does not yet have any
!     entries in this column.
!     Elements no larger than  parmlu(3)  are treated as zero.
!     klast  will be set so that the last row to be affected
!     (in pivotal order) is row  ip(klast).

!     09 May 1988: First f77 version.
!     ------------------------------------------------------------------

    parameter        ( zero = 0.0d+0 )

    small  = parmlu(3)
    vnorm  = zero
    klast  = 0

    do 200 k  = 1, nrank
        i      = ip(k)
        if (abs( v(i) ) <= small) go to 200
        klast  = k
        vnorm  = vnorm  +  abs( v(i) )
        leni   = lenr(i)

    !        Compress row file if necessary.

        minfre = leni + 1
        nfree  = lena - lenL - lrow
        if (nfree < minfre) then
            call lu1rec( m, .TRUE. , luparm, lrow, lena, &
            a, indr, lenr, locr )
            nfree  = lena - lenL - lrow
            if (nfree < minfre) go to 970
        end if

    !        Move row  i  to the end of the row file,
    !        unless it is already there.
    !        No need to move if there is a gap already.

        if (leni == 0) locr(i) = lrow + 1
        lr1    = locr(i)
        lr2    = lr1 + leni - 1
        if (lr2    ==   lrow) go to 150
        if (indr(lr2+1) == 0) go to 180
        locr(i) = lrow + 1

        do 140 l = lr1, lr2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
        140 END DO

        150 lr2     = lrow
        lrow    = lrow + 1

    !        Add the element of  v.

        180 lr2       = lr2 + 1
        a(lr2)    = v(i)
        indr(lr2) = jadd
        lenr(i)   = leni + 1
        lenU      = lenU + 1
    200 END DO

!     Normal exit.

    inform = 0
    go to 990

!     Not enough storage.

    970 inform = 7

    990 return

    end ! subroutine lu7add

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7cyc( kfirst, klast, ip )

    integer ::            ip(klast)

!     ------------------------------------------------------------------
!     lu7cyc performs a cyclic permutation on the row or column ordering
!     stored in ip, moving entry kfirst down to klast.
!     If kfirst .ge. klast, lu7cyc should not be called.
!     Sometimes klast = 0 and nothing should happen.

!     09 May 1988: First f77 version.
!     ------------------------------------------------------------------

    if (kfirst < klast) then
        ifirst = ip(kfirst)

        do 100 k = kfirst, klast - 1
            ip(k) = ip(k + 1)
        100 END DO

        ip(klast) = ifirst
    end if

    end ! subroutine lu7cyc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7elm( m, n, jelm, v, &
    lena, luparm, parmlu, &
    lenL, lenU, lrow, nrank, &
    a, indc, indr, ip, iq, lenr, locc, locr, &
    inform, diag )

    implicit           double precision (a-h,o-z)
    integer ::            luparm(30)
    double precision ::   parmlu(30), a(lena), v(m)
    integer ::            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
    integer ::            locc(n), locr(m)

!     ------------------------------------------------------------------
!     lu7elm  eliminates the subdiagonal elements of a vector  v(*),
!     where  L*v = y  for some vector y.
!     If  jelm > 0,  y  has just become column  jelm  of the matrix  A.
!     lu7elm  should not be called unless  m  is greater than  nrank.

!     inform = 0 if y contained no subdiagonal nonzeros to eliminate.
!     inform = 1 if y contained at least one nontrivial subdiagonal.
!     inform = 7 if there is insufficient storage.

!     09 May 1988: First f77 version.
!                  No longer calls lu7for at end.  lu8rpc, lu8mod do so.
!     ------------------------------------------------------------------

    parameter        ( zero = 0.0d+0 )

    small  = parmlu(3)
    nrank1 = nrank + 1
    diag   = zero

!     Compress row file if necessary.

    minfre = m - nrank
    nfree  = lena - lenL - lrow
    if (nfree >= minfre) go to 100
    call lu1rec( m, .TRUE. , luparm, lrow, lena, a, indr, lenr, locr )
    nfree  = lena - lenL - lrow
    if (nfree < minfre) go to 970

!     Pack the subdiagonals of  v  into  L,  and find the largest.

    100 vmax   = zero
    kmax   = 0
    l      = lena - lenL + 1

    do 200 k = nrank1, m
        i       = ip(k)
        vi      = abs( v(i) )
        if (vi <= small) go to 200
        l       = l - 1
        a(l)    = v(i)
        indc(l) = i
        if (vmax >= vi ) go to 200
        vmax    = vi
        kmax    = k
        lmax    = l
    200 END DO

    if (kmax == 0) go to 900

!     ------------------------------------------------------------------
!     Remove  vmax  by overwriting it with the last packed  v(i).
!     Then set the multipliers in  L  for the other elements.
!     ------------------------------------------------------------------

    imax       = ip(kmax)
    vmax       = a(lmax)
    a(lmax)    = a(l)
    indc(lmax) = indc(l)
    l1         = l + 1
    l2         = lena - lenL
    lenL       = lenL + (l2 - l)

    do 300 l = l1, l2
        a(l)    = - a(l) / vmax
        indr(l) =   imax
    300 END DO

!     Move the row containing vmax to pivotal position nrank + 1.

    ip(kmax  ) = ip(nrank1)
    ip(nrank1) = imax
    diag       = vmax

!     ------------------------------------------------------------------
!     If jelm is positive, insert  vmax  into a new row of  U.
!     This is now the only subdiagonal element.
!     ------------------------------------------------------------------

    if (jelm > 0) then
        lrow       = lrow + 1
        locr(imax) = lrow
        lenr(imax) = 1
        a(lrow)    = vmax
        indr(lrow) = jelm
    end if

    inform = 1
    go to 990

!     No elements to eliminate.

    900 inform = 0
    go to 990

!     Not enough storage.

    970 inform = 7

    990 return

    end ! subroutine lu7elm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7for( m, n, kfirst, klast, &
    lena, luparm, parmlu, &
    lenL, lenU, lrow, &
    a, indc, indr, ip, iq, lenr, locc, locr, &
    inform, diag )

    implicit           double precision (a-h,o-z)
    integer ::            luparm(30)
    double precision ::   parmlu(30), a(lena)
    integer ::            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
    integer ::            locc(n), locr(m)

!     ------------------------------------------------------------------
!     lu7for  (forward sweep) updates the LU factorization  A = L*U
!     when row  iw = ip(klast)  of  U  is eliminated by a forward
!     sweep of stabilized row operations, leaving  ip * U * iq  upper
!     triangular.

!     The row permutation  ip  is updated to preserve stability and/or
!     sparsity.  The column permutation  iq  is not altered.

!     kfirst  is such that row  ip(kfirst)  is the first row involved
!     in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero
!     in row  iw  in pivotal order.)  If  kfirst  is unknown it may be
!     input as  1.

!     klast   is such that row  ip(klast)  is the row being eliminated.
!     klast   is not altered.

!     lu7for  should be called only if  kfirst .le. klast.
!     If  kfirst = klast,  there are no nonzeros to eliminate, but the
!     diagonal element of row  ip(klast)  may need to be moved to the
!     front of the row.

!     On entry,  locc(*)  must be zero.

!     On exit:
!     inform = 0  if row iw has a nonzero diagonal (could be small).
!     inform = 1  if row iw has no diagonal.
!     inform = 7  if there is not enough storage to finish the update.

!     On a successful exit (inform le 1),  locc(*)  will again be zero.

!        Jan 1985: Final f66 version.
!     09 May 1988: First f77 version.
!     ------------------------------------------------------------------

    parameter        ( zero = 0.0d+0 )

    double precision ::   Ltol
    logical ::            swappd

    Ltol   = parmlu(2)
    small  = parmlu(3)
    uspace = parmlu(6)
    kbegin = kfirst
    swappd = .FALSE. 

!     We come back here from below if a row interchange is performed.

    100 iw     = ip(klast)
    lenw   = lenr(iw)
    if (lenw   ==   0  ) go to 910
    lw1    = locr(iw)
    lw2    = lw1 + lenw - 1
    jfirst = iq(kbegin)
    if (kbegin >= klast) go to 700

!     Make sure there is room at the end of the row file
!     in case row  iw  is moved there and fills in completely.

    minfre = n + 1
    nfree  = lena - lenL - lrow
    if (nfree < minfre) then
        call lu1rec( m, .TRUE. , luparm, lrow, lena, &
        a, indr, lenr, locr )
        lw1    = locr(iw)
        lw2    = lw1 + lenw - 1
        nfree  = lena - lenL - lrow
        if (nfree < minfre) go to 970
    end if

!     Set markers on row  iw.

    do 120 l = lw1, lw2
        j       = indr(l)
        locc(j) = l
    120 END DO


!     ==================================================================
!     Main elimination loop.
!     ==================================================================
    kstart = kbegin
    kstop  = min( klast, n )

    do 500 k  = kstart, kstop
        jfirst = iq(k)
        lfirst = locc(jfirst)
        if (lfirst == 0) go to 490

    !        Row  iw  has its first element in column  jfirst.

        wj     = a(lfirst)
        if (k == klast) go to 490

    !        ---------------------------------------------------------------
    !        We are about to use the first element of row  iv
    !               to eliminate the first element of row  iw.
    !        However, we may wish to interchange the rows instead,
    !        to preserve stability and/or sparsity.
    !        ---------------------------------------------------------------
        iv     = ip(k)
        lenv   = lenr(iv)
        lv1    = locr(iv)
        vj     = zero
        if (lenv      ==   0   ) go to 150
        if (indr(lv1) /= jfirst) go to 150
        vj     = a(lv1)
        if (            swappd               ) go to 200
        if (Ltol * abs( wj )  <  abs( vj )) go to 200
        if (Ltol * abs( vj )  <  abs( wj )) go to 150
        if (            lenv  <=  lenw     ) go to 200

    !        ---------------------------------------------------------------
    !        Interchange rows  iv  and  iw.
    !        ---------------------------------------------------------------
        150 ip(klast) = iv
        ip(k)     = iw
        kbegin    = k
        swappd    = .TRUE. 
        go to 600

    !        ---------------------------------------------------------------
    !        Delete the eliminated element from row  iw
    !        by overwriting it with the last element.
    !        ---------------------------------------------------------------
        200 a(lfirst)    = a(lw2)
        jlast        = indr(lw2)
        indr(lfirst) = jlast
        indr(lw2)    = 0
        locc(jlast)  = lfirst
        locc(jfirst) = 0
        lenw         = lenw - 1
        lenU         = lenU - 1
        if (lrow == lw2) lrow = lrow - 1
        lw2          = lw2  - 1

    !        ---------------------------------------------------------------
    !        Form the multiplier and store it in the  L  file.
    !        ---------------------------------------------------------------
        if (abs( wj ) <= small) go to 490
        amult   = - wj / vj
        l       = lena - lenL
        a(l)    = amult
        indr(l) = iv
        indc(l) = iw
        lenL    = lenL + 1

    !        ---------------------------------------------------------------
    !        Add the appropriate multiple of row  iv  to row  iw.
    !        We use two different inner loops.  The first one is for the
    !        case where row  iw  is not at the end of storage.
    !        ---------------------------------------------------------------
        if (lenv == 1) go to 490
        lv2    = lv1 + 1
        lv3    = lv1 + lenv - 1
        if (lw2 == lrow) go to 400

    !        ...............................................................
    !        This inner loop will be interrupted only if
    !        fill-in occurs enough to bump into the next row.
    !        ...............................................................
        do 350 lv = lv2, lv3
            jv     = indr(lv)
            lw     = locc(jv)
            if (lw > 0) then

            !              No fill-in.

                a(lw)  = a(lw)  +  amult * a(lv)
                if (abs( a(lw) ) <= small) then

                !                 Delete small computed element.

                    a(lw)     = a(lw2)
                    j         = indr(lw2)
                    indr(lw)  = j
                    indr(lw2) = 0
                    locc(j)   = lw
                    locc(jv)  = 0
                    lenU      = lenU - 1
                    lenw      = lenw - 1
                    lw2       = lw2  - 1
                end if
            else

            !              Row  iw  doesn't have an element in column  jv  yet
            !              so there is a fill-in.

                if (indr(lw2+1) /= 0) go to 360
                lenU      = lenU + 1
                lenw      = lenw + 1
                lw2       = lw2  + 1
                a(lw2)    = amult * a(lv)
                indr(lw2) = jv
                locc(jv)  = lw2
            end if
        350 END DO

        go to 490

    !        Fill-in interrupted the previous loop.
    !        Move row  iw  to the end of the row file.

        360 lv2      = lv
        locr(iw) = lrow + 1

        do 370 l = lw1, lw2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
            locc(j)    = lrow
        370 END DO

        lw1    = locr(iw)
        lw2    = lrow

    !        ...............................................................
    !        Inner loop with row  iw  at the end of storage.
    !        ...............................................................
        400 do 450 lv = lv2, lv3
            jv     = indr(lv)
            lw     = locc(jv)
            if (lw > 0) then

            !              No fill-in.

                a(lw)  = a(lw)  +  amult * a(lv)
                if (abs( a(lw) ) <= small) then

                !                 Delete small computed element.

                    a(lw)     = a(lw2)
                    j         = indr(lw2)
                    indr(lw)  = j
                    indr(lw2) = 0
                    locc(j)   = lw
                    locc(jv)  = 0
                    lenU      = lenU - 1
                    lenw      = lenw - 1
                    lw2       = lw2  - 1
                end if
            else

            !              Row  iw  doesn't have an element in column  jv  yet
            !              so there is a fill-in.

                lenU      = lenU + 1
                lenw      = lenw + 1
                lw2       = lw2  + 1
                a(lw2)    = amult * a(lv)
                indr(lw2) = jv
                locc(jv)  = lw2
            end if
        450 END DO

        lrow   = lw2

    !        The  k-th  element of row  iw  has been processed.
    !        Reset  swappd  before looking at the next element.

        490 swappd = .FALSE. 
    500 END DO

!     ==================================================================
!     End of main elimination loop.
!     ==================================================================

!     Cancel markers on row  iw.

    600 lenr(iw) = lenw
    if (lenw == 0) go to 910
    do 620 l = lw1, lw2
        j       = indr(l)
        locc(j) = 0
    620 END DO

!     Move the diagonal element to the front of row  iw.
!     At this stage,  lenw gt 0  and  klast le n.

    700 do 720 l = lw1, lw2
        ldiag = l
        if (indr(l) == jfirst) go to 730
    720 END DO
    go to 910

    730 diag        = a(ldiag)
    a(ldiag)    = a(lw1)
    a(lw1)      = diag
    indr(ldiag) = indr(lw1)
    indr(lw1)   = jfirst

!     If an interchange is needed, repeat from the beginning with the
!     new row  iw,  knowing that the opposite interchange cannot occur.

    if (swappd) go to 100
    inform = 0
    go to 950

!     Singular.

    910 diag   = zero
    inform = 1

!     Force a compression if the file for  U  is much longer than the
!     no. of nonzeros in  U  (i.e. if  lrow  is much bigger than  lenU).
!     This should prevent memory fragmentation when there is far more
!     memory than necessary  (i.e. when  lena  is huge).

    950 limit  = uspace * lenU + m + n + 1000
    if (lrow > limit) then
        call lu1rec( m, .TRUE. , luparm, lrow, lena, &
        a, indr, lenr, locr )
    end if
    go to 990

!     Not enough storage.

    970 inform = 7

!     Exit.

    990 return

    end ! subroutine lu7for

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7rnk( m, n, jsing, &
    lena, luparm, parmlu, &
    lenL, lenU, lrow, nrank, &
    a, indc, indr, ip, iq, lenr, locc, locr, &
    inform, diag )

    implicit           double precision (a-h,o-z)
    integer ::            luparm(30)
    double precision ::   parmlu(30), a(lena)
    integer ::            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
    integer ::            locc(n), locr(m)

!     ------------------------------------------------------------------
!     lu7rnk (check rank) assumes U is currently nrank by n
!     and determines if row nrank contains an acceptable pivot.
!     If not, the row is deleted and nrank is decreased by 1.

!     jsing is an input parameter (not altered).  If jsing is positive,
!     column jsing has already been judged dependent.  A substitute
!     (if any) must be some other column.

!     -- Jul 1987: First version.
!     09 May 1988: First f77 version.
!     ------------------------------------------------------------------

    parameter        ( zero = 0.0d+0 )

    Utol1    = parmlu(4)
    diag     = zero

!     Find Umax, the largest element in row nrank.

    iw       = ip(nrank)
    lenw     = lenr(iw)
    if (lenw == 0) go to 400
    l1       = locr(iw)
    l2       = l1 + lenw - 1
    Umax     = zero
    lmax     = l1

    do 100 l = l1, l2
        if (Umax < abs( a(l) )) then
            Umax   =  abs( a(l) )
            lmax   =  l
        end if
    100 END DO

!     Find which column that guy is in (in pivotal order).
!     Interchange him with column nrank, then move him to be
!     the new diagonal at the front of row nrank.

    diag   = a(lmax)
    jmax   = indr(lmax)

    do 300 kmax = nrank, n
        if (iq(kmax) == jmax) go to 320
    300 END DO

    320 iq(kmax)  = iq(nrank)
    iq(nrank) = jmax
    a(lmax)   = a(l1)
    a(l1)     = diag
    indr(lmax)= indr(l1)
    indr(l1)  = jmax

!     See if the new diagonal is big enough.

    if (Umax <= Utol1) go to 400
    if (jmax == jsing) go to 400

!     ------------------------------------------------------------------
!     The rank stays the same.
!     ------------------------------------------------------------------
    inform = 0
    return

!     ------------------------------------------------------------------
!     The rank decreases by one.
!     ------------------------------------------------------------------
    400 inform = -1
    nrank  = nrank - 1
    if (lenw > 0) then

    !        Delete row nrank from U.

        lenU     = lenU - lenw
        lenr(iw) = 0
        do 420 l = l1, l2
            indr(l) = 0
        420 END DO

        if (l2 == lrow) then

        !           This row was at the end of the data structure.
        !           We have to reset lrow.
        !           Preceding rows might already have been deleted, so we
        !           have to be prepared to go all the way back to 1.

            do 450 l = 1, l2
                if (indr(lrow) > 0) go to 900
                lrow  = lrow - 1
            450 END DO
        end if
    end if

    900 return

    end ! subroutine lu7rnk

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu7zap( m, n, jzap, kzap, &
    lena, lenU, lrow, nrank, &
    a, indr, ip, iq, lenr, locr )

    implicit           double precision (a-h,o-z)
    double precision ::   a(lena)
    integer ::            indr(lena), ip(m), iq(n), lenr(m)
    integer ::            locr(m)

!     ------------------------------------------------------------------
!     lu7zap  eliminates all nonzeros in column  jzap  of  U.
!     It also sets  kzap  to the position of  jzap  in pivotal order.
!     Thus, on exit we have  iq(kzap) = jzap.

!     -- Jul 1987: nrank added.
!     10 May 1988: First f77 version.
!     ------------------------------------------------------------------

    do 100 k  = 1, nrank
        i      = ip(k)
        leni   = lenr(i)
        if (leni == 0) go to 90
        lr1    = locr(i)
        lr2    = lr1 + leni - 1
        do 50 l = lr1, lr2
            if (indr(l) == jzap) go to 60
        50 END DO
        go to 90

    !        Delete the old element.

        60 a(l)      = a(lr2)
        indr(l)   = indr(lr2)
        indr(lr2) = 0
        lenr(i)   = leni - 1
        lenU      = lenU - 1

    !        Stop if we know there are no more rows containing  jzap.

        90 kzap   = k
        if (iq(k) == jzap) go to 800
    100 END DO

!     nrank must be smaller than n because we haven't found kzap yet.

    do 200 k = nrank+1, n
        kzap  = k
        if (iq(k) == jzap) go to 800
    200 END DO

!     See if we zapped the last element in the file.

    800 if (lrow > 0) then
        if (indr(lrow) == 0) lrow = lrow - 1
    end if

    end ! subroutine lu7zap
!***********************************************************************

!     File  lusol8a.f

!     lu8rpc

!     Sparse LU update: Replace Column
!     LUSOL's sparse implementation of the Bartels-Golub update.

! 01 May 2002: Derived from LUSOL's original lu8a.f file.
! 01 May 2002: Current version of lusol8a.f.
! 15 Sep 2004: Test nout. gt. 0 to protect write statements.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine lu8rpc( mode1, mode2, m, n, jrep, v, w, &
    lena, luparm, parmlu, &
    a, indc, indr, ip, iq, &
    lenc, lenr, locc, locr, &
    inform, diag, vnorm )

    implicit           double precision(a-h,o-z)
    integer ::            luparm(30)
    double precision ::   parmlu(30), a(lena), v(m), w(n)
    integer ::            indc(lena), indr(lena), ip(m), iq(n)
    integer ::            lenc(n), lenr(m)
    integer ::            locc(n), locr(m)

!     ------------------------------------------------------------------
!     lu8rpc  updates the LU factorization  A = L*U  when column  jrep
!     is replaced by some vector  a(new).

!     lu8rpc  is an implementation of the Bartels-Golub update,
!     designed for the case where A is rectangular and/or singular.
!     L is a product of stabilized eliminations (m x m, nonsingular).
!     P U Q is upper trapezoidal (m x n, rank nrank).

!     If  mode1 = 0,  the old column is taken to be zero
!                     (so it does not have to be removed from  U).

!     If  mode1 = 1,  the old column need not have been zero.

!     If  mode2 = 0,  the new column is taken to be zero.
!                     v(*)  is not used or altered.

!     If  mode2 = 1,  v(*)  must contain the new column  a(new).
!                     On exit,  v(*)  will satisfy  L*v = a(new).

!     If  mode2 = 2,  v(*)  must satisfy  L*v = a(new).

!     The array  w(*)  is not used or altered.

!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform ne 7), this will again be true.

!     On exit:
!     inform = -1  if the rank of U decreased by 1.
!     inform =  0  if the rank of U stayed the same.
!     inform =  1  if the rank of U increased by 1.
!     inform =  2  if the update seemed to be unstable
!                  (diag much bigger than vnorm).
!     inform =  7  if the update was not completed (lack of storage).
!     inform =  8  if jrep is not between 1 and n.

!     -- Jan 1985: Original F66 version.
!     -- Jul 1987: Modified to maintain U in trapezoidal form.
!     10 May 1988: First f77 version.
!     16 Oct 2000: Added test for instability (inform = 2).
!     ------------------------------------------------------------------

    logical ::            singlr
    parameter        ( zero = 0.0d+0 )

    nout   = luparm(1)
    lprint = luparm(2)
    nrank  = luparm(16)
    lenL   = luparm(23)
    lenU   = luparm(24)
    lrow   = luparm(25)
    Utol1  = parmlu(4)
    Utol2  = parmlu(5)
    nrank0 = nrank
    diag   = zero
    vnorm  = zero
    if (jrep < 1) go to 980
    if (jrep > n) go to 980

!     ------------------------------------------------------------------
!     If mode1 = 0, there are no elements to be removed from  U
!     but we still have to set  krep  (using a backward loop).
!     Otherwise, use lu7zap to remove column  jrep  from  U
!     and set  krep  at the same time.
!     ------------------------------------------------------------------
    if (mode1 == 0) then
        krep   = n + 1

        10 krep   = krep - 1
        if (iq(krep) /= jrep) go to 10
    else
        call lu7zap( m, n, jrep, krep, &
        lena, lenU, lrow, nrank, &
        a, indr, ip, iq, lenr, locr )
    end if

!     ------------------------------------------------------------------
!     Insert a new column of u and find klast.
!     ------------------------------------------------------------------

    if (mode2 == 0) then
        klast  = 0
    else
        if (mode2 == 1) then

        !           Transform v = a(new) to satisfy  L*v = a(new).

            call lu6sol( 1, m, n, v, w, lena, luparm, parmlu, &
            a, indc, indr, ip, iq, &
            lenc, lenr, locc, locr, inform )
        end if

    !        Insert into  U  any nonzeros in the top of  v.
    !        row  ip(klast)  will contain the last nonzero in pivotal order.
    !        Note that  klast  will be in the range  ( 0, nrank ).

        call lu7add( m, n, jrep, v, &
        lena, luparm, parmlu, &
        lenL, lenU, lrow, nrank, &
        a, indr, ip, lenr, locr, &
        inform, klast, vnorm )
        if (inform == 7) go to 970
    end if

!     ------------------------------------------------------------------
!     In general, the new column causes U to look like this:

!                 krep        n                 krep  n

!                ....a.........          ..........a...
!                 .  a        .           .        a  .
!                  . a        .            .       a  .
!                   .a        .             .      a  .
!        P U Q =     a        .    or        .     a  .
!                    b.       .               .    a  .
!                    b .      .                .   a  .
!                    b  .     .                 .  a  .
!                    b   ......                  ..a...  nrank
!                    c                             c
!                    c                             c
!                    c                             c     m

!     klast points to the last nonzero "a" or "b".
!     klast = 0 means all "a" and "b" entries are zero.
!     ------------------------------------------------------------------

    if (mode2 == 0) then
        if (krep > nrank) go to 900
    else if (nrank < m) then

    !        Eliminate any "c"s (in either case).
    !        Row nrank + 1 may end up containing one nonzero.

        call lu7elm( m, n, jrep, v, &
        lena, luparm, parmlu, &
        lenL, lenU, lrow, nrank, &
        a, indc, indr, ip, iq, lenr, locc, locr, &
        inform, diag )
        if (inform == 7) go to 970

        if (inform == 1) then

        !           The nonzero is apparently significant.
        !           Increase nrank by 1 and make klast point to the bottom.

            nrank = nrank + 1
            klast = nrank
        end if
    end if

    if (nrank < n) then

    !        The column rank is low.
    
    !        In the first case, we want the new column to end up in
    !        position nrank, so the trapezoidal columns will have a chance
    !        later on (in lu7rnk) to pivot in that position.
    
    !        Otherwise the new column is not part of the triangle.  We
    !        swap it into position nrank so we can judge it for singularity.
    !        lu7rnk might choose some other trapezoidal column later.

        if (krep < nrank) then
            klast     = nrank
        else
            iq(krep ) = iq(nrank)
            iq(nrank) = jrep
            krep      = nrank
        end if
    end if

!     ------------------------------------------------------------------
!     If krep .lt. klast, there are some "b"s to eliminate:

!                  krep

!                ....a.........
!                 .  a        .
!                  . a        .
!                   .a        .
!        P U Q =     a        .  krep
!                    b.       .
!                    b .      .
!                    b  .     .
!                    b   ......  nrank

!     If krep .eq. klast, there are no "b"s, but the last "a" still
!     has to be moved to the front of row krep (by lu7for).
!     ------------------------------------------------------------------

    if (krep <= klast) then

    !        Perform a cyclic permutation on the current pivotal order,
    !        and eliminate the resulting row spike.  krep becomes klast.
    !        The final diagonal (if any) will be correctly positioned at
    !        the front of the new krep-th row.  nrank stays the same.

        call lu7cyc( krep, klast, ip )
        call lu7cyc( krep, klast, iq )

        call lu7for( m, n, krep, klast, &
        lena, luparm, parmlu, &
        lenL, lenU, lrow, &
        a, indc, indr, ip, iq, lenr, locc, locr, &
        inform, diag )
        if (inform == 7) go to 970
        krep   = klast

    !        Test for instability (diag much bigger than vnorm).

        singlr = vnorm < Utol2 * abs( diag )
        if ( singlr ) go to 920
    end if

!     ------------------------------------------------------------------
!     Test for singularity in column krep (where krep .le. nrank).
!     ------------------------------------------------------------------

    diag   = zero
    iw     = ip(krep)
    singlr = lenr(iw) == 0

    if ( .NOT. singlr) then
        l1     = locr(iw)
        j1     = indr(l1)
        singlr = j1 /= jrep

        if ( .NOT. singlr) then
            diag   = a(l1)
            singlr = abs( diag ) <= Utol1          .OR. &
            abs( diag ) <= Utol2 * vnorm
        end if
    end if

    if ( singlr  .AND.  krep < nrank ) then

    !        Perform cyclic permutations to move column jrep to the end.
    !        Move the corresponding row to position nrank
    !        then eliminate the resulting row spike.

        call lu7cyc( krep, nrank, ip )
        call lu7cyc( krep, n    , iq )

        call lu7for( m, n, krep, nrank, &
        lena, luparm, parmlu, &
        lenL, lenU, lrow, &
        a, indc, indr, ip, iq, lenr, locc, locr, &
        inform, diag )
        if (inform == 7) go to 970
    end if

!     Find the best column to be in position nrank.
!     If singlr, it can't be the new column, jrep.
!     If nothing satisfactory exists, nrank will be decreased.

    if ( singlr  .OR.  nrank < n ) then
        jsing  = 0
        if ( singlr ) jsing = jrep

        call lu7rnk( m, n, jsing, &
        lena, luparm, parmlu, &
        lenL, lenU, lrow, nrank, &
        a, indc, indr, ip, iq, lenr, locc, locr, &
        inform, diag )
    end if

!     ------------------------------------------------------------------
!     Set inform for exit.
!     ------------------------------------------------------------------

    900 if (nrank == nrank0) then
        inform =  0
    else if (nrank < nrank0) then
        inform = -1
        if (nrank0 == n) then
            if (nout > 0  .AND.  lprint >= 0) &
            write(nout, 1100) jrep, diag
        end if
    else
        inform =  1
    end if
    go to 990

!     Instability.

    920 inform = 2
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1200) jrep, diag
    go to 990

!     Not enough storage.

    970 inform = 7
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1700) lena
    go to 990

!     jrep  is out of range.

    980 inform = 8
    if (nout > 0  .AND.  lprint >= 0) &
    write(nout, 1800) m, n, jrep

!     Exit.

    990 luparm(10) = inform
    luparm(15) = luparm(15) + 1
    luparm(16) = nrank
    luparm(23) = lenL
    luparm(24) = lenU
    luparm(25) = lrow
    return

    1100 format(/ ' lu8rpc  warning.  Singularity after replacing column.', &
    '    jrep =', i8, '    diag =', 1p, e12.2 )
    1200 format(/ ' lu8rpc  warning.  Instability after replacing column.', &
    '    jrep =', i8, '    diag =', 1p, e12.2 )
    1700 format(/ ' lu8rpc  error...  Insufficient storage.', &
    '    lena =', i8)
    1800 format(/ ' lu8rpc  error...  jrep  is out of range.', &
    '    m =', i8, '    n =', i8, '    jrep =', i8)

    end ! subroutine lu8rpc
