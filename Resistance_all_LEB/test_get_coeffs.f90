program test
    use ACCESS_Constants, only: rk=>dp, npts, ninteg
    use DryDep, only: GetDryDepExCoeffs

    implicit none

    integer :: i, j
    real(rk) :: rb(npts, ninteg)
    real(rk) :: rc(npts, ninteg)
    real(rk) :: rm(npts, ninteg)
    real(rk) :: rs(npts, ninteg)

    print *, "rs before call:"  ! (uninitialized)
    print "(*(g0))", ((rs(i,j), " ", j = 1, ninteg), new_line("A"), i = 1, npts)
    call GetDryDepExCoeffs(rb,rc,rm,rs)
    print *, "rs after call:"
    print "(*(g0))", ((rs(i,j), " ", j = 1, ninteg), new_line("A"), i = 1, npts)

end program test
