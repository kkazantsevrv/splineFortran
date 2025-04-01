module CubicSplineModule
    implicit none

    integer, parameter :: wp = kind(1.0d0), nn = 10
    
    type CubicSpline
        real(wp), dimension(nn) :: x, y, M
        real(wp), dimension(nn-1) :: h
    contains
        procedure :: init
        procedure :: div
        procedure :: assemble
        procedure :: sinta
        ! procedure :: sintall
        ! procedure :: spl
    end type CubicSpline

contains
    subroutine linear_solve(a, b, c, d, x, c1, c2)
        real(wp), intent(in) :: a(:), b(:), c(:), d(:)
        real(wp), intent(out) :: x(:)
        real(wp), intent(in) :: c1, c2
        real(wp), allocatable :: aa(:), bb(:), cc(:), dd(:)
        integer :: n, i

        ! Проверка размеров массивов
        if (size(a) < 2 .or. size(b) < 2 .or. size(c) < 2 .or. size(d) < 2) then
            error stop "Arrays must have at least 2 elements"
        end if

        ! Выделяем подмассивы (игнорируем первый и последний элементы)
        n = size(d) - 2
        allocate(aa(n), bb(n), cc(n), dd(n))

        aa = a(2:n+1)
        bb = b(2:n+1)
        cc = c(2:n+1)
        dd = d(2:n+1)

        ! Прямой ход прогонки
        do i = 2, n
            dd(i) = dd(i) - dd(i-1) / bb(i-1) * aa(i-1)
            bb(i) = bb(i) - cc(i-1) / bb(i-1) * aa(i-1)
        end do

        ! Обратный ход прогонки
        do i = n-1, 1, -1
            dd(i) = dd(i) - dd(i+1) / bb(i+1) * cc(i)
        end do

        ! Решение для внутренних точек
        x(2:n+1) = dd(1:n) / bb(1:n)

        ! Граничные условия
        x(1) = (d(1) - c(1)*x(2) - c1*x(3)) / b(1)
        x(n+2) = (d(n+2) - a(n+1)*x(n+1) - c2*x(n)) / b(n+2)
        deallocate(aa, bb, cc, dd)
    end subroutine linear_solve

    subroutine init(this, x, y)
        class(CubicSpline), intent(inout) :: this
        real(wp), intent(in) :: x(nn), y(nn)
        this%x = x
        this%y = y

        this%h = this%x(2:) - this%x(:nn-1)
    end subroutine init

    function div(this, i, j, k) result(res)
        class(CubicSpline), intent(in) :: this
        integer, intent(in) :: i, j, k
        real(wp) :: res
        res = ((this%y(k) - this%y(j))/(this%x(k) - this%x(j)) - &
                    (this%y(j) - this%y(i))/(this%x(j) - this%x(i))) / &
                    (this%x(k) - this%x(i))
    end function div

    subroutine assemble(this)
        class(CubicSpline), intent(inout) :: this
        integer :: i, n
        real(wp), dimension(nn-1) :: mu1, lam1, h
        real(wp), dimension(nn) :: x, y, d, b, M
        real(wp), dimension(nn-1) :: a, c
        real(wp), dimension(nn-2) :: hsum, hdel2, hdel1, mu, lam
        real(wp) :: c1, c2

        n = nn
        x = this%x
        y = this%y
        h = this%h
        d = 0.0_wp
        mu = h(:n-1)/(h(:n-1) + h(2:))
        lam = 1.0_wp - mu

        mu1(:n-2) = mu
        lam1(2:) = lam
        mu1(n-1) = -1.0_wp
        lam1(1) = -1.0_wp
        mu1(n-2) = mu(n-2) - lam(n-2)
        mu1(1) = 0.0_wp
        lam1(n-1) = 0.0_wp
        lam1(2) = lam(1) - mu(1)

        hsum = h(:n-2) + h(2:)
        hdel2 = h(2:)
        hdel1 = h(:n-2)

        b = 2.0_wp
        c = lam1
        a = mu1
        c1 = mu(1)
        c2 = lam(n-2)
        b(1) = lam(1)
        b(2) = 1.0_wp + lam(1)
        b(n) = mu(n-2)
        b(nn-1) = 1 + mu(nn-2)

        do i = 1, n-2
            d(i+1) = 6.0_wp/hsum(i)*((y(i+2) - y(i+1))/hdel2(i) - (y(i+1) - y(i))/hdel1(i))
        end do

        d(1) = 0.0_wp
        d(n) = 0.0_wp
        d(2) = d(2)*lam(1)
        d(n-1) = d(n-1)*mu(n-2)
        call linear_solve(a, b, c, d, M, c1, c2)
        this%M = M
    end subroutine assemble 

    subroutine sinta(this, ip, ls1)
        class(CubicSpline), intent(in) :: this
        integer, intent(in) :: ip
        real(wp), intent(out) :: ls1(nn)
        integer :: i, n
        real(wp) :: ls(nn-1)
        n = nn

        ls = 0.0_wp
        ls1 = 0.0_wp
        do i = 1, n-1
            ls(i) = -this%h(i)**3/24.0_wp*(this%M(i) + this%M(i+1))
        end do
        if ( ip /= 1 .and. ip /= n ) then
            ls(ip-1) = ls(ip-1) + this%h(ip-1)/2.0_wp
            ls(ip) = ls(ip) + this%h(ip)/2.0_wp
        end if
        if (ip == 1) then
            ls(1) = ls(1) + this%h(1)/2.0_wp
        end if
        if (ip == n) then
            ls(n-1) = ls(n-1) + this%h(n-1)/2.0_wp
        end if
        ls1(2:) = ls
        do i = 1, n-1
            ls1(i+1) = ls1(i+1) + ls1(i)
        end do
        print*, ls1
    end subroutine sinta

    subroutine sintall(this, ip, ls1)
        class(CubicSpline), intent(in) :: this
        integer, intent(in) :: ip
        real(wp), intent(out) :: ls1(nn)
        integer :: i, n
        real(wp) :: ls(nn-1)
        n = nn

        ls = 0.0_wp
        ls1 = 0.0_wp
        do i = 1, n-1
            ls(i) = -this%h(i)**3/24.0_wp*(this%M(i) + this%M(i+1))
        end do
        if ( ip /= 1 .and. ip /= n ) then
            ls(ip-1) = ls(ip-1) + this%h(ip-1)/2.0_wp
            ls(ip) = ls(ip) + this%h(ip)/2.0_wp
        end if
        if (ip == 1) then
            ls(1) = ls(1) + this%h(1)/2.0_wp
        end if
        if (ip == n) then
            ls(n-1) = ls(n-1) + this%h(n-1)/2.0_wp
        end if
        ls1(2:) = ls
        do i = 1, n-1
            ls1(i+1) = ls1(i+1) + ls1(i)
        end do
        print*, ls1
    end subroutine sinta
end module CubicSplineModule

program test_spline
    use CubicSplineModule
    implicit none

    type(CubicSpline) :: spline
    real(wp) :: x(nn), y(nn)
    real(wp) :: point, interpolated_value, ls(nn)
    integer :: i
    do i = 1, nn
        x(i) = ((i-1) * 1.0_wp / (nn-1))**2
        y(i) = 0.0_wp
    end do
    do i = 1, nn
        y(i) = 1.0_wp
        call spline%init(x, y)
        call spline%assemble()
        call spline%sinta(i, ls)
        y(i) = 0.0_wp
    end do

end program test_spline