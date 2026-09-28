!=======================================================================
!  IPS  -  N corpos (coordenadas heliocentricas canonicas + J2/J4)
!          integrados com RADAU (RA15) + FFT das series temporais
!
!  Versao enxuta dos programas ips_thalassa.f / ips_n_naiad.f /
!  ips_j_am.f / ssatfftm.f. Tudo que mudava entre eles (planeta,
!  satelites, ressonancia, faixa de a e e) agora fica nos arquivos:
!
!     ips.in        -> configuracao (planeta, integracao, variacao, saida)
!     planetas.in   -> satelites: massa e elementos orbitais
!
!  Compilar:   gfortran -O2 -o ips ips.f90
!              gfortran -O2 -fopenmp -o ips ips.f90   (paralelo, opcional)
!  Rodar:      ./ips            (le ips.in)
!              ./ips outro.in
!
!  Unidades internas: distancia = raio equatorial do planeta,
!  massa = massa do planeta, tempo = dia.
!=======================================================================
module sistema
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp), pi2 = 2.0_dp*pi
  real(dp), parameter :: conv = pi/180.0_dp
  integer :: N, neq                          ! n. de corpos, n. de equacoes (6N)
  real(dp) :: G, dj2, dj4
  real(dp), parameter :: dm0 = 1.0_dp        ! massa do planeta central
  real(dp), allocatable :: dm(:), ami(:), dmimi(:)
contains

  ! reduz angulo para [0, 2pi)
  pure real(dp) function mod2pi(x)
    real(dp), intent(in) :: x
    mod2pi = modulo(x, pi2)
  end function

end module sistema

!=======================================================================
module dinamica
  use sistema
  implicit none
contains

  ! Equacoes de movimento canonicas. x = (x,y,z,px,py,pz) de cada corpo.
  ! f1(1:3) = dx/dt (velocidade heliocentrica), f1(4:6) = dp/dt.
  subroutine force(x, f1)
    real(dp), intent(in)  :: x(neq)
    real(dp), intent(out) :: f1(neq)
    real(dp) :: ptot(3), ri(3), d(3), r3, d3, gj(3)
    integer :: i, j, l, m

    ptot = 0.0_dp
    do j = 1, N
      ptot = ptot + x(6*j-2:6*j)
    end do

    do i = 1, N
      l  = 6*(i-1)
      ri = x(l+1:l+3)
      ! dx/dt = p_i/m_i + (soma dos p)/m0
      f1(l+1:l+3) = dmimi(i)*x(l+4:l+6) + (ptot - x(l+4:l+6))/dm0

      r3 = (ri(1)**2 + ri(2)**2 + ri(3)**2)**1.5_dp
      f1(l+4:l+6) = -G*dm(i)*dm0*ri/r3

      do j = 1, N
        if (j == i) cycle
        m  = 6*(j-1)
        d  = ri - x(m+1:m+3)
        d3 = (d(1)**2 + d(2)**2 + d(3)**2)**1.5_dp
        f1(l+4:l+6) = f1(l+4:l+6) - G*dm(i)*dm(j)*d/d3
      end do

      call achat(ri, dm(i), gj)
      f1(l+4:l+6) = f1(l+4:l+6) - gj
    end do
  end subroutine force

  ! Gradiente dos termos de achatamento J2 e J4 (raio do planeta = 1).
  subroutine achat(c, dmm, gj)
    real(dp), intent(in)  :: c(3), dmm
    real(dp), intent(out) :: gj(3)
    real(dp) :: r, z, r5, r7, r9, r11, dp1(3), dp2(3), dp3(3), dp4(3), cj2, cj4

    z   = c(3)
    r   = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
    r5  = r**5; r7 = r**7; r9 = r**9; r11 = r**11
    dp1 = -5.0_dp*c/r7         ! d(r^-5)/dc
    dp2 = -3.0_dp*c/r5         ! d(r^-3)/dc
    dp3 = -9.0_dp*c/r11        ! d(r^-9)/dc
    dp4 = -7.0_dp*c/r9         ! d(r^-7)/dc

    cj2 = G*dm0*dmm*dj2/2.0_dp
    cj4 = G*dm0*dmm*dj4/8.0_dp

    gj = cj2*(3.0_dp*z*z*dp1 - dp2) &
       + cj4*(35.0_dp*z**4*dp3 - 30.0_dp*z*z*dp4 + 3.0_dp*dp1)
    gj(3) = gj(3) + cj2*6.0_dp*z/r5 &
          + cj4*(140.0_dp*z**3/r9 - 60.0_dp*z/r7)
  end subroutine achat

  ! Elementos -> estado canonico (posicoes heliocentricas e momenta).
  ! el(i,:) = (a, e, I, w, Om, M) em unidades internas / radianos.
  subroutine elem_para_estado(el, x)
    real(dp), intent(in)  :: el(N,6)
    real(dp), intent(out) :: x(neq)
    real(dp) :: r(3), v(3,N), mv(3), msum
    integer :: i, l

    mv = 0.0_dp; msum = 0.0_dp
    do i = 1, N
      l = 6*(i-1)
      call orbxyz(el(i,:), ami(i), r, v(:,i))
      x(l+1:l+3) = r
      mv   = mv + dm(i)*v(:,i)
      msum = msum + dm(i)
    end do
    ! Momenta: resolve v_i = p_i/m_i + P/m0 (antes feito com gaussj)
    mv = mv/(1.0_dp + msum/dm0)            ! = P (momento total)
    do i = 1, N
      l = 6*(i-1)
      x(l+4:l+6) = dm(i)*(v(:,i) - mv/dm0)
    end do
  end subroutine elem_para_estado

  ! Estado canonico -> elementos (a, e, I, w, Om, M) de cada corpo.
  subroutine estado_para_elem(x, el)
    real(dp), intent(in)  :: x(neq)
    real(dp), intent(out) :: el(N,6)
    real(dp) :: ptot(3), v(3)
    integer :: i, j, l

    ptot = 0.0_dp
    do j = 1, N
      ptot = ptot + x(6*j-2:6*j)
    end do
    do i = 1, N
      l = 6*(i-1)
      v = dmimi(i)*x(l+4:l+6) + (ptot - x(l+4:l+6))/dm0
      call xyzorb(x(l+1:l+3), v, ami(i), el(i,:))
    end do
  end subroutine estado_para_elem

  ! Elementos -> posicao e velocidade (relativas ao planeta)
  subroutine orbxyz(el, mu, r, v)
    real(dp), intent(in)  :: el(6), mu
    real(dp), intent(out) :: r(3), v(3)
    real(dp) :: a, e, di, w, om, u, f, rr, rp, rfp, cfw, sfw, coom, siom, coi, sini

    a = el(1); e = el(2); di = el(3); w = el(4); om = el(5)
    call kepler(el(6), e, u, f)

    rr  = a*(1.0_dp - e*cos(u))
    rp  = sqrt(mu*a)*e*sin(u)/rr
    rfp = sqrt(mu*a*(1.0_dp - e*e))/rr
    cfw = cos(w+f);  sfw = sin(w+f)
    coom = cos(om);  siom = sin(om)
    coi = cos(di);   sini = sin(di)

    v(1) = rp*cfw*coom - rfp*sfw*coom - rp*sfw*siom*coi - rfp*cfw*siom*coi
    v(2) = rp*cfw*siom - rfp*sfw*siom + rp*sfw*coom*coi + rfp*cfw*coom*coi
    v(3) = rp*sfw*sini + rfp*cfw*sini

    r(1) = rr*(cfw*coom - sfw*coi*siom)
    r(2) = rr*(cfw*siom + sfw*coi*coom)
    r(3) = rr*sini*sfw
  end subroutine orbxyz

  ! Equacao de Kepler (Newton-Raphson): anomalia media -> excentrica e verdadeira
  subroutine kepler(am, e, u, f)
    real(dp), intent(in)  :: am, e
    real(dp), intent(out) :: u, f
    real(dp) :: u0
    integer :: it

    u0 = am + e*sin(am)
    do it = 1, 35
      u = u0 - (u0 - e*sin(u0) - am)/(1.0_dp - e*cos(u0))
      if (abs(u-u0) <= 1.0e-12_dp .and. abs(u - e*sin(u) - am) <= 1.0e-12_dp) exit
      u0 = u
    end do
    f = mod2pi(atan2(sqrt(1.0_dp - e*e)*sin(u), cos(u) - e))
  end subroutine kepler

  ! Posicao e velocidade -> elementos (Brouwer). Angulos em [0, 2pi).
  subroutine xyzorb(r, v, mu, el)
    real(dp), intent(in)  :: r(3), v(3), mu
    real(dp), intent(out) :: el(6)
    real(dp) :: rr, v2, a, esinu, ecosu, e, u, f, h, cosi, di, om, w, alat
    real(dp) :: coslat, sinlat

    rr = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    v2 = v(1)**2 + v(2)**2 + v(3)**2
    a  = mu/(2.0_dp*mu/rr - v2)
    esinu = (r(1)*v(1) + r(2)*v(2) + r(3)*v(3))/sqrt(mu*a)
    ecosu = rr*v2/mu - 1.0_dp
    e  = sqrt(esinu**2 + ecosu**2)
    u  = mod2pi(atan2(esinu, ecosu))
    f  = mod2pi(atan2(sqrt(1.0_dp - e*e)*sin(u), cos(u) - e))

    h    = sqrt(mu*a*(1.0_dp - e*e))
    cosi = (r(1)*v(2) - r(2)*v(1))/h
    if (cosi >= 1.0_dp) then
      ! orbita plana: Om = 0 e w vira longitude do pericentro
      di = 0.0_dp
      om = 0.0_dp
      alat = atan2(r(2), r(1))
    else
      di = acos(cosi)
      om = mod2pi(atan2(r(2)*v(3) - r(3)*v(2), -(r(3)*v(1) - r(1)*v(3))))
      coslat = r(1)*cos(om) + r(2)*sin(om)
      sinlat = (r(2)*cos(om) - r(1)*sin(om))*cos(di) + r(3)*sin(di)
      alat = atan2(sinlat, coslat)
    end if
    if (e <= 1.0e-11_dp) then
      ! circular: w indefinido, uso w = 0 e M = argumento da latitude
      w = 0.0_dp
      u = mod2pi(alat)
    else
      w = mod2pi(alat - f)
    end if

    el = [a, e, di, w, om, mod2pi(u - e*sin(u))]
  end subroutine xyzorb

end module dinamica

!=======================================================================
!  Integrador RADAU de Everhart (ordem 15), mesmo algoritmo do original.
!  Integra x' = F(x) de t = 0 ate tf (NCLASS = 1).
!=======================================================================
module radau
  use sistema
  use dinamica, only: force
  implicit none
contains

  subroutine ra15(x, tf, ll, nv)
    integer,  intent(in)    :: ll, nv
    real(dp), intent(inout) :: x(nv)
    real(dp), intent(in)    :: tf
    real(dp) :: f1(nv), fj(nv), y(nv), b(7,nv), g(7,nv), e(7,nv), bd(7,nv)
    real(dp) :: c(21), d(21), r(21), w(7)
    integer,  parameter :: nw(8) = [0, 0, 1, 3, 6, 10, 15, 21]
    real(dp), parameter :: h(8) = [0.0_dp, &
      .05626256053692215_dp, .18024069173689236_dp, .35262471711316964_dp, &
      .54715362633055538_dp, .73421017721541053_dp, .88532094683909577_dp, &
      .97752061356128750_dp]
    real(dp), parameter :: sr = 1.4_dp
    logical  :: nsf, nper
    real(dp) :: dir, pw, ss, tp, t, tval, s, a, hv, temp, gk, tm, q
    integer  :: n, k, l, la, lb, lc, ld, le, ncount, ns, ni, m, j, jd

    nper = .false.
    nsf  = .false.
    dir  = sign(1.0_dp, tf)
    pw   = 1./9.                   ! (precisao simples, como no original)

    do n = 2, 8
      w(n-1) = 1.0_dp/n
    end do
    b  = 0.0_dp
    bd = 0.0_dp

    c(1) = -h(2)
    d(1) =  h(2)
    r(1) = 1.0_dp/(h(3) - h(2))
    la = 1
    lc = 1
    do k = 3, 7
      lb = la
      la = lc + 1
      lc = nw(k+1)
      c(la) = -h(k)*c(lb)
      c(lc) = c(la-1) - h(k)
      d(la) = h(2)*d(lb)
      d(lc) = -c(lc)
      r(la) = 1.0_dp/(h(k+1) - h(2))
      r(lc) = 1.0_dp/(h(k+1) - h(k))
      do l = 4, k
        ld = la + l - 3
        le = lb + l - 4
        c(ld) = c(le) - h(k)*c(le+1)
        d(ld) = d(le) + h(l-1)*d(le+1)
        r(ld) = 1.0_dp/(h(k+1) - h(l-1))
      end do
    end do

    ss = 10.**(-ll)                ! (precisao simples, como no original)
    tp = 0.1_dp*dir
    if (tp/tf > 0.5_dp) tp = 0.5_dp*tf
    ncount = 0

    ! ---- primeira sequencia
1000 ns = 0
    ni = 6
    tm = 0.0_dp
    call force(x, f1)

    ! ---- inicio de cada sequencia
2000 do k = 1, nv
      g(1,k) = b(1,k) + d(1)*b(2,k) + d(2)*b(3,k) + d(4)*b(4,k) + d(7)*b(5,k) + d(11)*b(6,k) + d(16)*b(7,k)
      g(2,k) = b(2,k) + d(3)*b(3,k) + d(5)*b(4,k) + d(8)*b(5,k) + d(12)*b(6,k) + d(17)*b(7,k)
      g(3,k) = b(3,k) + d(6)*b(4,k) + d(9)*b(5,k) + d(13)*b(6,k) + d(18)*b(7,k)
      g(4,k) = b(4,k) + d(10)*b(5,k) + d(14)*b(6,k) + d(19)*b(7,k)
      g(5,k) = b(5,k) + d(15)*b(6,k) + d(20)*b(7,k)
      g(6,k) = b(6,k) + d(21)*b(7,k)
      g(7,k) = b(7,k)
    end do

    t    = tp
    tval = abs(t)

    do m = 1, ni
      do j = 2, 8
        jd = j - 1
        s  = h(j)
        do k = 1, nv
          a = w(3)*b(3,k) + s*(w(4)*b(4,k) + s*(w(5)*b(5,k) + s*(w(6)*b(6,k) + s*w(7)*b(7,k))))
          y(k) = x(k) + t*s*(f1(k) + s*(w(1)*b(1,k) + s*(w(2)*b(2,k) + s*a)))
        end do
        call force(y, fj)

        do k = 1, nv
          temp = g(jd,k)
          gk   = (fj(k) - f1(k))/s
          if (j == 2) then
            g(1,k) = gk
          else
            gk = (gk - g(1,k))*r(nw(jd)+1)
            do l = 2, jd-1
              gk = (gk - g(l,k))*r(nw(jd)+l)
            end do
            g(jd,k) = gk
          end if
          temp = g(jd,k) - temp
          b(jd,k) = b(jd,k) + temp
          do l = 1, j-2
            b(l,k) = b(l,k) + c(nw(j-1)+l)*temp
          end do
        end do
      end do

      hv = 0.0_dp
      do k = 1, nv
        hv = max(hv, abs(b(7,k)))
      end do
      hv = hv*w(7)/(tval*tval*tval*tval*tval*tval*tval)
    end do

    ! controle do tamanho do passo na primeira sequencia
    if (.not. nsf) then
      tp = (ss**pw)/(hv**pw)*dir
      if (tp/t <= 1.0_dp) then
        tp = 0.8_dp*tp
        ncount = ncount + 1
        if (ncount > 10) return
        go to 1000
      end if
      nsf = .true.
    end if

    do k = 1, nv
      x(k) = x(k) + t*(f1(k) + b(1,k)*w(1) + b(2,k)*w(2) + b(3,k)*w(3) + b(4,k)*w(4) &
                     + b(5,k)*w(5) + b(6,k)*w(6) + b(7,k)*w(7))
    end do
    tm = tm + t
    ns = ns + 1

    if (nper) return
    call force(x, f1)

    tp = dir*(ss**pw)/(hv**pw)
    if (tp/t > sr) tp = t*sr
    if (dir*(tm+tp) >= dir*tf - 1.0e-8_dp) then
      tp   = tf - tm
      nper = .true.
    end if

    ! preve os novos B para a proxima sequencia
    q = tp/t
    do k = 1, nv
      if (ns /= 1) bd(:,k) = b(:,k) - e(:,k)
      e(1,k) = q*(b(1,k) + 2.0_dp*b(2,k) + 3.0_dp*b(3,k) + 4.0_dp*b(4,k) + 5.0_dp*b(5,k) &
                  + 6.0_dp*b(6,k) + 7.0_dp*b(7,k))
      e(2,k) = q**2*(b(2,k) + 3.0_dp*b(3,k) + 6.0_dp*b(4,k) + 10.0_dp*b(5,k) &
                  + 15.0_dp*b(6,k) + 21.0_dp*b(7,k))
      e(3,k) = q**3*(b(3,k) + 4.0_dp*b(4,k) + 10.0_dp*b(5,k) + 20.0_dp*b(6,k) + 35.0_dp*b(7,k))
      e(4,k) = q**4*(b(4,k) + 5.0_dp*b(5,k) + 15.0_dp*b(6,k) + 35.0_dp*b(7,k))
      e(5,k) = q**5*(b(5,k) + 6.0_dp*b(6,k) + 21.0_dp*b(7,k))
      e(6,k) = q**6*(b(6,k) + 7.0_dp*b(7,k))
      e(7,k) = q**7*b(7,k)
      b(:,k) = e(:,k) + bd(:,k)
    end do

    ni = 2
    go to 2000
  end subroutine ra15

end module radau

!=======================================================================
module espectro
  use sistema
  implicit none
contains

  ! FFT complexa (Numerical Recipes), data(2*nn) = (re, im, re, im, ...)
  subroutine four1(data, nn, isign)
    integer,  intent(in)    :: nn, isign
    real(dp), intent(inout) :: data(2*nn)
    integer  :: i, istep, j, m, mmax, n
    real(dp) :: theta, wi, wpi, wpr, wr, wtemp, tempi, tempr

    n = 2*nn
    j = 1
    do i = 1, n, 2
      if (j > i) then
        tempr = data(j);  tempi = data(j+1)
        data(j) = data(i); data(j+1) = data(i+1)
        data(i) = tempr;   data(i+1) = tempi
      end if
      m = n/2
      do while (m >= 2 .and. j > m)
        j = j - m
        m = m/2
      end do
      j = j + m
    end do
    mmax = 2
    do while (n > mmax)
      istep = 2*mmax
      theta = pi2/(isign*mmax)
      wpr = -2.0_dp*sin(0.5_dp*theta)**2
      wpi = sin(theta)
      wr = 1.0_dp
      wi = 0.0_dp
      do m = 1, mmax, 2
        do i = m, n, istep
          j = i + mmax
          tempr = wr*data(j) - wi*data(j+1)
          tempi = wr*data(j+1) + wi*data(j)
          data(j)   = data(i) - tempr
          data(j+1) = data(i+1) - tempi
          data(i)   = data(i) + tempr
          data(i+1) = data(i+1) + tempi
        end do
        wtemp = wr
        wr = wr*wpr - wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
      end do
      mmax = istep
    end do
  end subroutine four1

  ! Espectro da serie e seus picos. Devolve periodo de cada pico e a
  ! amplitude relativa (em % da maior amplitude com periodo > pmin).
  subroutine picos(serie, nn, dt, pmin, per, rel, np)
    integer,  intent(in)  :: nn
    real(dp), intent(in)  :: serie(nn), dt, pmin
    real(dp), allocatable, intent(out) :: per(:), rel(:)
    integer,  intent(out) :: np
    real(dp), allocatable :: data(:), amp(:)
    real(dp) :: ampmax
    integer  :: k

    allocate(data(2*nn), amp(0:nn/2+1))
    data(1:2*nn:2) = serie
    data(2:2*nn:2) = 0.0_dp
    call four1(data, nn, 1)

    amp = 0.0_dp
    do k = 1, nn/2
      amp(k) = 2.0_dp*sqrt(data(2*k+1)**2 + data(2*k+2)**2)/nn
    end do

    ! picos = maximos locais do espectro (frequencia k/(nn*dt))
    allocate(per(nn/2), rel(nn/2))
    np = 0
    ampmax = -1.0_dp
    do k = 1, nn/2
      if (amp(k) > amp(k-1) .and. amp(k+1) < amp(k)) then
        np = np + 1
        per(np) = 1.0_dp/(real(k,dp)/(nn*dt))
        rel(np) = amp(k)
        if (per(np) > pmin .and. amp(k) > ampmax) ampmax = amp(k)
      end if
    end do
    rel(1:np) = rel(1:np)/(ampmax/100.0_dp)
  end subroutine picos

end module espectro

!=======================================================================
program ips
  use sistema
  use dinamica
  use radau
  use espectro
  implicit none

  integer, parameter :: nser = 5
  character(len=4), parameter :: nomeser(nser) = ['res1', 'res2', 'inc ', 'a   ', 'e   ']

  type :: resultado
    logical :: pronto = .false.
    real(dp) :: a_km, e
    integer  :: np(nser)
    real(dp), allocatable :: per(:,:), rel(:,:)
  end type

  character(len=256) :: arqin, arqpla, linha
  character(len=20), allocatable :: nome(:)
  real(dp), allocatable :: el0(:,:), lim(:)
  type(resultado), allocatable :: res(:)
  integer, allocatable :: un(:,:)
  real(dp) :: req, mpla, razao, dt, a1, a2, e1, e2, pmin, ua, esc, grav
  integer  :: pot, nn, ll, ivar, na, ne, ia, ib, iobs, nlim, kres(2,6)
  integer  :: i, k, npts, iprox, u, ulog

  ! ------------------------------------------------ leitura do ips.in
  arqin = 'ips.in'
  if (command_argument_count() >= 1) call get_command_argument(1, arqin)
  open(newunit=u, file=arqin, status='old', action='read')
  call prox(u, linha); arqpla = primeira_palavra(linha)
  call prox(u, linha); read(linha, *) req
  call prox(u, linha); read(linha, *) mpla
  call prox(u, linha); read(linha, *) razao
  call prox(u, linha); read(linha, *) dj2, dj4
  call prox(u, linha); read(linha, *) dt
  call prox(u, linha); read(linha, *) pot
  call prox(u, linha); read(linha, *) ll
  call prox(u, linha); read(linha, *) ivar
  call prox(u, linha); read(linha, *) a1, a2, na
  call prox(u, linha); read(linha, *) e1, e2, ne
  call prox(u, linha); read(linha, *) ia, ib
  call prox(u, linha); read(linha, *) kres(1,:)
  call prox(u, linha); read(linha, *) kres(2,:)
  call prox(u, linha); read(linha, *) iobs
  call prox(u, linha); read(linha, *) pmin
  call prox(u, linha); read(linha, *) nlim
  allocate(lim(nlim))
  read(linha, *) nlim, lim
  close(u)

  ! ------------------------------------------------ leitura do planetas.in
  open(newunit=u, file=arqpla, status='old', action='read')
  N = 0
  do
    call prox(u, linha, fim=k)
    if (k /= 0) exit
    N = N + 1
  end do
  rewind(u)
  neq = 6*N
  allocate(nome(N), el0(N,6), dm(N), ami(N), dmimi(N))

  ! G em (raio do planeta)^3 / (massa do planeta * dia^2), como no original
  ua   = 1.49597870691e11_dp                ! metros
  esc  = ua/(req*1000.0_dp)
  grav = 4.0_dp*pi*pi*esc*esc*esc
  G    = grav/(razao*365.25_dp*365.25_dp)

  do i = 1, N
    call prox(u, linha)
    nome(i) = primeira_palavra(linha)
    k = len_trim(primeira_palavra(linha))
    read(linha(k+1:), *) dm(i), el0(i,:)
    dm(i)    = dm(i)/mpla                   ! massa em unidades do planeta
    el0(i,1) = el0(i,1)/req                 ! a em raios do planeta
    el0(i,3:6) = el0(i,3:6)*conv            ! graus -> radianos
    ami(i)   = G*(dm0 + dm(i))
    dmimi(i) = (dm0 + dm(i))/(dm0*dm(i))
  end do
  close(u)

  if (any([ivar, ia, ib, iobs] < 1) .or. any([ivar, ia, ib, iobs] > N)) &
    stop 'ips.in: indice de corpo fora de 1..N'

  nn = 2**pot
  npts = (na+1)*(ne+1)

  ! mensagens vao para ips.log (nada na tela; da pra rodar com &)
  open(newunit=ulog, file='ips.log', status='replace', action='write')
  write(ulog,'(a)') ' Inicio: '//agora()
  write(ulog,'(a,i0,a)')   ' Corpos: ', N, '  ('//trim(juntar(nome))//')'
  write(ulog,'(a,a)')      ' Corpo variado: ', trim(nome(ivar))
  write(ulog,'(a,i0,a,f0.2,a)') ' Pontos por integracao: ', nn, '  (', nn*dt/365.25_dp, ' anos)'
  write(ulog,'(a,i0,a,i0,a)')   ' Grade: ', na+1, ' valores de a x ', ne+1, ' valores de e'

  ! ------------------------------------------------ arquivos de saida
  allocate(un(nser, nlim))
  do k = 1, nser
    do i = 1, nlim
      open(newunit=un(k,i), file=trim(nomeser(k))//'_'//trim(fmtlim(lim(i)))//'.dat', &
           status='replace', action='write')
    end do
  end do

  ! ------------------------------------------------ varredura em a e e
  allocate(res(npts))
  iprox = 1
  !$omp parallel do schedule(dynamic,1)
  do k = 1, npts
    call rodar_ponto(k)
  end do
  !$omp end parallel do

  do k = 1, nser
    do i = 1, nlim
      close(un(k,i))
    end do
  end do
  write(ulog,'(a)') ' Fim: '//agora()
  close(ulog)

contains

  ! Integra um ponto da grade e faz a FFT das 5 series
  subroutine rodar_ponto(k)
    integer, intent(in) :: k
    real(dp), allocatable :: s(:,:), pk(:), rl(:)
    real(dp) :: el(N,6), x(neq)
    integer  :: j, js, np

    el = el0
    el(ivar,1) = a1 + (a2-a1)*real(mod(k-1, na+1), dp)/max(na,1)
    el(ivar,2) = e1 + (e2-e1)*real((k-1)/(na+1), dp)/max(ne,1)
    res(k)%a_km = el(ivar,1)
    res(k)%e    = el(ivar,2)
    el(ivar,1)  = el(ivar,1)/req

    allocate(s(nn, nser))
    call elem_para_estado(el, x)
    call estado_para_elem(x, el)
    call amostra(el, s(1,:))
    do j = 2, nn
      call ra15(x, dt, ll, neq)
      call estado_para_elem(x, el)
      call amostra(el, s(j,:))
    end do

    allocate(res(k)%per(nn/2, nser), res(k)%rel(nn/2, nser))
    do js = 1, nser
      call picos(s(:,js), nn, dt, pmin, pk, rl, np)
      res(k)%np(js) = np
      res(k)%per(1:np, js) = pk(1:np)
      res(k)%rel(1:np, js) = rl(1:np)
    end do
    deallocate(s)

    !$omp critical
    res(k)%pronto = .true.
    call escrever_prontos()
    !$omp end critical
  end subroutine rodar_ponto

  ! Series que vao para a FFT: res1, res2, I, a, e
  subroutine amostra(el, v)
    real(dp), intent(in)  :: el(N,6)
    real(dp), intent(out) :: v(nser)
    v(1) = angulo_res(el, kres(1,:))
    v(2) = angulo_res(el, kres(2,:))
    v(3) = el(iobs,3)
    v(4) = el(iobs,1)
    v(5) = el(iobs,2)
  end subroutine amostra

  ! phi = k1*lamA + k2*lamB + k3*varpiA + k4*varpiB + k5*OmA + k6*OmB
  real(dp) function angulo_res(el, kc)
    real(dp), intent(in) :: el(N,6)
    integer,  intent(in) :: kc(6)
    real(dp) :: lama, lamb, pia, pib
    pia  = el(ia,4) + el(ia,5)
    pib  = el(ib,4) + el(ib,5)
    lama = el(ia,6) + pia
    lamb = el(ib,6) + pib
    angulo_res = mod2pi(kc(1)*lama + kc(2)*lamb + kc(3)*pia + kc(4)*pib &
                      + kc(5)*el(ia,5) + kc(6)*el(ib,5))
  end function angulo_res

  ! Escreve, em ordem, os pontos da grade que ja terminaram
  subroutine escrever_prontos()
    integer :: js, il, j
    do while (iprox <= npts)
      if (.not. res(iprox)%pronto) exit
      associate (r => res(iprox))
        do js = 1, nser
          do il = 1, nlim
            do j = 1, r%np(js)
              if (r%rel(j,js) > lim(il)) &
                write(un(js,il), '(f14.4,1x,es14.6,1x,es16.8)') r%a_km, r%e, r%per(j,js)
            end do
            flush(un(js,il))
          end do
        end do
        write(ulog,'(a,i0,a,i0,a,f12.4,a,es12.5,a)') ' ponto ', iprox, '/', npts, &
             '   a(km) =', r%a_km, '   e =', r%e, '   '//agora()
        flush(ulog)
        deallocate(r%per, r%rel)
      end associate
      iprox = iprox + 1
    end do
  end subroutine escrever_prontos

  ! Le a proxima linha util (pula vazias e comentarios com ! ou #)
  subroutine prox(u, linha, fim)
    integer, intent(in) :: u
    character(len=*), intent(out) :: linha
    integer, intent(out), optional :: fim
    integer :: ios
    do
      read(u, '(a)', iostat=ios) linha
      if (ios /= 0) then
        if (present(fim)) then
          fim = ios
          return
        end if
        write(*,*) 'Erro: fim inesperado do arquivo de entrada.'
        stop 1
      end if
      linha = adjustl(replace_tab(linha))
      if (len_trim(linha) == 0) cycle
      if (linha(1:1) == '!' .or. linha(1:1) == '#') cycle
      exit
    end do
    if (present(fim)) fim = 0
  end subroutine prox

  ! data e hora, p/ o ips.log
  function agora() result(s)
    character(len=19) :: s
    integer :: v(8)
    call date_and_time(values=v)
    write(s,'(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') v(1:3), v(5:7)
  end function

  function replace_tab(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: k
    t = s
    do k = 1, len(t)
      if (t(k:k) == char(9)) t(k:k) = ' '
    end do
  end function

  function primeira_palavra(s) result(p)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: p
    integer :: k
    p = adjustl(s)
    k = index(p, ' ')
    if (k > 0) p = p(1:k-1)
  end function

  function juntar(v) result(s)
    character(len=*), intent(in) :: v(:)
    character(len=512) :: s
    integer :: k
    s = v(1)
    do k = 2, size(v)
      s = trim(s)//' '//trim(v(k))
    end do
  end function

  ! 0.5 -> '0.5', 8.0 -> '8'
  function fmtlim(x) result(s)
    real(dp), intent(in) :: x
    character(len=16) :: s
    integer :: k
    write(s, '(f16.3)') x
    s = adjustl(s)
    k = len_trim(s)
    do while (s(k:k) == '0')
      s(k:k) = ' '
      k = k - 1
    end do
    if (s(k:k) == '.') s(k:k) = ' '
  end function

end program ips
