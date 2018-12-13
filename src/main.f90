      program ssph

        implicit none

        call drive

      contains

      subroutine drive

        implicit none

        integer(4)              :: r
        integer(4)              :: t
        integer(4), parameter   :: ng=10000
        integer(4), parameter   :: nt=7
        real(8)                 :: aj
        real(8)                 :: del
        real(8)                 :: eng0
        real(8)                 :: emin
        real(8)                 :: gamn
        real(8)                 :: gamg
        real(8)                 :: sigg(ng+1)
        real(8)                 :: siga(ng+1,nt)
        real(8)                 :: spi
        real(8)                 :: tmp(nt)

      ! set temperatures

        tmp(1)=1.0d-06
        tmp(2)=100.0d0
        tmp(3)=300.0d0
        tmp(4)=500.0d0
        tmp(5)=1000.0d0
        tmp(6)=2000.0d0
        tmp(7)=10000.0d0

      ! plot psi line shape

        call pltpsi(nt,tmp)

      ! run single-level breit-wigner model

        siga=0.0d0
        do r=1,3
          do t=1,nt
            call respar(r,eng0,spi,aj,gamn,gamg)
            call seteng(r,ng,eng0,emin,del)
            call reconr(ng,del,eng0,emin,spi,aj,gamn,gamg,tmp(t),sigg)
            siga(:,t)=sigg(:)
          enddo
          call pltsig(r,ng,nt,del,emin,siga)
        enddo

      end subroutine drive

      subroutine pltpsi(nt,tmp)

        use physics

        implicit none

        integer(4), intent(in)  :: nt
        real(8),    intent(in)  :: tmp(nt)

        integer(4)              :: t
        integer(4)              :: n
        integer(4), parameter   :: nx=10000
        real(8)                 :: aimw
        real(8)                 :: awri
        real(8)                 :: ax
        real(8)                 :: delta
        real(8)                 :: dx
        real(8)                 :: eng
        real(8)                 :: eng0
        real(8)                 :: gamt
        real(8)                 :: psi(nx+1,nt)
        real(8)                 :: rew
        real(8)                 :: rpi
        real(8)                 :: tbk
        real(8)                 :: theta
        real(8)                 :: x
        real(8)                 :: xmin
        real(8)                 :: xmax
        real(8)                 :: y
        character(7)            :: datafile
        character(10)           :: plotfile
        character(12)           :: pdf_file
        logical(4)              :: flag

        rpi=dsqrt(pi)

        xmin=-40.0d0
        xmax= 40.0d0
        dx=(xmax-xmin)/dble(nx)

      ! set a few physical parameters

        gamt=0.023d0
        eng0=6.673491d+0
        awri=amassu/amassn

      ! loop over temp. and x-values

        psi=0.0d0
        do t=1,nt
          do n=1,nx+1
            x=xmin+(n-1)*dx
            eng=eng0+x*0.5d0*gamt
            tbk=tmp(t)*bk
            delta=dsqrt(4.0d0*tbk*eng/awri)
            theta=gamt/delta
            ax=theta*x/2.0d0
            y=theta/2.0d0
            call wofz(ax,y,rew,aimw,flag)
            if (flag) then
              write(0,*) ' Overflow occured during wofz call.'
              stop
            endif
            psi(n,t)=rpi*theta*rew/2.0d0
          enddo
        enddo

        datafile='psi.dat'
        open(unit=1,file=datafile,action='write',status='unknown')
        do n=1,nx+1
          x=xmin+(n-1)*dx
          write(1,'(8(es12.5))') x,(psi(n,t),t=1,nt)
        enddo
        close(1)

        plotfile='plot-psi.p'
        pdf_file='plot-psi.pdf'
        open(unit=1,file=plotfile,action='write',status='unknown')
        write(1,'(a)') 'set autoscale'
        write(1,'(a)') 'unset logscale'
        write(1,'(a)') 'unset label'
        write(1,'(a)') 'set xtic auto'
        write(1,'(a)') 'set ytic auto'
        write(1,'(a)') 'set title "Psi vs x"'
        write(1,'(a)') 'set xlabel "x" enhanced'
        write(1,'(a)') 'set ylabel "psi" enhanced'
        write(1,'(a)') 'plot    "' // trim(datafile) // &
                       '" using 1:2 title "T=0" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:3 title "T=100" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:4 title "T=300" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:5 title "T=500" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:6 title "T=1000" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:7 title "T=2000" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:8 title "T=10000" with lines '
        write(1,'(a)') 'set terminal pdfcairo enhanced color'
        write(1,'(a)') 'set output "' // trim(pdf_file)  // '"'
        write(1,'(a)') 'replot'
        write(1,'(a)') 'set terminal x11'
        close(1)

      end subroutine pltpsi

      subroutine pltsig(r,ng,nt,del,emin,siga)

        implicit none

        integer(4), intent(in)  :: r
        integer(4), intent(in)  :: ng
        integer(4), intent(in)  :: nt
        real(8),    intent(in)  :: del
        real(8),    intent(in)  :: emin
        real(8),    intent(in)  :: siga(ng+1,nt)

        integer(4)              :: ig
        integer(4)              :: t
        real(8)                 :: eng
        character(1)            :: res
        character(8)            :: datafile
        character(7)            :: plotfile
        character(9)            :: pdf_file

        write(res,'(i1)') r
        datafile='res'//res//'.dat'
        open(unit=1,file=datafile,action='write',status='unknown')
        do ig=ng+1,1,-1
          eng=emin+(ig-1)*del
          write(1,'(8(es12.5))') eng,(siga(ig,t),t=1,nt)
        enddo
        eng=emin+ng*del
        close(1)

        plotfile='plot'//res//'.p'
        pdf_file='plot'//res//'.pdf'
        open(unit=1,file=plotfile,action='write',status='unknown')
        write(1,'(a)') 'set autoscale'
        write(1,'(a)') 'unset logscale'
        write(1,'(a)') 'unset label'
        write(1,'(a)') 'set xtic auto'
        write(1,'(a)') 'set ytic auto'
        if (r == 1) then
          write(1,'(a)') 'set title "U-238 resonance at 6.67 eV"'
        elseif (r == 2) then
          write(1,'(a)') 'set title "U-238 resonance at 20.87 eV"'
        elseif (r == 3) then
          write(1,'(a)') 'set title "U-238 resonance at 36.68 eV"'
        endif
        write(1,'(a)') 'set xlabel "Energy (eV)" enhanced'
        write(1,'(a)') 'set ylabel "siga (b)" enhanced'
        write(1,'(a,es12.5,a,es12.5,a)') 'set xr [',emin,':',eng,']'
        write(1,'(a)') 'plot    "' // trim(datafile) // &
                       '" using 1:2 title "T=0" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:3 title "T=100" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:4 title "T=300" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:5 title "T=500" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:6 title "T=1000" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:7 title "T=2000" with lines, \'
        write(1,'(a)') '        "' // trim(datafile) // &
                       '" using 1:8 title "T=10000" with lines '
        write(1,'(a)') 'set terminal pdfcairo enhanced color'
        write(1,'(a)') 'set output "' // trim(pdf_file)  // '"'
        write(1,'(a)') 'replot'
        write(1,'(a)') 'set terminal x11'
        close(1)

      end subroutine pltsig

      subroutine respar(r,eng0,spi,aj,gamn,gamg)

        implicit none

        integer(4), intent(in)  :: r
        real(8),    intent(out) :: eng0
        real(8),    intent(out) :: spi
        real(8),    intent(out) :: aj
        real(8),    intent(out) :: gamn
        real(8),    intent(out) :: gamg

        ! data from /work/CMS/libgen/xs_evals/endf7r1_data/neutrons/n-092_U_238.endf
        !
        ! 9.223800+4 2.360058+2          0          0          1          09237 2151    1
        ! 9.223800+4 1.000000+0          0          0          2          09237 2151    2
        ! 1.000000-5 2.000000+4          1          3          0          09237 2151    3
        ! 0.000000+0 9.480000-1          0          0          2          29237 2151    4
        ! 2.360060+2 9.480000-1          0          0       5556        9269237 2151    5
        ! ...
        ! 6.673491+0 5.000000-1 1.475792-3 2.300000-2 0.000000+0 9.990000-99237 2151   28
        ! 2.087152+1 5.000000-1 1.009376-2 2.286379-2 5.420000-8 0.000000+09237 2151   29
        ! 3.668212+1 5.000000-1 3.354568-2 2.300225-2 0.000000+0 9.770000-99237 2151   30
        !
        ! resonance range from 1.000000e-5 eV to 2.000000e+4 eV
        !   target spin is 0.000000e+0
        !   scatter length is 9.480000e-1 (potential scatter is 11.29340d0)
        ! 

        spi =0.0d0
        aj  =0.5d0

        if (r == 1) then
        ! 
        ! resonance at 6.673491e+0 eV with spin 5.000000e-1
        !   neutron width is 1.475792e-3 eV (gamn)
        !   gamma   width is 2.300000e-2 eV (gamg)
        ! 
          eng0=6.673491d+0
          gamn=1.475792d-3
          gamg=2.300000d-2
        elseif (r == 2) then
        ! 
        ! resonance at 2.087152e+1 eV with spin 5.000000e-1
        !   neutron width is 1.009376e-2 eV (gamn)
        !   gamma   width is 2.286379e-2 eV (gamg)
        ! 
          eng0=2.087152d+1
          gamn=1.009376d-2
          gamg=2.286379d-2
        elseif (r == 3) then
        ! 
        ! resonance at 3.668212e+1 eV with spin 5.000000e-1
        !   neutron width is 3.354568e-2 eV (gamn)
        !   gamma   width is 2.300225e-2 eV (gamg)
        !
          eng0=3.668212d+1
          gamn=3.354568d-2
          gamg=2.300225d-2
        endif

      end subroutine respar

      subroutine seteng(r,ng,eng0,emin,del)

        implicit none

        integer(4), intent(in)  :: r
        integer(4), intent(in)  :: ng
        real(8),    intent(in)  :: eng0
        real(8),    intent(out) :: emin
        real(8),    intent(out) :: del

        real(8)                 :: emax
        real(8)                 :: gamr
        real(8),    parameter   :: x=20.0d0

      ! total resonance widths from data

        if (r == 1) then
          gamr=1.475792d-3+2.300000d-2
        elseif (r == 2) then
          gamr=1.009376d-2+2.286379d-2
        elseif (r == 3) then
          gamr=3.354568d-2+2.300225d-2
        endif

        emin=eng0-0.5d0*gamr*x
        emax=eng0+0.5d0*gamr*x

        if (emin < 0.0d0 .or. emax < 0.0d0) then
          write(0,*) ' Emin or Emax less than zero. Check 2 eV spacing.'
          stop
        endif

        del=(emax-emin)/dble(ng)

      end subroutine seteng

      subroutine reconr(ng,del,eng0,emin,spi,aj,gamn,gamg,tmp,sigg)

        implicit none

        integer(4), intent(in)  :: ng
        real(8),    intent(in)  :: del
        real(8),    intent(in)  :: eng0
        real(8),    intent(in)  :: emin
        real(8),    intent(in)  :: spi
        real(8),    intent(in)  :: aj
        real(8),    intent(in)  :: gamn
        real(8),    intent(in)  :: gamg
        real(8),    intent(out) :: sigg(ng+1)
        real(8),    intent(in)  :: tmp

        integer(4)              :: ig
        real(8)                 :: eng

        do ig=ng+1,1,-1
          eng=emin+(ig-1)*del
          call slbw(eng,eng0,spi,aj,gamn,gamg,tmp,sigg(ig))
        enddo

      end subroutine reconr

      subroutine slbw(eng,eng0,spi,aj,gamn,gamg,tmp,sigg)

        use physics

        implicit none

        real(8),    intent(in)  :: eng   ! eV
        real(8),    intent(in)  :: eng0  ! eV
        real(8),    intent(in)  :: spi   ! spin i
        real(8),    intent(in)  :: aj    ! spin j
        real(8),    intent(in)  :: gamn
        real(8),    intent(in)  :: gamg
        real(8),    intent(in)  :: tmp   ! K
        real(8),    intent(out) :: sigg

        real(8)                 :: aimw
        real(8)                 :: arat
        real(8)                 :: awri
        real(8)                 :: ax
        real(8)                 :: delta
        real(8)                 :: e0
        real(8)                 :: ex
        real(8)                 :: gamt
        real(8)                 :: gj
        real(8)                 :: m0
        real(8)                 :: psi
        real(8)                 :: rew
        real(8)                 :: rpi
        real(8)                 :: r
        real(8)                 :: spifac
        real(8)                 :: tbk
        real(8)                 :: xi
        real(8)                 :: y
        logical(4)              :: flag
        logical(4)              :: use_njoy

      ! decide on slbw njoy implementation

        use_njoy=.false.

      ! set constants

        rpi=dsqrt(pi)

        awri=amassu/amassn
        arat=awri/(awri+1.0d0)

      ! calculate spin factor

        spifac=1.0d0/(2.0d0*spi+1.0d0)
        gj=(2.0d0*aj+1.0d0)*spifac/2.0d0

      ! convert neutron mass to reduced mass

        m0=amassn*arat

      ! convert neutron energy from lab to center of mass system
      ! NOTE: incident energy approximated by resonance energy

        if (use_njoy) then
          e0=eng *arat
        else
          e0=eng0*arat
        endif

      ! calculate r in barn-eV units

        r=gj*2.0d0*pi*hbar*hbar*1.0d+24/(m0*amu*ev*e0)

      ! calculate total energy level width

        gamt=gamn+gamg

      ! use Bethe-Placzek approximation for delta: eng~eng0

        ex=2.0d0*(eng-eng0)/gamt
        tbk=tmp*bk
        if (use_njoy) then
          delta=dsqrt(4.0d0*tbk*eng /awri)
        else
          delta=dsqrt(4.0d0*tbk*eng0/awri)
        endif
        xi=gamt/delta

      ! calculate doppler broadening

        ax=xi*ex/2.0d0
        y=xi/2.0d0
        call wofz (ax,y,rew,aimw,flag)
        if (flag) then
          write(0,*) ' Overflow occured during wofz call.'
          stop
        endif
        psi=xi*rpi*rew/2.0d0

      ! final resonance cross section
      ! dsqrt(eng0/eng) may be included or approximated as unity

        if (use_njoy) then
          sigg=gamn*gamg/(gamt*gamt)                *r*psi
        else
          sigg=gamn*gamg/(gamt*gamt)*dsqrt(eng0/eng)*r*psi
        endif

      end subroutine slbw

      subroutine wofz (xi, yi, u, v, flag)
      !
      !      algorithm 680, collected algorithms from acm.
      !      this work published in transactions on mathematical software,
      !      vol. 16, no. 1, pp. 47.
      !
      !  given a complex number z = (xi,yi), this subroutine computes
      !  the value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
      !  where erfc is the complex complementary error-function and i
      !  means sqrt(-1).
      !  the accuracy of the algorithm for z in the 1st and 2nd quadrant
      !  is 14 significant digits; in the 3rd and 4th it is 13 significant
      !  digits outside a circular region with radius 0.126 around a zero
      !  of the function.
      !  all real variables in the program are double precision.
      !
      !
      !  the code contains a few compiler-dependent parameters :
      !     rmaxreal = the maximum value of rmaxreal equals the root of
      !                rmax = the largest number which can still be
      !                implemented on the computer in double precision
      !                floating-point arithmetic
      !     rmaxexp  = ln(rmax) - ln(2)
      !     rmaxgoni = the largest possible argument of a double precision
      !                goniometric function (dcos, dsin, ...)
      !  the reason why these parameters are needed as they are defined will
      !  be explained in the code by means of comments
      !
      !
      !  parameter list
      !     xi     = real      part of z
      !     yi     = imaginary part of z
      !     u      = real      part of w(z)
      !     v      = imaginary part of w(z)
      !     flag   = an error flag indicating whether overflow will
      !              occur or not; type logical;
      !              the values of this variable have the following
      !              meaning :
      !              flag=.false. : no error condition
      !              flag=.true.  : overflow will occur, the routine
      !                             becomes inactive
      !  xi, yi      are the input-parameters
      !  u, v, flag  are the output-parameters
      !
      !  furthermore the parameter factor equals 2/sqrt(pi)
      !
      !  the routine is not underflow-protected but any variable can be
      !  put to 0 upon underflow;
      !
      !  reference - gpm poppe, cmj wijers; more efficient computation of
      !  the complex error-function, acm trans. math. software.
      !
      !
      !
      !
      !
        implicit none
      !
        integer(4)              :: i
        integer(4)              :: j
        integer(4)              :: kapn
        integer(4)              :: n
        integer(4)              :: np1
        integer(4)              :: nu
        real(8)                 :: c
        real(8)                 :: daux
        real(8),    parameter   :: factor=1.12837916709551257388d0
        real(8)                 :: h
        real(8)                 :: h2
        real(8)                 :: qlambda
        real(8)                 :: qrho
        real(8),    parameter   :: rmaxexp=708.503061461606d0
        real(8),    parameter   :: rmaxgoni=3.53711887601422d+15
        real(8),    parameter   :: rmaxreal=0.5d+154
        real(8)                 :: rx
        real(8)                 :: ry
        real(8)                 :: sx
        real(8)                 :: sy
        real(8)                 :: tx
        real(8)                 :: ty
        real(8)                 :: u
        real(8)                 :: u1
        real(8)                 :: u2
        real(8)                 :: v
        real(8)                 :: v1
        real(8)                 :: v2
        real(8)                 :: w1
        real(8)                 :: x
        real(8)                 :: xabs
        real(8)                 :: xabsq
        real(8)                 :: xaux
        real(8)                 :: xi
        real(8)                 :: xquad
        real(8)                 :: xsum
        real(8)                 :: y
        real(8)                 :: yabs
        real(8)                 :: yi
        real(8)                 :: yquad
        real(8)                 :: ysum
        logical(4)              :: a
        logical(4)              :: b
        logical(4)              :: flag
      !
        flag = .false.
      !
        xabs = dabs(xi)
        yabs = dabs(yi)
        x    = xabs/6.3
        y    = yabs/4.4
      !
      !
      !     the following if-statement protects
      !     qrho = (x**2 + y**2) against overflow
      !
        if ((xabs.gt.rmaxreal).or.(yabs.gt.rmaxreal)) goto 100
      !
        qrho = x**2 + y**2
      !
        xabsq = xabs**2
        xquad = xabsq - yabs**2
        yquad = 2*xabs*yabs
      !
        a     = qrho.lt.0.085264d0
      !
        if (a) then
      !
      !  if (qrho.lt.0.085264d0) then the faddeeva-function is evaluated
      !  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
      !  n is the minimum number of terms needed to obtain the required
      !  accuracy
      !
          qrho  = (1-0.85*y)*dsqrt(qrho)
          n     = idnint(6 + 72*qrho)
          j     = 2*n+1
          xsum  = 1.0/j
          ysum  = 0.0d0
          do 10 i=n, 1, -1
            j    = j - 2
            xaux = (xsum*xquad - ysum*yquad)/i
            ysum = (xsum*yquad + ysum*xquad)/i
            xsum = xaux + 1.0/j
   10     continue
          u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0
          v1   =  factor*(xsum*xabs - ysum*yabs)
          daux =  dexp(-xquad)
          u2   =  daux*dcos(yquad)
          v2   = -daux*dsin(yquad)
      !
          u    = u1*u2 - v1*v2
          v    = u1*v2 + v1*u2
      !
        else
      !
      !  if (qrho.gt.1.o) then w(z) is evaluated using the laplace
      !  continued fraction
      !  nu is the minimum number of terms needed to obtain the required
      !  accuracy
      !
      !  if ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
      !  by a truncated taylor expansion, where the laplace continued fraction
      !  is used to calculate the derivatives of w(z)
      !  kapn is the minimum number of terms in the taylor expansion needed
      !  to obtain the required accuracy
      !  nu is the minimum number of terms of the continued fraction needed
      !  to calculate the derivatives with the required accuracy
      !
      !
          if (qrho.gt.1.0) then
            h    = 0.0d0
            kapn = 0
            qrho = dsqrt(qrho)
            nu   = idint(3 + (1442/(26*qrho+77)))
          else
            qrho = (1-y)*dsqrt(1-qrho)
            h    = 1.88*qrho
            h2   = 2*h
            kapn = idnint(7  + 34*qrho)
            nu   = idnint(16 + 26*qrho)
          endif
      !
          b = (h.gt.0.0)
      !
          if (b) qlambda = h2**kapn
      !
          rx = 0.0
          ry = 0.0
          sx = 0.0
          sy = 0.0
      !
          do 11 n=nu, 0, -1
            np1 = n + 1
            tx  = yabs + h + np1*rx
            ty  = xabs - np1*ry
            c   = 0.5/(tx**2 + ty**2)
            rx  = c*tx
            ry  = c*ty
            if ((b).and.(n.le.kapn)) then
              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda/h2
            endif
   11     continue
      !
          if (h.eq.0.0) then
            u = factor*rx
            v = factor*ry
          else
            u = factor*sx
            v = factor*sy
          end if
      !
          if (yabs.eq.0.0) u = dexp(-xabs**2)
      !
        end if
      !
      !
      !
      !  evaluation of w(z) in the other quadrants
      !
      !
        if (yi.lt.0.0) then
      !
          if (a) then
            u2    = 2*u2
            v2    = 2*v2
          else
            xquad =  -xquad
      !
      !
      !   the following if-statement protects 2*exp(-z**2)
      !   against overflow
      !
            if ((yquad.gt.rmaxgoni).or. &
                (xquad.gt.rmaxexp)) goto 100
      !
            w1 =  2*dexp(xquad)
            u2  =  w1*dcos(yquad)
            v2  = -w1*dsin(yquad)
          end if
      !
          u = u2 - u
          v = v2 - v
          if (xi.gt.0.0) v = -v
        else
          if (xi.lt.0.0) v = -v
        end if
      !
        return
      !
    100 flag = .true.
        return
      !
      end subroutine wofz

      end program
