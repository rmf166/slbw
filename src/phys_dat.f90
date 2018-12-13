      module physics

      ! Provides pi and physics constants taken from CODATA'89 as
      ! given on the NIST site; namely, bk (Boltzmann's constant),
      ! amassn (the neutron mass in amu), amu (the amu value itself),
      ! hbar (Planck's constant), ev (the conversion to eV), and
      ! clight (the speed of light).

        implicit none

        real(8),    parameter   :: pi=3.14159265358979d0
        real(8),    parameter   :: bk=8.617385d-05        ! eV/K
        real(8),    parameter   :: amassn=1.008664904d0   ! u
        real(8),    parameter   :: amassu=238.050788423d0 ! u
        real(8),    parameter   :: amu=1.6605402d-24      ! g/u
        real(8),    parameter   :: hbar=1.05457266d-27    ! g-cm**2/s
        real(8),    parameter   :: ev=1.60217733d-12      ! g-cm**2/s**2
        real(8),    parameter   :: clight=2.99792458d+10

      end module physics
