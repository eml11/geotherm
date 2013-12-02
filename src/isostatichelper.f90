      module isostatichelper
      use netcdf
      use mathmodule
      use modeldomainmodule
      implicit none

      contains
      subroutine get_basedensity(filename,density,domain)
      type (modeldomain) domain
      double precision :: density(:,:)
      double precision :: data_ar(domain%n,domain%m)
      integer i,ncid,varid,retval
      character(len=256) :: filename
      
      density = 0D0

      !print *, domain%minerals
      !print *, domain%minerals
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      do i=1,domain%minerals
        retval = nf90_inq_varid(ncid, &
        &domain%mineralarray(i)%mineralname, varid)
        retval = nf90_get_var(ncid, varid, data_ar)
        !print *, data_ar
        !print *, data_ar
        where (data_ar.GT.1)
          data_ar=0
        elsewhere
          data_ar=data_ar
        end where
        density = density + data_ar*domain%mineralarray(i)%density
        !print *, MAXVAL(domain%mineralarray(i)%density)
        !print *, maxval(data_ar)
      enddo
      !print *, density
      retval = nf90_close(ncid)

      end subroutine

      subroutine domaindensityintegral(density,domain,retrn_ar )
      type (modeldomain) domain
      double precision :: retrn_ar(domain%n,domain%m)
      double precision :: cnst(domain%n,domain%m)
      double precision :: density(domain%n,domain%m)

      call array_integral2dydim(density,retrn_ar, &
      &domain%incriment(2),domain%n,domain%m)

      call extend_ardimension(retrn_ar(:,1),cnst,domain%m)
      retrn_ar = retrn_ar - cnst

      end subroutine

      subroutine mineraldensityintegral(density,domain,retrn_ar)
      type (modeldomain) domain
      double precision :: retrn_ar(domain%n,domain%m)
      double precision :: cnst(domain%n,domain%m)
      double precision :: density(domain%n,domain%m)

      call array_integral2dydim(density,retrn_ar, &
      &domain%incriment(2),domain%n,domain%m)

      call extend_ardimension(retrn_ar(:,1),cnst,domain%m)
      retrn_ar = retrn_ar - cnst

      end subroutine

      subroutine compute_deltay(domain, &
      &integral1,integral2,retrn_ar)
      type (modeldomain) domain
      double precision :: density_0(domain%n)
      double precision :: integral1(:)
      double precision :: integral2(:)
      double precision :: retrn_ar(:)
      double precision :: bulkmodulus(domain%n)
      double precision, parameter :: g = 9.81     
 
      bulkmodulus = domain%bulkmodulus(:,domain%m)
      density_0 = domain%density(:,domain%m)

      retrn_ar = ((bulkmodulus)/(density_0*g)) - &
      &((-integral2)/(density_0)) - &
      &((bulkmodulus)/(density_0*g)) * &
      &(1 + ((integral1*g)/(bulkmodulus))) * &
      &10 ** ((-integral2 + integral1)/(bulkmodulus))


      end subroutine

      subroutine write_output(outfile,deltay)
      character(len = 256) outfile
      double precision :: deltay(:)
      integer i,sz

      OPEN(1,file = outfile)
      sz = SIZE(deltay)

      WRITE(1,*) "possition,topographic_change(m)"

      do i=1,sz
        WRITE(1,*) i,",",deltay(1)
      enddo

      end subroutine

      end module
