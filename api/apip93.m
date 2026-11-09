c----------------------------------------------------------------------c

      subroutine readpost (fname)
      implicit none
      character*(*) fname
Use(P93dat)   # atn,atw,nt,nr,nn

c     local variables --
      integer ios, nget

c     procedures --
      external xerrab, remark, gallot, readpost1

c----------------------------------------------------------------------c
c     Read impurity emissivity and charge state from POST93 data files
c----------------------------------------------------------------------c

      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark("**** data file mist.dat not found --")
         call remark(" ")
         call remark("**** Data files for various impurities are available;")
         call remark("**** check uedge/in/api or contact authors")
         call remark(" ")
         call remark("**** For UEDGE, the data file must be re-named mist.dat")
         call xerrab("")
      endif

c     read array dimensions
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9001) atn, atw, nt, nr, nn

c     allocate storage for arrays
      call gallot("P93dat",0)

c     read data arrays
      call readpost1 (nget)

      close (nget)

 9000 format()
 9001 format(5(1x,i2/))

      return
      end

c----------------------------------------------------------------------c

      subroutine readpost1 (nget)
      implicit none
      integer nget
Use(P93dat)   # nt,nr,nn,tdatm,rdatm,ndatm,emdatm,z1datm,z2datm
Use(Physical_constants2)  # ev2

c     local variables --
      integer jn, jr, jt

      do jn=1,nn
         do jr=1,nr
            do jt=1,nt
               read(nget,9010) tdatm(jt,jr,jn), rdatm(jt,jr,jn), 
     &                         ndatm(jt,jr,jn), emdatm(jt,jr,jn), 
     &                         z1datm(jt,jr,jn), z2datm(jt,jr,jn)
            enddo
         enddo
      enddo

c     convert to SI units:
      do jt=1,nt
         do jr=1,nr
            do jn=1,nn
               tdatm(jt,jr,jn)=tdatm(jt,jr,jn)*ev2
               ndatm(jt,jr,jn)=ndatm(jt,jr,jn)*1.0e+06
               emdatm(jt,jr,jn)=emdatm(jt,jr,jn)*1.0e-06
            enddo
         enddo
      enddo

 9010 format(6(1x,e12.5))

      return
      end

c----------------------------------------------------------------------c

      subroutine splinem
      implicit none
Use(P93dat)    # nt,nr,nn
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata
      external gallot, splinem1

c     Construct 3-dimensional B-spline representation for impurity
c     emissivity and charge state versus e-temperature, ng/ne and
c     n*tau (data from POST '93 tables)

c     Allocate arrays for spline fitting --
      nxdata_api=nt		# temperature
      nydata_api=nr		# density ratio
      nzdata=nn		# n*tau
      nwork2 = kyords_api*kzords + 3*max(kxords_api,kyords_api,kzords) + kzords + 2
      nwork3 = nxdata_api*nydata_api*nzdata + 2*max( kxords_api*(nxdata_api+1),
     &           kyords_api*(nydata_api+1), kzords*(nzdata+1) )
      call gallot("Imslwrk",0)

      call splinem1

      return
      end

c----------------------------------------------------------------------c

      subroutine splinem1
      implicit none
Use(P93dat)    # tdatm,rdatm,ndatm,emdatm,z1datm,z2datm
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,fdata_api,ldf_api,mdf,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,emcoef,
               # z1coef,z2coef
      integer i,j,k
      external B3INT

c     Define data arrays --
      do i=1,nxdata_api
        xdata_api(i)=log10( tdatm(i,1,1) )
      enddo
      do j=1,nydata_api
        ydata_api(j)=log10( rdatm(1,j,1) )
      enddo
      do k=1,nzdata
        zdata(k)=log10( ndatm(1,1,k) )
      enddo

c     Define the order of the spline fit
c      kxords_api=4		# cubic in x=log10(temperature)
c      kyords_api=4		# cubic in y=log10(density ratio)
c      kzords=4		# cubic in z=log10(n*tau)

      ldf_api=nxdata_api
      mdf=nydata_api
      iflagi = 1	# let B3INT choose knots

c     Compute the coefficients --
c     first, for emissivity:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=log10( emdatm(i,j,k) )
            emcoef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            emcoef, ldf_api, mdf, work3, iflagi)
c     next, for average Z:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=z1datm(i,j,k)
            z1coef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            z1coef, ldf_api, mdf, work3, iflagi)
c     next, for average Z**2:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=z2datm(i,j,k)
            z2coef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            z2coef, ldf_api, mdf, work3, iflagi)

      return
      end

c----------------------------------------------------------------------c

      subroutine loadpostmul (nzsp, fnames, inew)    #..zml
      implicit none
      integer nzsp, inew
      character*(*) fnames(*)
Use(P93dat)          # atn
Use(Multi_fix_frc)   # atn_mul
Use(Multi_fix_frc_ld)  # ipostld,nzspffld,fnameld

c     Read and spline the Post fixed-fraction tables for nzsp species,
c     but only when needed: tables are kept in memory between exmain
c     calls instead of being re-read and re-splined every time.
c       inew = 1: reload if nothing is loaded yet, or nzsp, a file name,
c                 or a file's atomic number has changed since the last load
c       inew = 2: always reload (previous behaviour)
c       inew = 0: never reload; tables must already be in memory

      integer is, iatn, ireload
      character*256 msg
      external xerrab, remark, readpost, splinem, savemulffparam
      external peekpostatn

c     fnameld is dimensioned nzspffmx in api.v
      if (nzsp .lt. 1 .or. nzsp .gt. size(fnameld)) then
         write(msg,'(a,i3,a,i3)') '*** loadpostmul: nzspff =', nzsp,
     .        ' must be between 1 and nzspffmx =', size(fnameld)
         call xerrab(trim(msg))
      endif

      if (inew .eq. 0) then
         if (ipostld .eq. 0 .or. nzspffld .ne. nzsp) call xerrab(
     .     '*** loadpostmul: newpost=0 but tables for nzspff not loaded')
         return
      endif

      ireload = 0
      if (inew .eq. 2 .or. ipostld .eq. 0 .or. nzspffld .ne. nzsp)
     .   ireload = 1
      do is = 1, nzsp
         if (ireload .eq. 0) then
            if (trim(fnames(is)) .ne. trim(fnameld(is))) ireload = 1
         endif
         if (ireload .eq. 0) then
c           same name, but the file may differ (e.g. another case dir)
            call peekpostatn(fnames(is), iatn)
            if (iatn .ne. atn_mul(is)) ireload = 1
         endif
      enddo
      if (ireload .eq. 0) return

      do is = 1, nzsp
         call readpost(fnames(is))
         call splinem
         call savemulffparam(is)
         fnameld(is) = fnames(is)
         write(msg,'(a,i2,a,i3,2a)') '*** fixed-fraction species', is,
     .        ': atomic number', atn, ' from ', trim(fnames(is))
         call remark(trim(msg))
      enddo
      nzspffld = nzsp
      ipostld = 1

      return
      end

c----------------------------------------------------------------------c

      subroutine peekpostatn (fname, iatn)    #..zml
      implicit none
      character*(*) fname
      integer iatn

c     Return the atomic number in the header of a Post table without
c     reading the table itself.

      integer ios, nget, iatw, intt, inr, inn
      external xerrab
      data nget /56/

      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) call xerrab('*** peekpostatn: cannot open '
     .                            //trim(fname))
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9001) iatn, iatw, intt, inr, inn
      close (nget)

 9000 format()
 9001 format(5(1x,i2/))

      return
      end

c----------------------------------------------------------------------c

      subroutine savemulffparam (izspff)    #..zml
      implicit none
      integer izspff
Use(P93dat)   # nt,nr,nn,tdatm,rdatm,ndatm,emdatm,z1datm,z2datm
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,fdata_api,ldf_api,mdf,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,emcoef,
               # z1coef,z2coef
Use(Multi_fix_frc)   # atn_mul,atw_mul,emdatm_mul,z1datm_mul,z2datm_mul,
                     # emcoef_mul,z1coef_mul,z2coef_mul,tgrid_mul,rgrid_mul,ngrid_mul
Use(Multi_fix_frc_ld)  # ntld,nrld,nnld
      integer i
      real reltol
      character*256 msg
      external gallot, xerrab
      data reltol /1.e-6/

c     All species share one set of spline knots (those of the last table
c     splined), so every table must have species 1's grid.  Check before
c     storing anything.
      if (izspff .lt. 2) then
         call gallot("Multi_fix_frc",0)
         ntld = nt
         nrld = nr
         nnld = nn
         do i = 1, nt
            tgrid_mul(i) = tdatm(i,1,1)
         enddo
         do i = 1, nr
            rgrid_mul(i) = rdatm(1,i,1)
         enddo
         do i = 1, nn
            ngrid_mul(i) = ndatm(1,1,i)
         enddo
      else
         if (nt .ne. ntld .or. nr .ne. nrld .or. nn .ne. nnld) then
            write(msg,'(a,i2,a,3i4,a,3i4)') '*** savemulffparam: species',
     .        izspff, ' table size', nt, nr, nn,
     .        ' differs from species 1', ntld, nrld, nnld
            call xerrab(trim(msg))
         endif
         do i = 1, nt
            if (abs(tdatm(i,1,1)-tgrid_mul(i)) .gt.
     .          reltol*abs(tgrid_mul(i))) call xerrab(
     .        '*** savemulffparam: temperature grid differs from species 1')
         enddo
         do i = 1, nr
            if (abs(rdatm(1,i,1)-rgrid_mul(i)) .gt.
     .          reltol*abs(rgrid_mul(i))) call xerrab(
     .        '*** savemulffparam: n0/ne grid differs from species 1')
         enddo
         do i = 1, nn
            if (abs(ndatm(1,1,i)-ngrid_mul(i)) .gt.
     .          reltol*abs(ngrid_mul(i))) call xerrab(
     .        '*** savemulffparam: n*tau grid differs from species 1')
         enddo
      endif

      atn_mul(izspff) = atn
      atw_mul(izspff) = atw
      emdatm_mul(:,:,:,izspff) = emdatm(:,:,:)
      z1datm_mul(:,:,:,izspff) = z1datm(:,:,:)
      z2datm_mul(:,:,:,izspff) = z2datm(:,:,:)
      emcoef_mul(:,:,:,izspff) = emcoef(:,:,:)
      z1coef_mul(:,:,:,izspff) = z1coef(:,:,:)
      z2coef_mul(:,:,:,izspff) = z2coef(:,:,:)

      return
      end

c-----------------------------------------------------------------------

      real function emissbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,emcoef

      integer nxcoef,nycoef,nzcoef
      real vlogw,w,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     emissivity [Watts-m**3] versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vlogw=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             emcoef, ldf_api, mdf, icont, iworki, work2, iflagi)
      w=10**vlogw
      emissbs=w

      return
      end

c-----------------------------------------------------------------------

      real function emissbs_mul (vte,vnr,vnt,izspff)    #..zml
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots
Use(Multi_fix_frc) # emcoef_mul

      integer nxcoef,nycoef,nzcoef,izspff
      real vlogw,w,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     emissivity [Watts-m**3] versus :
c               vte = e-temperature [J]
c               vnr = ng/ne density ratio
c               vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      #emcoef(:,:,:) = emcoef_mul(:,:,:,izspff)

      icont = 0
      vlogw=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords,
     &             emcoef_mul(:,:,:,izspff), ldf_api, mdf, icont, iworki, work2, iflagi)
      w=10**vlogw
      emissbs_mul=w

      return
      end

c----------------------------------------------------------------------c

      real function z1avgbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,z1coef

      integer nxcoef,nycoef,nzcoef
      real vz1,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     charge state Z versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vz1=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             z1coef, ldf_api, mdf, icont, iworki, work2, iflagi)
      z1avgbs=vz1

      return
      end

c-----------------------------------------------------------------------

      real function z2avgbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,z2coef

      integer nxcoef,nycoef,nzcoef
      real vz2,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     charge state Z**2 versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vz2=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             z2coef, ldf_api, mdf, icont, iworki, work2, iflagi)
      z2avgbs=vz2

      return
      end

