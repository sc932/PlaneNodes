c==============================
      program main
      implicit none
c====== For Rands ==========
      integer seed
      double precision rands(2)
c====== for timing ==========
      integer starttime, endtime
c====== Parameters ==========
      integer boxsize
      integer waves
      integer res
      integer iter
      double precision wavelength
      double precision dres
      parameter (boxsize=256)
      parameter (waves=1000)
      parameter (res=2)
      parameter (dres=0.5) ! 1/res
      parameter (iter=10)
      parameter (wavelength = 1.0)
c====== Variables for code ========
      integer i,j,k, xi, yi
      double precision psir(boxsize*res,boxsize*res)
      double precision psi(boxsize*res,boxsize*res)
      double precision a(waves)
      double precision theta
      double precision delta
c====== Variables for data ========
      double precision fmean
      integer numMaxima(iter)
      integer tnumMaxima
      double precision maxima(iter,boxsize*boxsize/4)
      double precision tmaxima(boxsize*boxsize/4)
      integer maximaOrder(iter,boxsize*boxsize/4)
      integer tmaximaOrder(boxsize*boxsize/4)
      double precision extrema(iter)
      double precision meanExt
      double precision meanMaxima(iter)
      double precision meanNodes
      double precision meanestMaxima
      double precision std
      double precision meanEx
c====== End of declarations =======

c====== Output initialization =====
      open(unit = 2, file = 'pnR_L256_w1000_res2_iter10.dat')
      open(unit = 3, file = 'pnR_L256_w1000_res2_iter10_maxima.dat')
      open(unit = 4, file = 'pnR_L256_w1000_iter10_maximaFULL.dat')

c====== initialize variables ======
      fmean = 0.0
      meanExt = 0.0
      meanNodes = 0.0
      meanestMaxima = 0.0
      std = 0.0
      meanEx = 0.0
      do i = 1, (waves)
         a(i) = 0D0
      enddo
      do i = 1, (iter)
         numMaxima(i) = 0
         extrema(i) = 0D0
         meanMaxima(i) = 0D0
         do j = 1, (boxsize*boxsize/4)
            maxima(i,j) = 0D0
            maximaOrder(i,j) = boxsize*boxsize/4
         enddo
      enddo
c======= end initialization =======
      
      write(6,*) 'Program Starting'
      starttime = second()
      call zufalli(starttime)
c======= Start program ============
      do j = 1, iter

      if (j.eq.(iter/2)) then
            write(6,*) 'Half way done'
            write(6,*) (second() - starttime)
      endif

c========= Initialize the arrays ==
         do i = 1, (boxsize*res)
            do k = 1, (boxsize*res)
               psir(i,k) = 0D0
            enddo
         enddo
         
         call normalen(waves,a)
c========= Generate the waves ===
         do i = 1, (waves)
            call zufall(2,rands)
            theta = rands(1)*2*3.14159256
            delta = rands(2)*2*3.14159256
            do xi = 1, (boxsize*res)
               do yi = 1, (boxsize*res)
                  psir(xi,yi) = psir(xi,yi) +
     $            a(i)*cos((xi)*dres*cos(theta)/wavelength +
     $            (yi)*dres*sin(theta)/wavelength + delta)
               enddo
            enddo
         enddo

c========= Combine waves by complex conjugate =====
         do xi = 1, (boxsize*res)
            do yi = 1, (boxsize*res)
               psi(xi,yi) = psir(xi,yi)*psir(xi,yi)
               fmean = fmean + psi(xi,yi)
            enddo
         enddo

c========== Start local maximum finder ============
         do xi = 2, (boxsize*res-1)
            do yi = 2, (boxsize*res-1)
               if (psi(xi-1,yi-1).lt.psi(xi,yi)) then
               if (psi(xi-1,yi).lt.psi(xi,yi)) then
               if (psi(xi-1,yi+1).lt.psi(xi,yi)) then
               if (psi(xi,yi+1).lt.psi(xi,yi)) then
               if (psi(xi+1,yi+1).lt.psi(xi,yi)) then
               if (psi(xi+1,yi).lt.psi(xi,yi)) then
               if (psi(xi+1,yi-1).lt.psi(xi,yi)) then
               if (psi(xi,yi-1).lt.psi(xi,yi)) then
                  numMaxima(j) = numMaxima(j) + 1
                  maxima(j,numMaxima(j)) = psi(xi,yi)
               endif
               endif
               endif
               endif
               endif
               endif
               endif
               endif
            enddo
         enddo

c======= maxima sorter =======================================
         do i = 1, (numMaxima(j))
            tmaxima(i) = maxima(j,i)
         enddo
         call QSORTI(tmaximaOrder,numMaxima(j),tmaxima)
         do i = 1, (numMaxima(j))
            maximaOrder(j,i) = tmaximaOrder(i)
         enddo

c========== Find mean of maxima and extrema =============
         extrema(j) = maxima(j,1)
         do i = 1,(numMaxima(j))
            meanMaxima(j) = meanMaxima(j) + maxima(j,i)
            if (maxima(j,i).gt.extrema(j)) then
               extrema(j) = maxima(j,i)
            endif
            if (j.lt.10) then
               write(4,*) maxima(j,i)
            endif
         enddo
         if (numMaxima(j).gt.0) then
            meanMaxima(j) = meanMaxima(j)/numMaxima(j)
         endif
      enddo
      
c====== Standard dev of extreama finder ====================
      do i = 1, iter
         meanEx = meanEx + extrema(i)
      enddo
      meanEx = meanEx/iter
      do i = 1, iter
         std = std + (extrema(i) - meanEx)*(extrema(i) - meanEx)
      enddo
      std = sqrt(std/iter)

      do i = 1, iter
         meanestMaxima = meanestMaxima + meanMaxima(i)
         meanNodes = meanNodes + numMaxima(i)
         meanExt = meanExt + extrema(i)
      enddo
      meanestMaxima = meanestMaxima/iter
      meanNodes = meanNodes/iter
      meanExt = meanExt/iter

      endtime = second()
      write(6,*) 'Computation done, ',iter,' iterations'
c========== Output results ===========================
      write(6,*) 'Area: ', boxsize
      write(6,*) 'Waves: ', waves
      write(6,*) 'iter: ', iter
      write(6,*) 'time: ', (endtime-starttime)
      write(6,*) 'mean extrema: ', meanExt
      write(6,*) 'mean maxima: ', meanestMaxima
      write(6,*) 'mean field: ', (fmean/(1000*boxsize*boxsize*res*res))
      write(6,*) 'mean nodes: ', meanNodes
c========== to file ====================
      write(6,*) 'Writing output...'
      write(2,*) 'Area: ', boxsize
      write(2,*) 'Waves: ', waves
      write(2,*) 'iter: ', iter
      write(2,*) 'time: ', (endtime-starttime)
      write(2,*) 'mean extrema: ', meanExt
      write(2,*) 'mean maxima: ', meanestMaxima
      write(2,*) 'mean field: ', (fmean/(1000*boxsize*boxsize*res*res))
      write(2,*) 'mean nodes: ', meanNodes
      do j = 1, (iter)
         write(2,*) 'iteration: ', j
         write(2,*) 'extrema: ', extrema(j)
         write(2,*) 'meanMaxima: ', meanMaxima(j)
         write(2,*) 'numMaxima: ', numMaxima(j)
         if(boxsize.gt.16) then
            write(3,*) 'START ', j
            do i = 1, 20
               write(3,*) maxima(j,maximaOrder(j,i))
            enddo
            write(3,*) 'MID ', j
            do i = (numMaxima(j)-20),(numMaxima(j))
               write(3,*) maxima(j,maximaOrder(j,i))
            enddo
            write(3,*) 'END ', j
         endif
      enddo

      close(2)
      close(3)

      end

* README for zufall random number package
* ------ --- ------ ------ ------ -------
* This package contains a portable random number generator set
* for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
* Poisson distributions. The basic module, the uniform generator,
* uses a lagged Fibonacci series generator:
* 
*               t    = u(n-273) + u(n-607)
*               u(n) = t - float(int(t))
* 
* where each number generated, u(k), is floating point. Since
* the numbers are floating point, the left end boundary of the
* range contains zero. This package is nearly portable except
* for the following. (1) It is written in lower case, (2) the
* test package contains a timer (second) which is not portable,
* and (3) there are cycle times (in seconds) in data statements 
* for NEC SX-3, Fujitsu VP2200, and Cray Y-MP. Select your 
* favorite and comment out the others. Replacement functions 
* for 'second' are included - comment out the others. Otherwise 
* the package is portable and returns the same set of floating 
* point numbers up to word precision on any machine. There are 
* compiler directives ($cdir for Cray, *vdir for SX-3, and VOCL 
* for Fujitsu VP2200) which should be otherwise ignored.
* 
* To compile this beast, note that all floating point numbers
* are declared 'double precision'. On Cray X-MP, Y-MP, and C-90
* machines, use the cft77 (cf77) option -dp to run this in 64
* bit mode (not 128 bit double).
* 
* External documentation, "Lagged Fibonacci Random Number Generators
* for the NEC SX-3," is to be published in the International
* Journal of High Speed Computing (1994). Otherwise, ask the
* author: 
* 
*          W. P. Petersen 
*          IPS, RZ F-5
*          ETHZ
*          CH 8092, Zurich
*          Switzerland
* 
* e-mail:  wpp@ips.ethz.ch.
* 
* The package contains the following routines:
* 
* ------------------------------------------------------
* UNIFORM generator routines:
* 
*       subroutine zufalli(seed)
*       integer seed
* c initializes common block containing seeds. if seed=0,
* c the default value is 1802.
* 
*       subroutine zufall(n,u)
*       integer n
*       double precision u(n)
* c returns set of n uniforms u(1), ..., u(n).
* 
*       subroutine zufallsv(zusave)
*       double precision zusave(608)
* c saves buffer and pointer in zusave, for later restarts
* 
*       subroutine zufallrs(zusave)
*       double precision zusave(608)
* c restores seed buffer and pointer from zusave
* ------------------------------------------------------
* 
* NORMAL generator routines:
* 
*       subroutine normalen(n,g)
*       integer n
*       double precision g(n)
* c returns set of n normals g(1), ..., g(n) such that
* c mean <g> = 0, and variance <g**2> = 1.
* 
*       subroutine normalsv(normsv)
*       double precision normsv(1634)
* c saves zufall seed buffer and pointer in normsv
* c buffer/pointer for normalen restart also in normsv
* 
*       subroutine normalrs(normsv)
*       double precision normsv(1634)
* c restores zufall seed buffer/pointer and 
* c buffer/pointer for normalen restart from normsv
* ------------------------------------------------------
* 
* POISSON generator routine:
* 
*       subroutine fische(n,mu,q)
*       integer n,q(n)
*       double precision mu
* c returns set of n integers q, with poisson
* c distribution, density p(q,mu) = exp(-mu) mu**q/q!
* c 
* c USE zufallsv and zufallrs for stop/restart sequence
* c

      subroutine zufallt(n,a)
      implicit none
      integer n
c
      double precision a(n)
      integer ia(20)
      double precision diff,t0,t1,t2,t3,second
      double precision svblk(608)
      double precision b(607),buff(607)
      double precision CYCLE
      integer i,ii,k,nits,ptr
      common /klotz0/buff,ptr
c
c clock cycle for machine
c
c cycle for SX-3:
c     data CYCLE/2.9E-9/
c cycle for Y-MP:
c     data CYCLE/6.0E-9/
c cycle for VP2200
c     data CYCLE/3.2E-9/
c cycle for SGI Indigo
c     data CYCLE/1.0E-8/
c cycle for Sparc 51
      data CYCLE/2.0E-8/
c
c number of iterations of test: nits
c
      nits = 128
c
      do 1 k=1,20
         ia(k) = 0
1     continue
      t0 = 100.
      do 2 k=1,nits
         t1 = second()
         t2 = second()
         t1 = t2 - t1
         t0 = min(t0,t1)
2     continue
      t1 = 100.
      do 3 k=1,nits
         t2 = second()
         call zufall(n,a)
         t3 = second()
         t2 = t3 - t2
         t1 = min(t2,t1)
         do 4 i=1,n
            ii     = int(a(i)*20.)+1
            ia(ii) = ia(ii) + 1
4        continue
c
c  last time, save klotz0 for save/resore test
c
         if(k.eq.nits-1)then
            call zufallsv(svblk)
         endif
c
3     continue
c
c  test save/restore sequence
c
      call zufallrs(svblk)
      call zufall(607,b)
      diff = 0.
      do 5 i=1,min(n,607)
         diff = diff + abs(b(i) - a(i))
5     continue
      if(diff.ne.0.) then
         print *,' ERROR in start/restart: diff = ',diff
      else
         print *,' zufall save/restore test OK'
      endif
c
      t1 = (t1 - t0)/float(n)
      print 100,t1
      print 200,t1/CYCLE
      print 300,(k,ia(k),k=1,20)
100   format(/1x,' Time/uniform = ',e12.3,' seconds')
200   format(1x,' Cps./uniform = ',e12.3)
300   format(/1x,' Uniform Histogram:',/,
     *        1x,' ------- ---------',/,
     *        20(/1x,' bin(',i2,') = ',i9))
      return
      end
c
      subroutine normalt(n,x)
      implicit none
      integer n
      double precision x(n),y(128),boxsv(1634)
      double precision x1,x2,x3,x4,x5,x6,xx2,xx4
      double precision diff,t0,t1,t2,t3,second
      integer bin(21),i,k,kk,nits
c
      double precision CYCLE
c
c clock cycle for machine
c
c cycle for SX-3:
c     data CYCLE/2.9E-9/
c cycle for Y-MP:
c     data CYCLE/6.0E-9/
c cycle for VP2200
c     data CYCLE/3.2E-9/
c cycle for SGI Indigo
c     data CYCLE/1.0E-8/
c cycle for Sparc 51
      data CYCLE/2.0E-8/
c
c number of iterations of test
c
      nits = 128
c
c initialize moments
c
      x1 = 0.
      x2 = 0.
      x3 = 0.
      x4 = 0.
      x5 = 0.
      x6 = 0.
c
      call normalen(n,x)
      do 1 i = 1,21
         bin(i) = 0
1     continue
c
      t0 = 10.
      do 2 k = 1,nits
         t1 = second()
         t2 = second()
         t1 = t2 - t1
         t0 = min(t1,t0)
2     continue
      t1 = 100.
      do 3 k = 1,nits
c
c  save seeds and pointers for save/restore test
c
         if(k.eq.nits) call normalsv(boxsv)
c
         t2 = second()
         call normalen(n,x)
         t3 = second()
         t2 = t3 - t2
         t1 = min(t1,t2)
c
         do 4 i=1,n
            kk = int(2.0*(x(i)+5.25)) + 1
            bin(kk) = bin(kk) + 1
4        continue
         do 5 i=1,n
            x1  = x1 + x(i)
            xx2 = x(i)*x(i)
            x2  = x2 + xx2
            x3  = x3 + xx2*x(i)
            xx4 = xx2*xx2
            x4  = x4 + xx4
            x5  = x5 + xx4*x(i)
            x6  = x6 + xx4*xx2
5        continue
c
c  restore previous seeds and pointers for save/restore test
c
         if(k.eq.nits) then
            call normalrs(boxsv)
            call normalen(128,y)
         endif
c
3     continue
c
c  save/restore check:
c
      do 6 i=1,128
         diff = diff + abs(y(i) - x(i))
6     continue
      if(diff.ne.0.) then
         print *,' ERROR in normalsv/normalrs: diff = ',diff
      else
         print *,' normalen save/restore test OK'
      endif
c
      x1 = x1/float(n*nits)
      x2 = x2/float(n*nits)
      x3 = x3/float(n*nits)
      x4 = x4/float(n*nits)
      x5 = x5/float(n*nits)
      x6 = x6/float(n*nits)
c
      t1 = (t1 - t0)/float(n)
      print 100,t1
      print 200,t1/CYCLE
      print 300,x1,x2,x3,x4,x5,x6
      print 400
      do 7 k=1,21
         print 500,k,bin(k)
7     continue
c
100   format(/1x,' Time/normal = ',e12.3,' seconds')
200   format(1x,' Cps./normal = ',e12.3)
300   format(/1x,' Moments:'/,
     *   4x,'Compare to:  (0.0)',18x,'(1.0)'/
     *   4x,' <x>    = ',e12.5,', <x**2> = ',e12.5,//
     *   4x,'Compare to:  (0.0)',18x,'(3.0)'/
     *   4x,' <x**3> = ',e12.5,', <x**4> = ',e12.5,//
     *   4x,'Compare to:  (0.0)',17x,'(15.0)'/
     *   4x,' <x**5> = ',e12.5,', <x**6> = ',e12.5)

400   format(/1x,' Histogram of gaussian distribution'/,
     *      1x,' --------- -- -------- ------------'/)
500   format(1x,' bin(',i2,') = ',i7)
c
      return
      end
c
      subroutine fischet(n,p)
      implicit none
      double precision mu,fp
      double precision p1,p2,p3,p4
      double precision x1,x2,x3,x4
      double precision t0,t1,t2,t3,second
      integer bin(20)
      integer i,k,kk,n,nits
      integer p(n)
c
      double precision CYCLE
c
c clock cycle for machine
c
c cycle for SX-3:
c     data CYCLE/2.9E-9/
c cycle for Y-MP:
c     data CYCLE/6.0E-9/
c cycle for VP2200
c     data CYCLE/3.2E-9/
c cycle for SGI Indigo
c     data CYCLE/1.0E-8/
c cycle for Sparc 51
      data CYCLE/2.0E-8/
c
      mu   = 2.0
      nits = 128
c
      do 1 k=1,20
         bin(k) = 0
1     continue
c
c moment comparison values
c
      p1 = mu
      p2 = mu + mu*mu
      p3 = mu + 3.*mu*mu + mu*mu*mu
      p4 = mu + 7.*mu*mu + 6.*mu*mu*mu + mu**4
c
      x1 = 0.
      x2 = 0.
      x3 = 0.
      x4 = 0.
c
      t0 = 10.
      do 2 k=1,nits
         t1 = second()
         t2 = second()
         t1 = t2 - t1
         t0 = min(t0,t1)
2     continue
      t1 = 10.
      do 3 k=1,nits
c
         t2 = second()
         call fische(n,mu,p)
         t3 = second()
         t2 = t3 - t2
         t1 = min(t1,t2)
c
         do 4 i=1,n
            kk = p(i)+1
            bin(kk) = bin(kk) + 1
4        continue
c
         do 5 i=1,n
            fp = float(p(i))
            x1 = x1 + fp
            x2 = x2 + fp*fp
            x3 = x3 + fp*fp*fp
            x4 = x4 + fp*fp*fp*fp
5        continue
c
3     continue
c
      x1 = x1/float(n*nits)
      x2 = x2/float(n*nits)
      x3 = x3/float(n*nits)
      x4 = x4/float(n*nits)
c
      t1 = (t1 - t0)/float(n)
      print 100,t1
      print 200,t1/CYCLE
      print 300,p1,p2,x1,x2,p3,p4,x3,x4
      print 400,mu
      do 6 k=1,20
         print 500,k,bin(k)
6     continue
c
100   format(1x,' Time/poisson = ',e12.3,' seconds ')
200   format(1x,' Cps./poisson = ',e12.3)
300   format(/1x,' Moments:'/,
     *   3x,'Compare: (',e12.5,')         (',e12.5,')'/
     *   3x,' <p>    = ',e12.5,', <p**2> = ',e12.5,//
     *   3x,'Compare: (',e12.5,')         (',e12.5,')'/
     *   3x,' <p**3> = ',e12.5,', <p**4> = ',e12.5/)
400   format(/1x,' Histogram of Poisson distribution: mu = ',f8.3/,
     *       1x,' --------- -- ------- ------------'/)
500   format(1x,' bin(',i2,') = ',i7)
      return
      end
c
c      double precision function second()
c
c  portable version of Cray function second(), comment out
c  entire function for Y-MP. For SX-3, Fujitsu VP2200,
c  or generic Unix, uncomment the version you want:
c
c NEC SX-3 version
c     double precision xx(2)
c     call clock(xx)
c
c VP2200 version
c     double precision xx(2)
c     call clockv(xx(2),xx(1),0,2)
c
c Generic Unix version
c      real etime
c      real time(2)
c      call etime(time)
c
c      second = time(1)
c
c      return
c      end
c
c ---------------- end of test programs -------------
c
      subroutine zufall(n,a)
      implicit none
c
c portable lagged Fibonacci series uniform random number
c generator with "lags" -273 und -607:
c
c       t    = u(i-273)+buff(i-607)  (floating pt.)
c       u(i) = t - float(int(t))
c
c W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
c
      double precision a(*)
      double precision buff(607)
      double precision t
      integer i,k,ptr,VL,k273,k607
      integer buffsz,nn,n,left,q,qq
      integer aptr,aptr0,bptr
c
      common /klotz0/buff,ptr
      data buffsz/607/
c
      aptr = 0
      nn   = n
c
1     continue
c
      if(nn .le. 0) return
c
c factor nn = q*607 + r
c
      q    = (nn-1)/607
      left = buffsz - ptr
c
      if(q .le. 1) then
c
c only one or fewer full segments
c
         if(nn .lt. left) then
            do 2 i=1,nn
               a(i+aptr) = buff(ptr+i)
2           continue
            ptr  = ptr + nn
            return
         else
            do 3 i=1,left
               a(i+aptr) = buff(ptr+i)
3           continue
            ptr  = 0
            aptr = aptr + left
            nn   = nn - left
c  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue
c
            goto 1
         endif
      else
c
c more than 1 full segment
c 
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
c 
c buff -> a(aptr0)
c 
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
*VOCL LOOP, TEMP(t)
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t)
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607
c
c a(aptr-607) -> a(aptr) for last of the q-1 segments
c
          aptr0 = aptr - 607
          VL    = 607
c
*vdir novector
          do 10 qq=1,q-2
             k273 = 334 + aptr0
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(a)
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue
c
c a(aptr0) -> buff, last segment before residual
c
          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
*VOCL LOOP, TEMP(t)
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      end
c
      subroutine zufalli(seed)
      implicit none
c
c  generates initial seed buffer by linear congruential
c  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
c  variable seed should be 0 < seed <31328
c
      integer seed
      integer ptr
      double precision s,t
      double precision buff(607)
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/
c
      if(seed.ne.0) ij = seed
c
      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do 1 ii=1,607
         s = 0.0
         t = 0.5
         do 2 jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s = s+t
            t = .5*t
2        continue
         buff(ii) = s
1     continue
      return
      end
c
      subroutine zufallsv(svblk)
      implicit none
c
c  saves common blocks klotz0, containing seeds and 
c  pointer to position in seed block. IMPORTANT: svblk must be
c  dimensioned at least 608 in driver. The entire contents
c  of klotz0 (pointer in buff, and buff) must be saved.
c
      double precision buff(607)
      integer ptr,i
      double precision svblk(*)
      common /klotz0/buff,ptr
c
      svblk(1) = ptr
      do 1 i=1,607
         svblk(i+1) = buff(i)
1     continue
c
      return
      end
      subroutine zufallrs(svblk)
      implicit none
c
c  restores common block klotz0, containing seeds and pointer
c  to position in seed block. IMPORTANT: svblk must be
c  dimensioned at least 608 in driver. The entire contents
c  of klotz0 must be restored.
c
      double precision buff(607)
      integer i,ptr
      double precision svblk(*)
      common /klotz0/buff,ptr
c
      ptr = svblk(1)
      do 1 i=1,607
         buff(i) = svblk(i+1)
1     continue
c
      return
      end
      subroutine normalen(n,x)
      implicit none
c
c Box-Muller method for Gaussian random numbers
c
      double precision x(*)
      double precision xbuff(1024)
      integer i,ptr,xptr,first
      integer buffsz,nn,n,left 
      common /klotz1/xbuff,first,xptr
      data buffsz/1024/
c
      nn   = n
      if(nn .le. 0) return
      if(first.eq.0)then
         call normal00
         first = 1
      endif
      ptr = 0
c
1     continue
      left = buffsz - xptr
      if(nn .lt. left) then
         do 2 i=1,nn
            x(i+ptr) = xbuff(xptr+i)
2        continue
         xptr = xptr + nn
         return
      else
         do 3 i=1,left
            x(i+ptr) = xbuff(xptr+i)
3        continue
         xptr = 0
         ptr  = ptr+left
         nn   = nn - left
         call normal00
         goto 1
      endif
      end
      subroutine normal00
      implicit none
      double precision pi,twopi
      parameter(pi=3.141592653589793)
      double precision xbuff(1024),r1,r2,t1,t2
      integer first,xptr,i
      common /klotz1/xbuff,first,xptr
c
      twopi = 2.*pi
      call zufall(1024,xbuff)
*VOCL LOOP, TEMP(r1,r2,t1,t2), NOVREC(xbuff)
      do 1 i=1,1024,2
         r1         = twopi*xbuff(i)
         t1         = cos(r1)
         t2         = sin(r1)
         r2         = sqrt(-2.*log(1.-xbuff(i+1)))
         xbuff(i)   = t1*r2
         xbuff(i+1) = t2*r2
1     continue
c
      return
      end
      subroutine normalsv(svbox)
      implicit none
c
c  saves common block klotz0 containing buffers
c  and pointers. IMPORTANT: svbox must be dimensioned at 
c  least 1634 in driver. The entire contents of blocks 
c  klotz0 (via zufallsv) and klotz1 must be saved.
c
      double precision buff(607)
      integer i,k,ptr
      double precision xbuff(1024)
      integer xptr,first
      double precision svbox(*)
      common /klotz0/buff,ptr
      common /klotz1/xbuff,first,xptr
c
      if(first.eq.0)then
         print *,' ERROR in normalsv, save of unitialized block'
      endif
c
c  save zufall block klotz0
c
      call zufallsv(svbox)
c
      svbox(609) = first
      svbox(610) = xptr
      k = 610
      do 1 i=1,1024
         svbox(i+k) = xbuff(i)
1     continue
c
      return
      end
      subroutine normalrs(svbox)
      implicit none
c
c  restores common blocks klotz0, klotz1 containing buffers
c  and pointers. IMPORTANT: svbox must be dimensioned at 
c  least 1634 in driver. The entire contents
c  of klotz0 and klotz1 must be restored.
c
      double precision buff(607)
      integer ptr
      double precision xbuff(1024)
      integer i,k,xptr,first
      double precision svbox(*)
      common /klotz0/buff,ptr
      common /klotz1/xbuff,first,xptr
c
c restore zufall blocks klotz0 and klotz1
c
      call zufallrs(svbox)
      first = svbox(609)
      if(first.eq.0)then
         print *,' ERROR in normalsv, restoration of unitialized block'
      endif
      xptr  = svbox(610)
      k = 610
      do 1 i=1,1024
         xbuff(i) = svbox(i+k)
1     continue
c
      return
      end
      subroutine fische(n,mu,p)
      implicit none
      integer p(*)
      integer indx(1024)
      integer n,i,ii,jj,k,left,nl0,nsegs,p0
      double precision u(1024),q(1024)
      double precision q0,pmu,mu
c
c Poisson generator for distribution function of p's:
c
c    q(mu,p) = exp(-mu) mu**p/p!
c
c initialize arrays, pointers
c
      if (n.le.0) return
c
      pmu = exp(-mu)
      p0  = 0
c
      nsegs = (n-1)/1024 
      left  = n - nsegs*1024
      nsegs = nsegs + 1
      nl0   = left
c
      do 2 k = 1,nsegs
c
         do 3 i=1,left
            indx(i)    = i
            p(p0+i)    = 0
            q(i)       = 1.0
3        continue
c
c Begin iterative loop on segment of p's
c
1        continue
c
c Get the needed uniforms
c
         call zufall(left,u)
c
         jj = 0
c
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(ii,q0), NOVREC(indx,p,q)
         do 4 i=1,left
            ii    = indx(i)
            q0    = q(ii)*u(i)
            q(ii) = q0
            if( q0.gt.pmu ) then
               jj       = jj + 1
               indx(jj) = ii
               p(p0+ii) = p(p0+ii) + 1
            endif
4        continue
c
c any left in this segment?
c
         left = jj
         if(left.gt.0)then
            goto 1
         endif
c
         p0    = p0 + nl0
         nl0   = 1024
         left  = 1024
c
2     continue
c
      return
      end
c
      block data
      implicit none
c
c globally accessable, compile-time initialized data
c
      integer ptr,xptr,first
      double precision buff(607),xbuff(1024)
      common /klotz0/buff,ptr
      common /klotz1/xbuff,first,xptr
      data ptr/0/,xptr/0/,first/0/
      end

C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.

!http://www.fortran.com/quick_sort1.f
      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      INTEGER X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      double precision A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END

