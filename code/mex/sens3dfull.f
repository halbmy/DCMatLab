      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C--------------------------------------------------------------
C     sens3dfull(x,y,z,pa,pm)
C
      integer plhs(*), prhs(*)
      integer mxCreateDoubleMatrix, mxGetPr
      integer xe_p, ye_p, ze_p, pa_p, pm_p, out_p
      integer*4 nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      integer m, n, lx, ly, lz, size
      integer i, j, k, nx, ny, nz, ns, ne
      real*8 xe(100), ye(100), ze(100), pa(3), pm(3)
      real*8 out(120000), iout, jout, kout
      real*8 dxm, dym, dzm, dxa, dya, dza, z1, z2
      real*8 x, y, z, tmp, ra, rm, ram, dx, dy, dz
      real*8 xg(20), wg(20)
C     Nodes and weights for Gauss-Legendre-Integration for 2,5 and 10
      data xg /0,0,0.21132486540519d0,0.78867513459481d0,0,
     &         0.04691007703067d0,0.23076534494716d0,0.5d0,
     &         0.76923465505284d0,0.95308992296933d0,0.01304673574141d0,
     &         0.06746831665551d0,0.16029521585017d0,0.28330230293537d0,
     &         0.42556283050918d0,0.57443716949082d0,0.71669769706463d0,
     &         0.83970478414983d0,0.93253168334449d0,0.98695326425859d0/
      data wg /0,0,0.5d0,0.5d0,0,0.11846340022799d0,0.23931433524315d0,
     &         0.28444444444444d0,0.23931433524315d0,0.11846340022799d0,
     &         0.03333565450117d0,0.07472567454220d0,0.10954279105033d0,
     &         0.13463332456072d0,0.14776211198671d0,0.14776211198671d0,
     &         0.13463332456072d0,0.10954279105033d0,0.07472567454220d0,
     &         0.03333565450117d0/
C     Check for proper number and size of arguments. 
      if(nrhs .ne. 5) then
         call mexErrMsgTxt('Five inputs arguments required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output argument required.')
      endif
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      lx=m*n-1
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      ly=m*n-1
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      lz=m*n-1
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('pa has to be a 3 vector')
      endif
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('pm has to be a 3 vector')
      endif
      plhs(1) = mxCreateDoubleMatrix(lx*ly*lz, 1, 0)
      xe_p = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(xe_p, xe, lx+1)
      ye_p = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(ye_p, ye, ly+1)
      ze_p = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(ze_p, ze, lz+1)
      pa_p = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(pa_p, pa, 3)
      pm_p = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(pm_p, pm, 3)
      out_p = mxGetPr(plhs(1))
      do 60 nz=1,lz
        dz=ze(nz+1)-ze(nz)
        n=5
        if(ze(nz).eq.0) then
          n=10
        endif
        if(ze(nz).gt.5) then
          n=2
        endif
        ns=n+1
        ne=2*n
        do 50 nx=1,lx
          dx=xe(nx+1)-xe(nx)
          do 40 ny=1,ly
            dy=ye(ny+1)-ye(ny)
            iout=0.0
            do 30 i=ns,ne
              x=xe(nx)+dx*xg(i)
              jout=0.0
              do 20 j=ns,ne
                y=ye(ny)+dy*xg(j)
                kout=0.0
                do 10 k=ns,ne
                  z=ze(nz)+dz*xg(k)
                  dxa=x-pa(1)
                  dxm=x-pm(1)
                  dya=y-pa(2)
                  dym=y-pm(2)
                  dza=z-pa(3)
                  dzm=z-pm(3)
                  ra=dxa**2+dya**2+dza**2
                  rm=dxm**2+dym**2+dzm**2
                  ram=sqrt(ra*rm)*ra*rm
                  tmp=dxa*dxm+dya*dym+dza*dzm
                  kout=kout+tmp/ram*wg(k)
10              continue
                jout=jout+kout*wg(j)
20            continue
              iout=iout+jout*wg(i)
30          continue
            out(nx+(ny-1)*lx+(nz-1)*lx*ly)=iout*dx*dy*dz/39.47841760435743
40        continue
50      continue
60    continue

      call mxCopyReal8ToPtr(out, out_p, lx*ly*lz)     
      return
      end
