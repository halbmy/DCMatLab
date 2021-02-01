      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C--------------------------------------------------------------
C     sens2dxz(x,z,pa,pm)
C
      integer plhs(*), prhs(*)
      integer mxCreateDoubleMatrix, mxGetPr
      integer xe_p, ze_p, pa_p, pm_p, out_p
      integer*4 nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      integer m, n, lx, lz, size
      integer i, j, k, nx, nz
      real*8 xe(2000), ze(100), pa(2), pm(2)
      real*8 out(20000), iout, jout, kout
      real*8 dxm, dzm, dxa, dza, ee, ek, t1, t2
      real*8 x, z, tmp, aq, bq, a, dx, dz
      real*8 xg(5), wg(5)
      real*8 yint
      data xg / 0.04691007703067d0,0.23076534494716d0,0.5d0,
     &          0.76923465505284d0,0.95308992296933d0 /
      data wg / 0.11846340022799d0,0.23931433524315d0,
     &          0.28444444444444d0,
     &          0.23931433524315d0,0.11846340022799d0 / 
C     Check for proper number and size of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgTxt('four inputs arguments required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output argument required.')
      endif
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      lx=m*n-1
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      lz=m*n-1
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      size = m*n
      if(size.ne.2) then
         call mexErrMsgTxt('pa has to be a 2 vector')
      endif
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      size = m*n
      if(size.ne.2) then
         call mexErrMsgTxt('pm has to be a 2 vector')
      endif
      plhs(1) = mxCreateDoubleMatrix(lx*lz, 1, 0)
      xe_p = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(xe_p, xe, lx+1)
      ze_p = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(ze_p, ze, lz+1)
      pa_p = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(pa_p, pa, 2)
      pm_p = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(pm_p, pm, 2)
      out_p = mxGetPr(plhs(1))
      n=5
      do 50 nx=1,lx
        dx=xe(nx+1)-xe(nx)
        do 40 nz=1,lz
          dz=ze(nz+1)-ze(nz)
          iout=0.0
          do 30 i=1,n
            x=xe(nx)+dx*xg(i)
            jout=0.0
            do 20 j=1,n
              z=ze(nz)+dz*xg(j)
              m=1
              tmp=yint(x,z,pa(1),pa(2),pm(1),pm(2))
C possible cases for subsurface electrodes              
              if(pa(2).ne.0) then
                tmp=tmp+yint(x,z,pa(1),-pa(2),pm(1),pm(2))
                m=m+1
                if(pm(2).ne.0) then
                  tmp=tmp+yint(x,z,pa(1),-pa(2),pm(1),-pm(2))
                  m=m+1
                endif
              endif
              if(pm(2).ne.0) then
                tmp=tmp+yint(x,z,pa(1),pa(2),pm(1),-pm(2))
                m=m+1
              endif
              jout=jout+tmp*wg(j)/m
20          continue
            iout=iout+jout*wg(i)
30        continue
          out(nx+(nz-1)*lx)=iout*dx*dz/39.478417604
40      continue
50    continue

      call mxCopyReal8ToPtr(out, out_p, lx*lz)     
      return
      end
      
      real*8 function yint(x,z,xa,za,xm,zm)
      real*8 x,z,xa,za,xm,zm,t1,t2,a
      real*8 dxa, dza, dxm, dzm, aq, bq, ee, ek, tmp
      dxa=x-xa
      dxm=x-xm
      dza=z-za
      dzm=z-zm
      aq=dxa**2+dza**2
      bq=dxm**2+dzm**2
      if(bq.gt.aq) then
        tmp=aq
        aq=bq
        bq=tmp
      endif
      a=sqrt(aq)
      if((aq-bq).lt.1e-10) then
        tmp=aq*a;
        t1=1.17809724509617/tmp/aq;
        t2=0.39269908169872/tmp;
      else
        call ellipke(1-bq/aq,ek,ee)
        t2=2/(a*(aq-bq)**2)
        t1=t2*((aq+bq)*ee-2*bq*ek)/bq
        t2=t2*((aq+bq)*ek-2*aq*ee)
      endif
      yint=(t1*(dxa*dxm+dza*dzm)+t2)      
      return
      end
      
      subroutine ellipke(m,ek,ee)
      real*8 m, ek, ee
      real*8 a0, b0, s0, a1, b1, c1, w1
      integer i1
      a0=1
      b0=sqrt(1-m)
      s0=m
      i1=0
100   a1=(a0+b0)/2
      b1=sqrt(a0*b0)
      c1=(a0-b0)/2
      i1=i1 + 1
      w1=2**i1*c1**2
      s0=s0+w1
      a0=a1
      b0=b1
      if(w1.gt.1e-15) then
        goto 100
      endif
      ek=1.57079632679490/a1
      ee=ek*(1-s0/2)
      if(m.eq.0) then
        ee=1
        ek=1e20
      endif
      return
      end
