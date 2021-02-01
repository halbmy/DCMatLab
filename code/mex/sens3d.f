      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C--------------------------------------------------------------
C     sens3d(p1,p2,pa,pm,n)
C
      integer plhs(*), prhs(*)
      integer mxCreateFull, mxGetPr
      integer p1_p, p2_p, pa_p, pm_p, out_p
      integer*4 nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      integer m, n, size, i, j, k, ns, ne
      real*8 p1(3), p2(3), pa(3), pm(3)
      real*8 out, jout, kout, dxm, dym, dzm, dxa, dya, dza
      real*8 x, y, z, tmp, ra, rm, dx, dy, dz
      real*8 xg(21), wg(21)
      data xg / 0.5,0.21132486540519,0.78867513459481,
     &          0.11270166537926,0.5,0.88729833462074,
     &          0.06943184420297,0.33000947820712,0.66999052179288,
     &          0.93056815579703,0.04691007703067,0.23076534494716,0.5,
     &          0.76923465505284,0.95308992296933,0.03376524289842,
     &          0.16939530676687,0.38069040695840,0.61930959304160,
     &          0.83060469323313,0.96623475710158/
      data wg / 1.0,0.5,0.5,0.27777773678784,0.44444444444444,
     &          0.27777773678784,0.17392737578631,0.32607180445171,
     &          0.32607180445171,0.17392737578631,0.11846340022799,
     &          0.23931433524315,0.28444444444444,0.23931433524315,
     &          0.11846340022799,0.08566221039227,0.18038078650699,
     &          0.23395693935088,0.23395693935088,0.18038078650699,
     &          0.08566221039227/
C     Check for proper number and size of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgTxt('Four inputs required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('p1 has to be a 3 vector')
      endif
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('p2 has to be a 3 vector')
      endif
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('pa has to be a 3 vector')
      endif
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      size = m*n
      if(size.ne.3) then
         call mexErrMsgTxt('pm has to be a 3 vector')
      endif
      
      plhs(1) = mxCreateFull(1, 1, 0)
      p1_p = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(p1_p, p1, 3)
      p2_p = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(p2_p, p2, 3)
      pa_p = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(pa_p, pa, 3)
      pm_p = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(pm_p, pm, 3)
      out_p = mxGetPr(plhs(1))
      
      n=1
      ne=0
      do 5 i=1,n
        ne=ne+i
5     continue
      ns=ne-n+1
      out=0.0
      dx=p2(1)-p1(1)
      dy=p2(2)-p1(2)
      dz=p2(3)-p1(3)
      do 30 i=ns,ne
        x=p1(1)+dx*xg(i)
        jout=0.0
        do 20 j=ns,ne
          y=p1(2)+dy*xg(j)
          kout=0.0
          do 10 k=ns,ne
            z=p1(3)+dz*xg(k)
            dxa=x-pa(1)
            dxm=x-pm(1)
            dya=y-pa(2)
            dym=y-pm(2)
            dza=z-pa(3)
            dzm=z-pm(3)
            ra=dxa**2+dya**2+dza**2
            rm=dxm**2+dym**2+dzm**2
            ra=ra*sqrt(ra)
            rm=rm*sqrt(rm)
            ram=ra*rm
            ram=sqrt(ram)*ram
            tmp=dxa*dxm+dya*dym+dza*dzm
            kout=kout+tmp/ra/rm*wg(k)
10        continue
          jout=jout+kout*wg(j)
20      continue
        out=out+jout*wg(i)
30    continue
      out=out*dx*dy*dz/39.47841760435743
      call mxCopyReal8ToPtr(out, out_p, 1)     
      return
      end
