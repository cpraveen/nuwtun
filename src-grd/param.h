      integer nparammax, nsurfmax, nsurf
      integer param_type
      parameter(nparammax=100)
      parameter(nsurfmax=100)
      integer ibeg (nsurfmax),
     1        jbeg (nsurfmax),
     2        iend (nsurfmax),
     3        jend (nsurfmax),
     4        blk  (nsurfmax),
     4        nspt (nsurfmax),
     5        nhhp (nsurfmax),
     6        nhhpo(nsurfmax)
      real    xw(nparammax,nsurfmax)
      real    xwb(nparammax,nsurfmax)

      common/iparams/nsurf, ibeg, jbeg, iend, jend, blk, nspt, nhhp,
     1               nhhpo, param_type
      common/rparams/xw, xwb

      integer nblkmax, nblks
      parameter(nblkmax=100)
      integer idim(nblkmax), jdim(nblkmax), ioffr(nblkmax)
      common/grids/nblks, idim, jdim, ioffr
