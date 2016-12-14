      subroutine setnthreadf(nthread, nthreadm)
c
      include 'NScluster_f.h'
c
      integer nthread, nthreadm, ncpu
      integer omp_get_num_procs
c
      ncpu = omp_get_num_procs()
      nthreadm = nthread
      if( nthreadm.eq.0 .or. nthreadm.gt.ncpu ) nthreadm = ncpu
      call OMP_SET_NUM_THREADS(nthreadm)
      return
      end
