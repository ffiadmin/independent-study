/* CLAPACK 3.0 BLAS wrapper macros
 * Feb 5, 2000
 */
#include "gsl_cblas.h"
#ifndef __BLASWRAP_H
#define __BLASWRAP_H

#ifndef NO_BLAS_WRAP
 
/* BLAS1 routines */
#define srotg_ cblas_srotg
#define drotg_ cblas_drotg
#define srotmg_ cblas_srotmg
#define drotmg_ cblas_drotmg
#define srot_ cblas_srot
#define drot_ cblas_drot
#define srotm_ cblas_srotm
#define drotm_ cblas_drotm
#define sswap_ cblas_sswap
#define dswap_ cblas_dswap
#define cswap_ cblas_cswap
#define zswap_ cblas_zswap
#define sscal_ cblas_sscal
#define dscal_ cblas_dscal
#define cscal_ cblas_cscal
#define zscal_ cblas_zscal
#define csscal_ cblas_csscal
#define zdscal_ cblas_zdscal
#define scopy_ cblas_scopy
#define dcopy_ cblas_dcopy
#define ccopy_ cblas_ccopy
#define zcopy_ cblas_zcopy
#define saxpy_ cblas_saxpy
#define daxpy_ cblas_daxpy
#define caxpy_ cblas_caxpy
#define zaxpy_ cblas_zaxpy
#define sdot_ cblas_sdot
#define ddot_ cblas_ddot
#define cdotu_ cblas_cdotu
#define zdotu_ cblas_zdotu
#define cdotc_ cblas_cdotc
#define zdotc_ cblas_zdotc
#define snrm2_ cblas_snrm2
#define dnrm2_ cblas_dnrm2
#define scnrm2_ cblas_scnrm2
#define dznrm2_ cblas_dznrm2
#define sasum_ cblas_sasum
#define dasum_ cblas_dasum
#define scasum_ cblas_scasum
#define dzasum_ cblas_dzasum
#define isamax_ cblas_isamax
#define idamax_ cblas_idamax
#define icamax_ cblas_icamax
#define izamax_ cblas_izamax
 
/* BLAS2 routines */
#define sgemv_ cblas_sgemv
#define dgemv_ cblas_dgemv
#define cgemv_ cblas_cgemv
#define zgemv_ cblas_zgemv
#define sgbmv_ cblas_sgbmv
#define dgbmv_ cblas_dgbmv
#define cgbmv_ cblas_cgbmv
#define zgbmv_ cblas_zgbmv
#define chemv_ cblas_chemv
#define zhemv_ cblas_zhemv
#define chbmv_ cblas_chbmv
#define zhbmv_ cblas_zhbmv
#define chpmv_ cblas_chpmv
#define zhpmv_ cblas_zhpmv
#define ssymv_ cblas_ssymv
#define dsymv_ cblas_dsymv
#define ssbmv_ cblas_ssbmv
#define dsbmv_ cblas_dsbmv
#define sspmv_ cblas_sspmv
#define dspmv_ cblas_dspmv
#define strmv_ cblas_strmv
#define dtrmv_ cblas_dtrmv
#define ctrmv_ cblas_ctrmv
#define ztrmv_ cblas_ztrmv
#define stbmv_ cblas_stbmv
#define dtbmv_ cblas_dtbmv
#define ctbmv_ cblas_ctbmv
#define ztbmv_ cblas_ztbmv
#define stpmv_ cblas_stpmv
#define dtpmv_ cblas_dtpmv
#define ctpmv_ cblas_ctpmv
#define ztpmv_ cblas_ztpmv
#define strsv_ cblas_strsv
#define dtrsv_ cblas_dtrsv
#define ctrsv_ cblas_ctrsv
#define ztrsv_ cblas_ztrsv
#define stbsv_ cblas_stbsv
#define dtbsv_ cblas_dtbsv
#define ctbsv_ cblas_ctbsv
#define ztbsv_ cblas_ztbsv
#define stpsv_ cblas_stpsv
#define dtpsv_ cblas_dtpsv
#define ctpsv_ cblas_ctpsv
#define ztpsv_ cblas_ztpsv
#define sger_ cblas_sger
#define dger_ cblas_dger
#define cgeru_ cblas_cgeru
#define zgeru_ cblas_zgeru
#define cgerc_ cblas_cgerc
#define zgerc_ cblas_zgerc
#define cher_ cblas_cher
#define zher_ cblas_zher
#define chpr_ cblas_chpr
#define zhpr_ cblas_zhpr
#define cher2_ cblas_cher2
#define zher2_ cblas_zher2
#define chpr2_ cblas_chpr2
#define zhpr2_ cblas_zhpr2
#define ssyr_ cblas_ssyr
#define dsyr_ cblas_dsyr
#define sspr_ cblas_sspr
#define dspr_ cblas_dspr
#define ssyr2_ cblas_ssyr2
#define dsyr2_ cblas_dsyr2
#define sspr2_ cblas_sspr2
#define dspr2_ cblas_dspr2
 
/* BLAS3 routines */
#define sgemm_ cblas_sgemm
#define dgemm_ cblas_dgemm
#define cgemm_ cblas_cgemm
#define zgemm_ cblas_zgemm
#define ssymm_ cblas_ssymm
#define dsymm_ cblas_dsymm
#define csymm_ cblas_csymm
#define zsymm_ cblas_zsymm
#define chemm_ cblas_chemm
#define zhemm_ cblas_zhemm
#define ssyrk_ cblas_ssyrk
#define dsyrk_ cblas_dsyrk
#define csyrk_ cblas_csyrk
#define zsyrk_ cblas_zsyrk
#define cherk_ cblas_cherk
#define zherk_ cblas_zherk
#define ssyr2k_ cblas_ssyr2k
#define dsyr2k_ cblas_dsyr2k
#define csyr2k_ cblas_csyr2k
#define zsyr2k_ cblas_zsyr2k
#define cher2k_ cblas_cher2k
#define zher2k_ cblas_zher2k
#define strmm_ cblas_strmm
#define dtrmm_ cblas_dtrmm
#define ctrmm_ cblas_ctrmm
#define ztrmm_ cblas_ztrmm
#define strsm_ cblas_strsm
#define dtrsm_ cblas_dtrsm
#define ctrsm_ cblas_ctrsm
#define ztrsm_ cblas_ztrsm

#endif /* NO_BLAS_WRAP */

#endif /* __BLASWRAP_H */
