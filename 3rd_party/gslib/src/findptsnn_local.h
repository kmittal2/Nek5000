#ifndef FINDPTSNN_LOCAL_H
#define FINDPTSNN_LOCAL_H

#if !defined(MEM_H) || !defined(FINDPTSNN_EL_H) || !defined(OBBOX_H)
#warning "findptsnn_local.h" requires "mem.h", "findpnn_el.h", "obbox.h"
#endif

#define findptsnn_local_setup_2   PREFIXED_NAME(findptsnn_local_setup_2)
#define findptsnn_local_free_2    PREFIXED_NAME(findptsnn_local_free_2 )
#define findptsnn_local_2         PREFIXED_NAME(findptsnn_local_2      )
#define findptsnn_local_eval_2    PREFIXED_NAME(findptsnn_local_eval_2 )
#define findptsnn_local_el_eval_2    PREFIXED_NAME(findptsnn_local_el_eval_2 )

struct findptsnn_local_hash_data_2 {
  uint hash_n;
  struct dbl_range bnd[2];
  double fac[2];
  uint *offset;
  uint max;
};

struct findptsnn_local_data_2 {
  unsigned ntot;
  const double *elx[2];
  const unsigned *nsid;
  struct obbox_2 *obb;
  struct findptsnn_local_hash_data_2 hd;
  struct findptsnn_el_data_2 fed;
  double tol;
  double *distint;
  const double *distfint; //This is field for distance for all GLL points
};

void findptsnn_local_setup_2(struct findptsnn_local_data_2 *const fd,
                           const double *const elx[2],
                           const unsigned *const nsid,
                           const double *const distfint,
                           const unsigned n[2], const uint nel,
                           const unsigned m[2], const double bbox_tol,
                           const uint max_hash_size,
                           const unsigned npt_max, const double newt_tol);
void findptsnn_local_free_2(struct findptsnn_local_data_2 *const fd);
void findptsnn_local_2(
        uint   *const  code_base   , const unsigned  code_stride   ,
        uint   *const    el_base   , const unsigned    el_stride   ,
        double *const     r_base   , const unsigned     r_stride   ,
        double *const dist2_base   , const unsigned dist2_stride   ,
  const double *const     x_base[2], const unsigned     x_stride[2],
  const uint    *const  session_id_base, const unsigned session_id_stride,
        double *const disti_base   , const unsigned disti_stride   ,
  const uint npt, struct findptsnn_local_data_2 *const fd,
  buffer *buf);
void findptsnn_local_eval_2(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const double *const in, struct findptsnn_local_data_2 *const fd);
void findptsnn_local_el_eval_2(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,const uint elnum,
  const double *const in, struct findptsnn_local_data_2 *const fd);

#define findptsnn_local_setup_3   PREFIXED_NAME(findptsnn_local_setup_3)
#define findptsnn_local_free_3    PREFIXED_NAME(findptsnn_local_free_3 )
#define findptsnn_local_3         PREFIXED_NAME(findptsnn_local_3      )
#define findptsnn_local_eval_3    PREFIXED_NAME(findptsnn_local_eval_3 )
#define findptsnn_local_el_eval_3    PREFIXED_NAME(findptsnn_local_el_eval_3 )

struct findptsnn_local_hash_data_3 {
  uint hash_n;
  struct dbl_range bnd[3];
  double fac[3];
  uint *offset;
  uint max;
};

struct findptsnn_local_data_3 {
  unsigned ntot;
  const double *elx[3];
  const unsigned *nsid;
  struct obbox_3 *obb;
  struct findptsnn_local_hash_data_3 hd;
  struct findptsnn_el_data_3 fed;
  double tol;
  double *distint;
  const double *distfint; //This is field for distance for all GLL points
};

void findptsnn_local_setup_3(struct findptsnn_local_data_3 *const fd,
                           const double *const elx[3],
                           const unsigned *const nsid,
                           const double *const distfint,
                           const unsigned n[3], const uint nel,
                           const unsigned m[3], const double bbox_tol,
                           const uint max_hash_size,
                           const unsigned npt_max, const double newt_tol);
void findptsnn_local_free_3(struct findptsnn_local_data_3 *const fd);
void findptsnn_local_3(
        uint   *const  code_base   , const unsigned  code_stride   ,
        uint   *const    el_base   , const unsigned    el_stride   ,
        double *const     r_base   , const unsigned     r_stride   ,
        double *const dist2_base   , const unsigned dist2_stride   ,
  const double *const     x_base[3], const unsigned     x_stride[3],
  const uint    *const  session_id_base, const unsigned session_id_stride,
        double *const disti_base   , const unsigned disti_stride   ,
  const uint npt, struct findptsnn_local_data_3 *const fd,
  buffer *buf);
void findptsnn_local_eval_3(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const double *const in, struct findptsnn_local_data_3 *const fd);
void findptsnn_local_el_eval_3(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,const uint elnum,
  const double *const in, struct findptsnn_local_data_2 *const fd);

#endif
