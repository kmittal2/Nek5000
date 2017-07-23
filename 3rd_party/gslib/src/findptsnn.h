#ifndef FINDPTSNN_H
#define FINDPTSNN_H

#if !defined(COMM_H)
#warning "findptsnn.h" requires "comm.h"
#endif

#define findptsnn_setup_2   PREFIXED_NAME(findptsnn_setup_2)
#define findptsnn_free_2    PREFIXED_NAME(findptsnn_free_2 )
#define findptsnn_2         PREFIXED_NAME(findptsnn_2      )
#define findptsnn_eval_2    PREFIXED_NAME(findptsnn_eval_2 )
#define findptsnn_setup_3   PREFIXED_NAME(findptsnn_setup_3)
#define findptsnn_free_3    PREFIXED_NAME(findptsnn_free_3 )
#define findptsnn_3         PREFIXED_NAME(findptsnn_3      )
#define findptsnn_eval_3    PREFIXED_NAME(findptsnn_eval_3 )

struct findptsnn_data_2;
struct findptsnn_data_3;

struct findptsnn_data_2 *findptsnn_setup_2(
  const struct comm *const comm,
  const double *const elx[2],
  const unsigned n[2], const uint nel,
  const unsigned m[2], const double bbox_tol,
  const uint local_hash_size, const uint global_hash_size,
  const unsigned npt_max, const double newt_tol, 
  const uint *const nsid, const double *const distfint);

struct findptsnn_data_3 *findptsnn_setup_3(
  const struct comm *const comm,
  const double *const elx[3],
  const unsigned n[3], const uint nel,
  const unsigned m[3], const double bbox_tol,
  const uint local_hash_size, const uint global_hash_size,
  const unsigned npt_max, const double newt_tol,
  const uint *const nsid, const double *const distfint);

void findptsnn_free_2(struct findptsnn_data_2 *fd);
void findptsnn_free_3(struct findptsnn_data_3 *fd);

void findptsnn_2(    uint   *const  code_base   , const unsigned  code_stride   ,
                   uint   *const  proc_base   , const unsigned  proc_stride   ,
                   uint   *const    el_base   , const unsigned    el_stride   ,
                   double *const     r_base   , const unsigned     r_stride   ,
                   double *const dist2_base   , const unsigned dist2_stride   ,
             const double *const     x_base[2], const unsigned     x_stride[2],
    const uint  *const  session_id_base , const unsigned  session_id_stride   ,
                   double *const disti_base   , const unsigned disti_stride   ,
                   uint   *const elsid_base   , const unsigned elsid_stride   ,
             const uint npt, struct findptsnn_data_2 *const fd);

void findptsnn_3(    uint   *const  code_base   , const unsigned  code_stride   ,
                   uint   *const  proc_base   , const unsigned  proc_stride   ,
                   uint   *const    el_base   , const unsigned    el_stride   ,
                   double *const     r_base   , const unsigned     r_stride   ,
                   double *const dist2_base   , const unsigned dist2_stride   ,
             const double *const     x_base[3], const unsigned     x_stride[3],
  const uint *const  session_id_base   , const unsigned    session_id_stride   ,
                   double *const disti_base   , const unsigned disti_stride   ,
                   uint   *const elsid_base   , const unsigned elsid_stride   ,
             const uint npt, struct findptsnn_data_3 *const fd);

void findptsnn_eval_2(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const double *const in, struct findptsnn_data_2 *const fd);
 
void findptsnn_eval_3(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const double *const in, struct findptsnn_data_3 *const fd);

#endif
