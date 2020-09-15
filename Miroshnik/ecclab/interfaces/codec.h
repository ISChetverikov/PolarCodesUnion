//=============================================================================
// Codec interface.
//
// Copyright 2001 and onwards Kirill Shabunov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//=============================================================================

#ifndef CODEC_H

#define CODEC_H

//-----------------------------------------------------------------------------
// Defines

// #define DEC_NEEDS_SIGMA
// #define DEC_NEEDS_CSNRN

//-----------------------------------------------------------------------------
// Prototypes.

int
cdc_init(
   char param_str[],
   void **cdc
);

int cdc_get_n(void *cdc);
int cdc_get_k(void *cdc);
int cdc_get_flip_num(void *cdc);
double cdc_get_flip_alpha(void *cdc);

#ifdef DEC_NEEDS_CSNRN
void
cdc_set_csnrn(
   void *cdc,
   int csnrn
);
#endif // DEC_NEEDS_CSNRN

#ifdef DEC_NEEDS_SIGMA
void
cdc_set_sg(
   void *cdc,
   double noise_sg
);
#endif // DEC_NEEDS_SIGMA

int
enc_bpsk(
   void *cdc,
   int x[],
   double y[]
);

#ifdef GCCDEC
int
enc_bpsk_gcc(
   void *cdc,
   int x[],
   double y[],
   int *outer_code,
   int *inner_code,
   int out_cd_len,
   int in_cd_len
);
#endif // GCCDEC

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int xd[]
);

#ifdef FLIPPING
int
dec_bpsk_flipping(
   void *cdc,
   double c_out[],
   int xd[],
   int Ti
);
#endif // FLIPPING

#ifdef LISTFLIPPING
int
dec_bpsk_list_flipping(
   void *cdc,
   double c_out[],
   int xd[],
   int Ti,
   double alpha,
   int *var1,
   int *var2,
   int *var3
);
#endif // LISTFLIPPING

#ifdef LISTFLIPPINGPRECALC
int
dec_bpsk_list_flipping_precalc(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T
);
#endif

#ifdef LISTFLIPPINGFAST
int
dec_bpsk_list_flipping_fast(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T
);
#endif

#ifdef GCCDEC
int
dec_bpsk_gcc(
   void *cdc,
   double c_out[],
   int x_dec[],
   int *inner_matrix_inv,
   int in_cd_len,
   int out_cd_len,
   int *cdc_opts,
   int *cdc_sizes,
   int **cdc_mats,
   int ***grid4,
   double *last_part_true_codeword
);
#endif // GCCDEC

void
cdc_close(
   void *cdc
);

#endif // #ifndef CODEC_H