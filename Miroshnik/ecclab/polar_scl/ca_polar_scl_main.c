//=============================================================================
// CRC-aided (CA) Polar SCL codec.
// Implements codec.h
// Needs channel sigma (define global DEC_NEEDS_SIGMA).
//
// Copyright 2019 and onwards Kirill Shabunov.
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

//-----------------------------------------------------------------------------
// Includes.

#include <stdlib.h>
#include <string.h>
#include "../common/typedefs.h"
#include "../common/std_defs.h"
#include "../interfaces/codec.h"
#include "../common/spf_par.h"
#include "../interfaces/ui_utils.h"
#include "../common/srm_utils.h"
#include "polar_scl_inner.h"

//-----------------------------------------------------------------------------
// Internal defines.

// #define DBG

// Max # of SNR values.
#define VARL_NUM_MAX        50

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
   decoder_type dc;
   int c_n; // Code length.
   int c_k; // Code dimension.
   int eff_k; // Effective code dimension (without CRC).
   int c_m; // RM m.
   int *crc;
   int crc_len;
   int csnrn;
   #ifdef LISTFLIPPING
   int flip_num;
   double flip_alpha;
   #endif // LISTFLIPPING
   #if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
   int flip_num;
   #endif
   double sg22; // 2 / sg^2 for calculation of epsilons.
   int lsiz[VARL_NUM_MAX];
} cdc_inst_type;


//-----------------------------------------------------------------------------
// Global data.

//-----------------------------------------------------------------------------
// Internal functions.

int get_membuf_size(decoder_type *dd) {
   int i;
   int flsiz = MAX(dd->peak_lsiz, dd->p_num) * FLSIZ_MULT;
   int mem_buf_size = 0;

   mem_buf_size += dd->c_n * sizeof(ylitem); // y_in
   mem_buf_size += flsiz * sizeof(int); // frind
   mem_buf_size += flsiz * sizeof(int); // lorder
   mem_buf_size += (dd->c_k * flsiz + 1) * sizeof(xlist_item); // xlist
   mem_buf_size += flsiz * sizeof(xlist_item *); // lind2xl
   mem_buf_size += dd->c_m * flsiz * sizeof(ylitem *); // ylist
   mem_buf_size += dd->c_m * sizeof(ylitem *) + (dd->c_n - 1) * MAX(dd->peak_lsiz, dd->p_num) * 4 * sizeof(ylitem); // yitems
   mem_buf_size += flsiz * sizeof(slitem); // slist
   mem_buf_size += flsiz * sizeof(int); // plist.
   mem_buf_size += dd->c_n * sizeof(xlitem); // xtmp.
   mem_buf_size += flsiz * sizeof(int); // parent.
#ifdef FLIPPING
   mem_buf_size += flsiz * sizeof(lv) + flsiz * dd->c_n * sizeof(double); // lv_array_cur
#endif // FLIPPING
#ifdef LISTFLIPPING
  mem_buf_size += dd->c_n * sizeof(double); // Malpha_cur
  mem_buf_size += flsiz * sizeof(slitem); // slist_alpha
#endif // LISTFLIPPING
#ifdef LISTFLIPPINGOPT
  mem_buf_size += dd->c_k * sizeof(int); // cur_inf_word
#endif
   return mem_buf_size;
}

int get_k_from_node_table(
  int c_n,
  uint32 *node_table
) {
  int i, k = 0;
  for (i = 0; i < c_n; i++) {
    k += node_table[i];
  }
  return k;
}

int *malloc_calc_crc(
  int *x,
  int x_len,
  int *crc,
  int crc_len,
  int *tmp
) {
  int i, j;

  if (tmp == NULL) {
    tmp = (int *)malloc(x_len * sizeof(int));
    if (tmp == NULL) {
      return NULL;
    }
  }

  memcpy(tmp, x, x_len * sizeof(int));
  //memset(tmp + x_len, 0, crc_len * sizeof(int));

  for (i = 0; i < x_len - crc_len; i++) {
    if (tmp[i]) {
      for (j = 0; j < crc_len; j++) {
        tmp[i + j + 1] ^= crc[j];
      }
    }
  }

  return tmp;
}

int add_crc(
  int *x,
  int x_len,
  int *crc,
  int crc_len,
  int *check
) {
  int *tmp;

  if (crc_len <= 0) {
    return 0;
  }

  tmp = malloc_calc_crc(x, x_len, crc, crc_len, NULL);
  if (tmp == NULL) {
    return -1;
  }

  memcpy(check, tmp + x_len - crc_len, crc_len * sizeof(int));

  free(tmp);

  return 0;
}

int get_best_in_list(
  int *x,
  int x_len,
  int lsiz,
  int *crc,
  int crc_len
) {
  int *curx, *tmp = NULL;
  int s, best_ind;
  int sum_len = x_len + crc_len;
  int i, j;

  if (crc_len <= 0) {
    return 0;
  }

  best_ind = 0;
  for (i = 0; i < lsiz; i++) {
    curx = x + i * sum_len;
    tmp = malloc_calc_crc(curx, x_len, crc, crc_len, tmp);
    if (tmp == NULL) {
      return 0;
    }
    s = 0;
    for (j = 0; j < crc_len; j++) {
      s += curx[x_len + j] ^ tmp[x_len - crc_len + j];
#ifdef DBG
      printf("i: %d, s: %d\n", i, s);
#endif
    }
    if (s == 0) {
      best_ind = i;
      break;
    }
  }

  free(tmp);

  return best_ind;
}

#if defined LISTFLIPPINGOPT

#define INF (1e10)

int find_top_pos(double *min_abs_llrs, int max_diff_pos, int c_n) {
   double cur_llr;
   int i, pos, top_pos;
   top_pos = 0;
   while (1) {
      cur_llr = INF;
      pos = 0;
      for (i = 0; i < c_n; i++) {
         if (min_abs_llrs[i] != -1 && cur_llr > min_abs_llrs[i]) {
            cur_llr = min_abs_llrs[i];
            pos = i;
         }
      }
      min_abs_llrs[pos] = INF;
      if (max_diff_pos == pos) {
         return top_pos;
      }
      top_pos++;
   }
}

void find_lowest_llrs(double *min_abs_llrs, int c_n, int T, int* res) {
   double cur_llr;
   int i, pos, j;
   for (j = 0; j < T; j++) {
      cur_llr = INF;
      pos = 0;
      for (i = 0; i < c_n; i++) {
         if (min_abs_llrs[i] != -1 && cur_llr > min_abs_llrs[i]) {
            cur_llr = min_abs_llrs[i];
            pos = i;
         }
      }
      min_abs_llrs[pos] = INF;

      res[j] = pos;
   }
}

#if defined LLRPRECALCCOMP
void find_top_metrics(double *metrics, int c_n, int T, int* res) {
   double cur_metric;
   int i, pos, j;
   for (j = 0; j < T; j++) {
      cur_metric = -INF;
      pos = 0;
      for (i = 0; i < c_n; i++) {
         if (cur_metric < metrics[i]) {
            cur_metric = metrics[i];
            pos = i;
         }
      }
      metrics[pos] = -INF;
      res[j] = pos;
   }
}
#endif

#if defined LISTFLIPPINGTHRESHOLD
void find_lowest_llrs_threshold(double *min_abs_llrs, int c_n, int T, int* res, double threshold) {
   double cur_llr;
   int i, pos, j, cnt_success;
   cnt_success = 0;
   for (j = 0; j < T; j++) {
      cur_llr = INF;
      pos = -1;
      for (i = 0; i < c_n; i++) {
         if (min_abs_llrs[i] != -1 && cur_llr > min_abs_llrs[i]) {
            cur_llr = min_abs_llrs[i];
            pos = i;
         }
      }
      if (min_abs_llrs[pos] > threshold || pos == -1) break; 
      min_abs_llrs[pos] = INF;

      res[j] = pos;
      cnt_success++;
   }
   for (i = cnt_success; i < T; i++) res[i] = -1;
}
#endif

#endif

#if defined LISTFLIPPING || defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
int check_crc_list(
  int *x,
  int x_len,
  int lsiz,
  int *crc,
  int crc_len
) {
  int *curx, *tmp = NULL;
  int s, best_ind;
  int sum_len = x_len + crc_len;
  int i, j;

  if (crc_len <= 0) {
    return -1;
  }

  best_ind = -1;
  for (i = 0; i < lsiz; i++) {
    curx = x + i * sum_len;
    tmp = malloc_calc_crc(curx, x_len, crc, crc_len, tmp);
    if (tmp == NULL) {
      return -1;
    }
    s = 0;
    for (j = 0; j < crc_len; j++) {
      s += curx[x_len + j] ^ tmp[x_len - crc_len + j];
    }
    if (s == 0) {
      best_ind = i;
      break;
    }
  }

  free(tmp);

  return best_ind;
}
#endif 

#ifdef LISTFLIPPING
void find_alpha_indices(double *Malpha, int* node_table, int c_n, int T, int *res) {
   double *tmp;
   int i, j, max_ind;
   double cur_max, min_alpha;
   tmp = (double *) malloc(sizeof(double) * c_n);
   memcpy(tmp, Malpha, sizeof(double) * c_n);
   min_alpha = tmp[0];
   for (i = 1; i < c_n; i++) {
      min_alpha = fmin(min_alpha, tmp[i]);
   }
   for (j = 0; j < T; j++) {
      cur_max = min_alpha-1;
      max_ind = 0;
      for (i = 0; i < c_n; i++) {
         if (node_table[i] && cur_max < tmp[i]) {
            cur_max = tmp[i];
            max_ind = i;
         }
      }
      res[j] = max_ind;
      tmp[max_ind] = min_alpha - 1; 
   }
   free(tmp);
   return;
}
#endif // LISTFLIPPING

#ifdef LISTFLIPPINGFAST
void find_subcodes_indices(int *subcodes_first_symb, double *subcodes_min_abs_llr, int* node_table, int c_n, int T, int *res) {
   double *tmp;
   int i, j, min_ind, len;
   double cur_min, max_llr;
   for (i = 0; i < c_n; i++) {
      if (subcodes_first_symb[i] == -2) {
         len = i;
         break;
      }
   }
   tmp = (double *) malloc(sizeof(double) * len);
   memcpy(tmp, subcodes_min_abs_llr, sizeof(double) * len);
   max_llr = tmp[0];
   for (i = 1; i < len; i++) {
      max_llr = fmax(max_llr, tmp[i]);
   }
   j = 0;
   while (j < T) {
      cur_min = max_llr + 1;
      min_ind = 0;
      for (i = 0; i < len; i++) {
         if (cur_min > tmp[i] && subcodes_first_symb[i] != -1) {
            cur_min = tmp[i];
            min_ind = i;
         }
      }
      res[j] = subcodes_first_symb[min_ind];
      tmp[min_ind] = max_llr + 1; 
      j++;
   }
   free(tmp);
   return;
}
#endif

#ifdef FLIPPING
int check_crc(
  int *x,
  int x_len,
  int lsiz,
  int *crc,
  int crc_len
) {
  int *tmp = NULL;
  int s;
  int i, j;

  if (crc_len <= 0) {
    return 0;
  }
  tmp = malloc_calc_crc(x, x_len, crc, crc_len, tmp);
  if (tmp == NULL) {
    return 1;
  }
  s = 0;
  for (j = 0; j < crc_len; j++) {
    s += x[x_len + j] ^ tmp[x_len - crc_len + j];
  }
  free(tmp);
  if (s == 0) {
    return 0;
  }
  return 1;
}

void find_lr_indices(double *llrs, int* node_table, int c_n, int T, int *res) {
   double *tmp;
   int i, j, min_ind;
   double max_llr, cur_min;
   tmp = (double *) malloc(sizeof(double) * c_n);
   memcpy(tmp, llrs, sizeof(double) * c_n);
   max_llr = 0;
   for (i = 0; i < c_n; i++) {
      if (node_table[i]) {
         max_llr = fmax(fabs(tmp[i]), max_llr);
      }
   }
   for (i = 0; i < c_n; i++) {
      if (!node_table[i]) {
         tmp[i] = max_llr + 1;
      }
   }
   for (j = 0; j < T; j++) {
      cur_min = max_llr;
      min_ind = 0;
      for (i = 0; i < c_n; i++) {
         if (node_table[i] && cur_min > fabs(tmp[i])) {
            cur_min = fabs(tmp[i]);
            min_ind = i;
         }
      }
      res[j] = min_ind;
      tmp[min_ind] = max_llr + 1; 
   }
   free(tmp);
   return;
}
#endif // FLIPPING

#ifdef DBG
void print_bin_vect(char rem[], int *x, int x_len) {
  int i;
  printf("%s: ", rem);
  for (i = 0; i < x_len; i++) {
    printf("%d", x[i]);
  }
  printf("\n");
}
#endif

#ifdef RETURNLLRS
void print_dec_info(cdc_inst_type *dd, int *x_candidates, lv *lv_array) {
   int cnt, i, j;
   printf("x_candidates:\n");
   for (i = 0; i < dd->dc.peak_lsiz; i++) {
      cnt = 0;
      j = 0;
      printf("Number %d\n", i+1);
      printf("L: %f\n", lv_array[i].l);
      printf("x, llrs:\n");
      while (cnt < dd->c_n) {
         printf("%d ", cnt);
         if (dd->dc.node_table[cnt]) {
            printf("%d ", x_candidates[i * dd->dc.peak_lsiz + j]);
            printf("%f \n", lv_array[i].llrs[cnt]);
            j++;
         }
         else {
            printf("f f\n");
         }
         cnt++;
      }
   }
}
#endif // RETURNLLRS

//-----------------------------------------------------------------------------
// codec interface functions.

int
cdc_init(
   char param_str[],
   void **cdc
)
{
   cdc_inst_type *dd = NULL;
   char *str1 = NULL, *token;
   int c_n = 0, c_k = 0, c_m = 0;
   int i, j, i1, i2, pn;
   int *p1, *px, *py;
   int p_info[10000], inf_coeff[2048];
   int p_num = 0;

   // Allocate decoder instance.
   dd = (cdc_inst_type *)malloc(sizeof(cdc_inst_type));
   if (dd == NULL) {
      show_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   memset(dd, 0, sizeof(cdc_inst_type));

   // Allocate temporary string for parsing.
   str1 = (char *)malloc(SP_STR_MAX);
   if (str1 == NULL) {
      err_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   str1[0] = 0;

   // ----- Parse parameters.
   strcpy(str1, param_str);

   token = strtok(str1, tk_seps_prepared);
   while (token != NULL) {

      TRYGET_INT_TOKEN(token, "c_m", c_m);
      if (strcmp(token, "info_bits_mask") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.node_table_len = strlen(token);
         dd->dc.node_table = (uint32 *)malloc(dd->dc.node_table_len * 4);
         if (dd->dc.node_table == NULL) {
            err_msg("cdc_init: Short of memory.");
            goto ret_err;
         }
         for (i = 0; i < dd->dc.node_table_len; i++)
            dd->dc.node_table[i] = (token[i] == '0') ? 0 : 1;
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "ca_polar_crc") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->crc_len = strlen(token);
         //dd->crc_len = strlen(token) + 1;
         dd->crc = (int *)malloc(dd->crc_len * sizeof(int));
         if (dd->crc == NULL) {
            err_msg("cdc_init: Short of memory.");
            goto ret_err;
         }
         for (i = 0; i < dd->crc_len; i++)
            dd->crc[i] = (token[i] == '0') ? 0 : 1;
         //for (i = 0; i < dd->crc_len-1; i++)
         //   dd->crc[i] = (token[i] == '0') ? 0 : 1;
         //dd->crc[dd->crc_len - 1] = 1;
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "list_size") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.peak_lsiz = atoi(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "permutations") == 0) {
         // TODO: Permutations are not usable with general Polar, but let's leave it here for now.
         if (c_m == 0) {
            err_msg("cdc_init(): parameters file: c_m must be specified before permutations.");
            goto ret_err;
         }
         token = strtok(NULL, tk_seps_prepared);
         if (token[0] == GROUP_CHAR_BEGIN) {
            token = strtok(NULL, tk_seps_prepared);
            i1 = 0;
            p_num = 0;
            while (token[0] != GROUP_CHAR_END) {
               for (i = 0; i < c_m; i++) {
                  p_info[i1++] = atoi(token);
                  token = strtok(NULL, tk_seps_prepared);
                  if (token == NULL) {
                     err_msg("cdc_init(): parameters file: permutations.");
                     goto ret_err;
                  }
               }
               p_num++;
            }
         }
         else {
            err_msg("cdc_init(): parameters file: permutations.");
            goto ret_err;
         }
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }

      #ifdef LISTFLIPPING
      if (strcmp(token, "flips") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.flip_num = atoi(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }

      if (strcmp(token, "alpha") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.flip_alpha = atof(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      #endif // LISTFLIPPING

      #if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
      if (strcmp(token, "flips") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.flip_num = atoi(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      #endif

      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   free(str1); str1 = NULL;

   dd->c_m = dd->dc.c_m = c_m;

   // Check if all necessary parameters are specified.
   if (c_m == 0) {
      err_msg("cdc_init: c_m is not specified.");
      goto ret_err;
   }
   if (dd->dc.node_table_len == 0) {
      err_msg("cdc_init: No info_bits_mask.");
      goto ret_err;
   }
   if (dd->dc.peak_lsiz == 0) {
      err_msg("cdc_init: No list size.");
      goto ret_err;
   }

   // Complete settings.
   dd->c_n = dd->dc.c_n = c_n = 1 << c_m;
   dd->c_k = dd->dc.c_k = c_k = get_k_from_node_table(c_n, dd->dc.node_table);
   dd->eff_k = c_k - dd->crc_len;
   //dd->eff_k = c_k - dd->crc_len + 1;
   dd->dc.ret_list = (dd->crc_len > 0) ? 1 : 0;
   for (i = 0; i < dd->dc.node_table_len; i++)
      dd->dc.node_table[i] = (dd->dc.node_table[i] == 0) ? 0 : dd->dc.peak_lsiz;

   #ifdef LISTFLIPPING
   dd->flip_num = dd->dc.flip_num;
   dd->flip_alpha = dd->dc.flip_alpha;
   #endif // LISTFLIPPING

   #if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
   dd->flip_num = dd->dc.flip_num;
   #endif

   // Generate permutations.
   // TODO: Permutations.
   if (p_num == 0) {
      // permutations parameter was not given => only one permutation
      // (0 1 ... m-1) (no permutation).
      for (i = 0; i < c_m; i++) p_info[i] = i;
      p_num = 1;
   }
   dd->dc.p_num = p_num;
   if (p_num > 0) {
      dd->dc.pxarr = (int *)malloc(p_num * c_k * sizeof(int));
      dd->dc.pyarr = (int *)malloc(p_num * c_n * sizeof(int));
      if ((dd->dc.pxarr == NULL) || (dd->dc.pyarr == NULL)) {
         err_msg("cdc_init: Short of memory.");
         goto ret_err;
      }
      px = dd->dc.pxarr;
      py = dd->dc.pyarr;
      p1 = p_info;
   }
   // Calculate inf_coeff[]. It contains all integers with weight <= c_r
   // in reverse lexicographic order.
   i1 = 0;
   for (i = c_n - 1; i >= 0; i--) {
      // i2 <-- weight(i);
      i2 = i & 1;
      for (j = 1; j < c_m; j++) i2 += (i >> j) & 1;
      if (i2 <= c_m) inf_coeff[i1++] = i;
   }
   // Generate permutation arrays.
   for (pn = 0; pn < p_num; pn++) {
      // Next permutation for y.
      for (i = 0; i < c_n; i++) {
         i1 = 0;
         for (j = 0; j < c_m; j++) i1 |= ((i >> j) & 1) << p1[j];
         py[i] = i1;
      }
      // for (i = 0; i < c_n; i++) msg_printf("%d ", py[i]); // DDD
      // msg_printf("\n"); // DDD
      py += c_n;
      // Next permutation for x.
      for (i = 0; i < c_k; i++) {
         i1 = 0;
         i2 = inf_coeff[i];
         for (j = 0; j < c_m; j++) i1 |= ((i2 >> j) & 1) << p1[j];
         for (j = 0; j < c_k; j++) if (inf_coeff[j] == i1) break;
         px[i] = j;
      }
      // for (i = 0; i < c_k; i++) msg_printf("%d ", px[i]); // DDD
      // msg_printf("\n"); // DDD
      px += c_k;
      p1 += c_m;
   }

#ifdef FORMAT_Q
   if (dec_tables_init()) {
      err_msg("cdc_init: error in dec_tables_init().");
      goto ret_err;
   }
#endif // FORMAT_Q

   // Allocate decoder memory buffer
   dd->dc.mem_buf = (uint8 *)malloc(get_membuf_size(&(dd->dc)));
   if (dd->dc.mem_buf == NULL) {
      err_msg("cdc_init: Short of memory.");
      goto ret_err;
   }

   (*cdc) = dd;

   return 0;

ret_err:

   if (str1 != NULL) free(str1);
   if (dd != NULL) {
      if (dd->dc.mem_buf != NULL) free(dd->dc.mem_buf);
      if (dd->dc.pxarr != NULL) free(dd->dc.pxarr);
      if (dd->dc.pyarr != NULL) free(dd->dc.pyarr);
      if (dd->dc.node_table != NULL) free(dd->dc.node_table);
      free(dd);
      (*cdc) = NULL;
   }
   return 1;
}

void
cdc_close(
   void *cdc
)
{
   cdc_inst_type *dd;

   if (cdc == NULL) return;
   dd = (cdc_inst_type *)cdc;
   if (dd->dc.mem_buf != NULL) free(dd->dc.mem_buf);
   if (dd->crc != NULL) free(dd->crc);
   if (dd->dc.node_table != NULL) free(dd->dc.node_table);
   if (dd->dc.pxarr != NULL) free(dd->dc.pxarr);
   if (dd->dc.pyarr != NULL) free(dd->dc.pyarr);
   free(dd);
#ifdef FORMAT_Q
   dec_tables_free();
#endif // FORMAT_Q
}

int cdc_get_n(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->c_n;
}

int cdc_get_k(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->eff_k;
}

#if defined LISTFLIPPING || defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
int cdc_get_flip_num(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->flip_num;
}
#endif

#ifdef LISTFLIPPING
double cdc_get_flip_alpha(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->flip_alpha;
}
#endif // LISTFLIPPING

void
cdc_set_sg(
   void *cdc,
   double noise_sg
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   dd->sg22 = 2.0 / (noise_sg * noise_sg);
}

int
enc_bpsk(
   void *cdc,
   int x[],
   double y[]
)
{
   cdc_inst_type *dd;
   int *x_crc;

   dd = (cdc_inst_type *)cdc;

   if (dd->crc_len > 0) {
      x_crc = (int *)malloc(dd->c_k * sizeof(int));
      if (x_crc == NULL) {
         err_msg("enc_bpsk: Short of memory.");
         return -1;
      }
      memcpy(x_crc, x, dd->eff_k * sizeof(int));
      add_crc(x, dd->eff_k, dd->crc, dd->crc_len, x_crc + dd->eff_k);
#ifdef DBG
      print_bin_vect("Input x", x, dd->eff_k);
      print_bin_vect("Input x with crc", x_crc, dd->c_k);
#endif
      smrm_par0_enc_bpsk(dd->c_m, dd->c_m, dd->dc.node_table, x_crc, y);
      free(x_crc);
      return dd->eff_k;
   }
   else {
      return smrm_par0_enc_bpsk(dd->c_m, dd->c_m, dd->dc.node_table, x, y);
   }
}

#if defined GCCDEC
int
enc_bpsk_gcc(
   void *cdc,
   int x[],
   double y[],
   int *outer_matrix,
   int *inner_matrix,
   int out_cd_len,
   int in_cd_len
)
{
   cdc_inst_type *dd;
   int i, j, k, in_vecs_len, node_ind;
   int *in_vecs, *res;
   
   in_vecs_len = out_cd_len / in_cd_len;
   dd = (cdc_inst_type *)cdc;
   in_vecs = (int *)malloc(sizeof(int) * in_cd_len * in_vecs_len);
   res = (int *)malloc(sizeof(int) * out_cd_len);
   memset(in_vecs, 0, sizeof(int) * in_cd_len * in_vecs_len);
   memset(res, 0, sizeof(int) * out_cd_len);

   node_ind = 0;

   if (dd->crc_len > 0) {
      return -1;
   }
   else {
      for (i = 0; i < in_cd_len; i++) {
         for (j = 0; j < in_vecs_len; j++) {
            if (!dd->dc.node_table[i * in_vecs_len + j]) continue;
            for (k = 0; k < in_vecs_len; k++) {
               in_vecs[i*in_vecs_len + k] ^= x[node_ind] * outer_matrix[(i*in_vecs_len + j) * out_cd_len + k];
            }
            node_ind++;
         }
      }
      for (i = 0; i < in_vecs_len; i++) {
         for (j = 0; j < in_cd_len; j++) {
            for (k = 0; k < in_cd_len; k++) {
               res[i + j * in_vecs_len] ^= in_vecs[k * in_vecs_len + i] * inner_matrix[k * in_cd_len + j];
            }
         }
      }
   }
   for (i = 0; i < out_cd_len; i++) {
      if (res[i]) y[i] = 1.0;
      else y[i] = -1.0;
   }
   free(in_vecs);
   free(res);
   return 0;
}
#endif // GCCDEC

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int x_dec[]
)
{
   cdc_inst_type *dd;
   double *dec_input;
   int *x_candidates;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind;
   double d1;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   if (dd->dc.ret_list) {
      mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   }
   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   if (dd->dc.ret_list) {
      x_candidates = (int *)mem_buf_ptr;
      mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);
   }

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;

   if (dd->dc.ret_list) {
      polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, -1, NULL, NULL, NULL, NULL, NULL);
      best_ind = get_best_in_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
#ifdef DBG
      printf("best_ind: %d\n", best_ind);
#endif
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      polar_dec(&(dd->dc), dec_input, x_dec, NULL, NULL, -1, -1, -1, NULL, NULL, NULL, NULL, NULL);
   }

   free(mem_buf);

   return 0;
}

#ifdef FLIPPING
int
dec_bpsk_flipping(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T
)
{
   cdc_inst_type *dd;
   double *dec_input;
   int *x_candidates, *lrdec;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j;
   double d1;
   int flsiz;
   lv *lv_array = NULL;

   dd = (cdc_inst_type *)cdc;

   if (dd->dc.peak_lsiz != 1) {
      err_msg("L != 1 for flipping");
      return -1;
   }

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   if (dd->dc.ret_list) {
      mem_buf_size += dd->c_k * sizeof(int); // x_candidates
   }
   flsiz = MAX(1, dd->dc.p_num) * FLSIZ_MULT;
   mem_buf_size += sizeof(lv) + dd->c_n * sizeof(double); // lv
   mem_buf_size += T * sizeof(int);
   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   if (dd->dc.ret_list) {
      x_candidates = (int *)mem_buf_ptr;
      mem_buf_ptr += dd->c_k * sizeof(int);
   }

   lv_array = (lv *)mem_buf_ptr;
   mem_buf_ptr += sizeof(lv);
   lv_array[0].llrs = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);
   lrdec = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   polar_dec(&(dd->dc), dec_input, x_dec, NULL, lv_array, -1, -1, -1, NULL, NULL, NULL, NULL, NULL);
#ifdef RETURNLLRS
   print_dec_info(dd, x_dec, lv_array);
#endif // RETURNLLRS
   if (T > 0 && check_crc(x_dec, dd->eff_k, 1, dd->crc, dd->crc_len)) {
      find_lr_indices(lv_array->llrs, dd->dc.node_table, dd->c_n, T, lrdec);
      for (i = 0; i < T; i++) {
         polar_dec(&(dd->dc), dec_input, x_dec, NULL, lv_array, lrdec[i], -1, -1, NULL, NULL, NULL, NULL, NULL);
         if (!check_crc(x_dec, dd->eff_k, 1, dd->crc, dd->crc_len)) {
            break;
         }
      }
   }
#ifdef RETURNLLRS
   print_dec_info(dd, x_dec, lv_array);
#endif // RETURNLLRS

   free(mem_buf);

   return 0;
}
#endif // FLIPPING

#ifdef LISTFLIPPING
int
dec_bpsk_list_flipping(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T,
   double alpha,
   int *var1,
   int *var2,
   int *var3
)
{
   cdc_inst_type *dd;
   double *dec_input, *Malpha;
   int *x_candidates, *Malpha_maxs;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind;
   double d1;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += dd->c_n * sizeof(double); // Malpha
   mem_buf_size += T * sizeof(int); // Malpha_maxs

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   Malpha = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, alpha, Malpha, NULL, NULL, NULL, NULL);
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      find_alpha_indices(Malpha, dd->dc.node_table, dd->c_n, T, Malpha_maxs);
      for (i = 0; i < T; i++) {
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, Malpha_maxs[i], alpha, Malpha, NULL, NULL, NULL, NULL);
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {

            for (j = 0; j < T; j++) {
               var1[Malpha_maxs[j]]++;
            }
            for (j = 0; j < i; j++) {
               var2[Malpha_maxs[j]]++;
            }
            var3[Malpha_maxs[i]]++;

            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}
 
#if defined LISTFLIPPINGOPT
int
dec_bpsk_list_flipping_opt(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T,
   double alpha,
   int *var1,
   int *var2,
   int *var3,
   int *x
)
{
   cdc_inst_type *dd;
   double *dec_input, *Malpha, *min_abs_llrs;
   int *x_candidates, *Malpha_maxs;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind, max_diff_pos, max_diff_pos_c, q;
   int x0;
   double d1;
   #ifdef PRECALCCOVER
   FILE *precalc_bits;
   int freq;
   #endif
   #ifdef LLRCOVER
   int top_pos;
   #endif

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += dd->c_n * sizeof(double); // Malpha
   mem_buf_size += T * sizeof(int); // Malpha_maxs
   mem_buf_size += dd->c_n * sizeof(double); // min_abs_llrs

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   Malpha = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   min_abs_llrs = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   for (i = 0; i < dd->eff_k; i++) x_dec[i] = 0;

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   x0 = x[0];
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, alpha, Malpha, NULL, NULL, x, min_abs_llrs);
   max_diff_pos = x[0];
   x[0] = x0;
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      #ifdef PRECALCCOVER
      precalc_bits = fopen("./precalc_bits/precalc_polar8_ca24_k128_L8.spf", "r");
      i = 0;
      while (!feof (precalc_bits)) {
         if (i == T) break;
         fscanf(precalc_bits, "%d %d", &Malpha_maxs[i], &freq);
         i++;
      }
      fclose(precalc_bits);
      for (i = 0; i < dd->c_n; i++) {
         if (Malpha_maxs[i] == max_diff_pos) {
            var1[i]++;
            break;
         }
      }
      #endif
      #ifdef LLRCOVER 
      top_pos = find_top_pos(min_abs_llrs, max_diff_pos, dd->c_n);
      var2[top_pos]++;
      #endif
      //T = 1;
      T = 0;
      Malpha_maxs[0] = max_diff_pos;
      for (q = 0; q < T; q++) {
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, Malpha_maxs[q], alpha, Malpha, NULL, NULL, x, min_abs_llrs);
         x[0] = x0;
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {
            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}

#if defined LLRPRECALCCOMP
int
dec_bpsk_list_flipping_comp(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T,
   double alpha,
   int *var1,
   int *var2,
   int *var3,
   int *x
)
{
   cdc_inst_type *dd;
   double *dec_input, *Malpha, *min_abs_llrs;
   int *x_candidates, *Malpha_maxs;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind, max_diff_pos, max_diff_pos_c, q;
   int x0;
   double d1;
   FILE *precalc_bits, *precalc_metric, *llr_metric;
   int freq;
   double *Malpha_maxs_metric, *min_abs_llrs_metric, *comp_metric;
   double a = 0.9;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += dd->c_n * sizeof(double); // Malpha
   mem_buf_size += dd->c_n * sizeof(int); // Malpha_maxs
   mem_buf_size += dd->c_n * sizeof(double); // min_abs_llrs
   mem_buf_size += dd->c_n * sizeof(double); // Malpha_maxs_metric
   mem_buf_size += dd->c_n * sizeof(double); // min_abs_llrs_metric
   mem_buf_size += dd->c_n * sizeof(double); // comp_metric

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   Malpha = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(int);

   min_abs_llrs = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs_metric = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   min_abs_llrs_metric = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   comp_metric = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   for (i = 0; i < dd->eff_k; i++) x_dec[i] = 0;

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   x0 = x[0];
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, alpha, Malpha, NULL, NULL, x, min_abs_llrs);
   max_diff_pos = x[0];
   x[0] = x0;
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      memset(comp_metric, 0, dd->c_n * sizeof(double));
      precalc_bits = fopen("./precalc_bits/precalc_polar7_ca24_k64_L2_upd2.spf", "r");
      i = 0;
      while (!feof (precalc_bits)) {
         fscanf(precalc_bits, "%d %d", &Malpha_maxs[i], &freq);
         i++;
      }
      fclose(precalc_bits);

      //for (i = 0; i < T; i++) printf("%d %d\n", i, Malpha_maxs[i]);

      for (i = 0; i < dd->c_n; i++) {
         comp_metric[Malpha_maxs[i]] += a * (dd->c_n - i);
      }

      find_lowest_llrs_threshold(min_abs_llrs, dd->c_n, dd->c_n, Malpha_maxs, 1e6);

      //for (i = 0; i < T; i++) printf("%d %d\n", i, Malpha_maxs[i]);

      for (i = 0; i < dd->c_n; i++) {
         if (Malpha_maxs[i] == -1) break;
         comp_metric[Malpha_maxs[i]] += (1 - a) * (dd->c_n - i);
      }

      find_top_metrics(comp_metric, dd->c_n, T, Malpha_maxs);

      //for (i = 0; i < T; i++) printf("%d %d\n", i, Malpha_maxs[i]);

      for (q = 0; q < T; q++) {
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, Malpha_maxs[q], alpha, Malpha, NULL, NULL, x, min_abs_llrs);
         x[0] = x0;
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {
            //printf("%d %d\n", q, Malpha_maxs[q]);
            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}
#endif

#if defined LISTFLIPPINGTHRESHOLD
int
dec_bpsk_list_flipping_threshold(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T,
   double alpha,
   int *var1,
   double *var2,
   double *var3,
   int *var4,
   int *x,
   double threshold
)
{
   cdc_inst_type *dd;
   double *dec_input, *Malpha, *min_abs_llrs;
   int *x_candidates, *Malpha_maxs;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind, max_diff_pos, max_diff_pos_c, q;
   int x0, top_pos, freq;
   double d1;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += dd->c_n * sizeof(double); // Malpha
   mem_buf_size += T * sizeof(int); // Malpha_maxs
   mem_buf_size += dd->c_n * sizeof(double); // min_abs_llrs

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   Malpha = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   min_abs_llrs = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   for (i = 0; i < dd->eff_k; i++) x_dec[i] = 0;

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   x0 = x[0];
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, alpha, Malpha, NULL, NULL, x, min_abs_llrs);
   max_diff_pos = x[0];
   x[0] = x0;
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (max_diff_pos == -1) {
      for (q = 0; q < dd->c_n; q++) {
         var2[q] += min_abs_llrs[q];
      }
   }
   else {
      for (q = 0; q < dd->c_n; q++) {
         if (q == max_diff_pos) {
            var3[q] += min_abs_llrs[q];
            var4[q]++;
            break;
         }
      }
   }
   find_lowest_llrs_threshold(min_abs_llrs, dd->c_n, T, Malpha_maxs, threshold);
   for (q = 0; q < dd->c_n; q++) {
      if (max_diff_pos == -1) break;
      if (Malpha_maxs[q] == max_diff_pos) {
         var1[q]++;
         break;
      }
   }
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      for (q = 0; q < T; q++) {
         if (Malpha_maxs[q] == -1) {
            break;
         }
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, Malpha_maxs[q], alpha, Malpha, NULL, NULL, x, min_abs_llrs);
         x[0] = x0;
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {
            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}
#endif
#endif
#endif // LISTFLIPPING

#ifdef LISTFLIPPINGPRECALC
int
dec_bpsk_list_flipping_precalc(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T
)
{
   cdc_inst_type *dd;
   double *dec_input, *Malpha;
   int *x_candidates, *Malpha_maxs;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind, freq;
   double d1;
   FILE* precalc_bits;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += dd->c_n * sizeof(double); // Malpha
   mem_buf_size += T * sizeof(int); // Malpha_maxs

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   Malpha = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   Malpha_maxs = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, -1, NULL, NULL, NULL, NULL, NULL);
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      precalc_bits = fopen("./precalc_bits/precalc_polar7_ca9_k64_L4_upd2.spf", "r");
      i = 0;
      while (!feof (precalc_bits)) {
         if (i == T) break;
         fscanf(precalc_bits, "%d %d", &Malpha_maxs[i], &freq);
         i++;
      }
      fclose(precalc_bits);
      for (i = 0; i < T; i++) {
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, Malpha_maxs[i], -1, NULL, NULL, NULL, NULL, NULL);
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {
            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}
#endif // LISTFLIPPINGPRECALC

#ifdef LISTFLIPPINGFAST
int
dec_bpsk_list_flipping_fast(
   void *cdc,
   double c_out[],
   int x_dec[],
   int T
)
{
   cdc_inst_type *dd;
   double *dec_input;
   int *x_candidates, *subcodes_mins;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int i, j, best_ind;
   double d1;
   int *subcodes_first_symb;
   double *subcodes_min_abs_llr;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   mem_buf_size += T * sizeof(int); // subcodes_mins
   mem_buf_size += dd->c_n * sizeof(int); // subcodes_first_symb
   mem_buf_size += dd->c_n * sizeof(double); // subcodes_min_abs_llr

   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   x_candidates = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);

   subcodes_mins = (int *)mem_buf_ptr;
   mem_buf_ptr += T * sizeof(int);

   subcodes_first_symb = (int *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(int);

   subcodes_min_abs_llr = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, -1, -1, NULL, subcodes_first_symb, subcodes_min_abs_llr, NULL, NULL);
   best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
   if (best_ind != -1) {
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      find_subcodes_indices(subcodes_first_symb, subcodes_min_abs_llr, dd->dc.node_table, dd->c_n, T, subcodes_mins);
      //for (i = 0; i < T; i++) printf("%d %d\n", i, subcodes_mins[i]);
      for (i = 0; i < T; i++) {
         polar_dec(&(dd->dc), dec_input, x_candidates, NULL, NULL, -1, subcodes_mins[i], -1, NULL, subcodes_first_symb, subcodes_min_abs_llr, NULL, NULL);
         best_ind = check_crc_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
         if (best_ind != -1) {
            memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
            break;
         }
      }
   }

   free(mem_buf);

   return 0;
}
#endif // LISTFLIPPINGFAST

#ifdef GCCDEC

#define INF (1e10)

int MLDec(double *llrs, int llrs_len, int *node_table, int *cdc_size, int *gen_mat, int *x_dec, int *ret_codeword) {
   double best_dist, cur_dist;
   int i, j, k, word, *codeword, best_word, cur_dec_ind;
   codeword = (int *)malloc(sizeof(int) * cdc_size[1]);

   best_dist = -INF;
   word = 0, best_word = 0;

   for (i = 0; i < pow(2, cdc_size[0]); i++) {
      memset(codeword, 0, sizeof(int) * cdc_size[1]);
      for (j = 0; j < cdc_size[0]; j++) {
         if (((word >> j) & 1)) {
            for (k = 0; k < cdc_size[1]; k++) {
               codeword[k] ^= gen_mat[j * cdc_size[1] + k];
            }
         }
      }

      cur_dist = 0;
      for (j = 0; j < llrs_len; j++) {
         cur_dist += -(2 * codeword[j] - 1) * llrs[j];
      }
      if (cur_dist > best_dist) {
         best_dist = cur_dist;
         best_word = word;
      }
      word++;
   }
   
   memset(codeword, 0, sizeof(int) * cdc_size[1]);
   for (j = 0; j < cdc_size[0]; j++) {
      if (((best_word >> j) & 1)) {
         for (k = 0; k < cdc_size[1]; k++) {
            codeword[k] ^= gen_mat[j * cdc_size[1] + k];
         }
      }
   }

   for (k = 0; k < cdc_size[1]; k++) {
      ret_codeword[k] = codeword[k];
   }

   for (i = 0; i < cdc_size[0]; i++) {
      x_dec[i] = ((best_word >> i) & 1);
   }

   free(codeword);
   return 0;
}

int SyndrDec(double *llrs, int llrs_len, int *node_table, int *cdc_size, int *gen_mat, int *x_dec, int ***grid) {
   double **grid_metric;
   int i, j, k, cur_symb, *dec_codeword;
   dec_codeword = (int *)malloc(sizeof(int) * llrs_len);
   grid_metric = (double **)malloc(sizeof(double *) * 4);
   for (i = 0; i < 4; i++) {
      grid_metric[i] = (double *)malloc(sizeof(double) * 33);
      for (j = 0; j < 33; j++) {
         grid_metric[i][j] = -INF;
      }
   }
   grid_metric[0][0] = 0;
   for (i = 1; i < 33; i++) {
      for (j = 0; j < 4; j++) {
         if (grid[j][i][4]) {
            for (k = 0; k < 4; k++) {
               if (grid[j][i][k] && k == j) {
                  grid_metric[j][i] = fmax(grid_metric[j][i], grid_metric[k][i-1]+llrs[i-1]);
               }
               else if (grid[j][i][k] && k != j) {
                  grid_metric[j][i] = fmax(grid_metric[j][i], grid_metric[k][i-1]-llrs[i-1]);
               }

            }
         }
      }
   }
   cur_symb = 0;
   for (i = 32; i > 0; i--) {
      for (k = 0; k < 4; k++) {
         if (grid[cur_symb][i][k]) {
            if (k == cur_symb && grid_metric[k][i-1] + llrs[i-1] == grid_metric[cur_symb][i]) {
               dec_codeword[i - 1] = 0;
               break;
            }
            else if (k != cur_symb && grid_metric[k][i-1] - llrs[i-1] == grid_metric[cur_symb][i]) {
               dec_codeword[i - 1] = 1;
               cur_symb = k;
               break;
            }
         }
      }
   }

   for (i = 0; i < 32; i++) {
      x_dec[i] = dec_codeword[i];
   }

   free(dec_codeword);
   for (i = 0; i < 4; i++) {
      free(grid_metric[i]);
   }
   free(grid_metric);
   
   return 0;
}

double logsumexp(double *array, int len) {
   int i;
   double cur_max, sum;
   cur_max = -INF;
   for (i = 0; i < len; i++) {
      if (array[i] > cur_max) {
         cur_max = array[i];
      }
   }
   sum = 0;
   for (i = 0; i < len; i++) {
      sum += exp(array[i] - cur_max);
   }
   sum /= len;
   return log(sum) + cur_max;
}

void bcjr2(double *input, double sigma, int input_len, int *cdc_size, int *gen_mat, double *output, int *output_hard) {
   int i, j, k, word, *codeword, best_word;
   double *likelihood, Py, *apost, P0, P1, max_apost;
   codeword = (int *)malloc(sizeof(int) * cdc_size[1]);
   likelihood = (double *)malloc(sizeof(double) * pow(2, cdc_size[0]));
   apost = (double *)malloc(sizeof(double) * pow(2, cdc_size[0]));
   for (i = 0; i < pow(2, cdc_size[0]); i++) {
      likelihood[i] = 0;
   }

   word = 0, Py = 0, best_word = 0, max_apost = -INF;

   for (i = 0; i < pow(2, cdc_size[0]); i++) {
      memset(codeword, 0, sizeof(int) * cdc_size[1]);
      for (j = 0; j < cdc_size[0]; j++) {
         if (((word >> j) & 1)) {
            for (k = 0; k < cdc_size[1]; k++) {
               codeword[k] ^= gen_mat[j * cdc_size[1] + k];
            }
         }
      }
      word++;
      for (k = 0; k < cdc_size[1]; k++) {
         codeword[k] = 1 - 2 * codeword[k];
      }
      for (j = 0; j < cdc_size[1]; j++) {
         likelihood[i] += -log(sqrt(2*3.14)*sigma) - (pow(fabs(input[j] - codeword[j]), 2) / (2*sigma*sigma));
      }
   }
   Py = logsumexp(likelihood, pow(2, cdc_size[0]));

   for (i = 0; i < pow(2, cdc_size[0]); i++) {
      apost[i] = likelihood[i] - Py - log(pow(2, cdc_size[0]));
      if (apost[i] > max_apost) {
         max_apost = apost[i];
         best_word = i;
      }
   }
   for (i = 0; i < cdc_size[0]; i++) {
      word = 0;
      P0 = P1 = 0;
      for (j = 0; j < pow(2, cdc_size[0]); j++) {
         if (((word >> (i)) & 1)) {
            P1 += exp(apost[j]);
         }
         else {
            P0 += exp(apost[j]);
         }
         word++;
      }
      output[i] = log(P0 / P1);
   }
   memset(codeword, 0, sizeof(int) * cdc_size[1]);
   for (j = 0; j < cdc_size[0]; j++) {
      if (((best_word >> j) & 1)) {
         output_hard[j] = -1;
      }
      else {
         output_hard[j] = 1;
      }
   }
   free(codeword);
   free(apost);
   free(likelihood);
}

// cdc - вспомогательные данные 
// c_out - входящие llr
// x_dec - результат декодирования
// inner_matrix_inv - обратная матрица внутреннего кода (совпадает с матрицей внутреннего кода)
// in_cd_len - размер матрицы внутреннего кода
// out_cd_len - длина кода
// cdc_opts - способы декодирования. 0 - ML декодер, 1 - синдромное декоирование
// cdc_sizes - размер матриц внешних кодов
// grid4 - сетка для синдромного декодирования
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
)
{
   cdc_inst_type *dd;
   double *in_dec_input, *out_dec_input, *inner_codeword, **inner_matrix_cds, *inner_dec_codeword;
   int i, j, k, in_vecs_len, cur_pos_dec, *inner_matrix_size, *inner_matrix, *inner_dec_codeword_int, *ret_codeword;
   in_vecs_len = out_cd_len / in_cd_len;

   dd = (cdc_inst_type *)cdc;

   in_dec_input = (double *)malloc(sizeof(double) * out_cd_len);
   out_dec_input = (double *)malloc(sizeof(double) * out_cd_len);
   inner_codeword = (double *)malloc(sizeof(double) * in_cd_len);
   inner_dec_codeword = (double *)malloc(sizeof(double) * in_cd_len);
   inner_dec_codeword_int = (int *)malloc(sizeof(int) * in_cd_len);
   inner_matrix_cds = (double **)malloc(sizeof(double *) * in_vecs_len);
   inner_matrix_size = (int *)malloc(sizeof(int) * 2);
   inner_matrix = (int *)malloc(sizeof(int) * in_cd_len * in_cd_len);
   ret_codeword = (int *)malloc(sizeof(int) * 32);
   for (i = 0; i < in_vecs_len; i++) {
      inner_matrix_cds[i] = (double *)malloc(sizeof(double) * in_vecs_len);
   }
   
   // in_dec_input - вход во внутренний декодер
   for (i = 0; i < dd->c_n; i++) in_dec_input[i] = c_out[i] * dd->sg22;

   // inner_matrix_cds - кодовое слово
   for (i = 0; i < in_cd_len; i++) {
      for (j = 0; j < in_vecs_len; j++) {
         inner_matrix_cds[j][i] = in_dec_input[i*in_vecs_len + j];
      }
   }

   cur_pos_dec = 0;
   for (i = 0; i < in_cd_len; i++) {
      // выкинем нужное количество верхних строк матрицы внутреннего кода
      for (j = i; j < in_cd_len; j++) {
         for (k = 0; k < in_cd_len; k++) {
            inner_matrix[(j-i) * in_cd_len + k] = inner_matrix_inv[j * in_cd_len + k];
         }
      }
      inner_matrix_size[0] = in_cd_len - i;
      inner_matrix_size[1] = in_cd_len;
      for (j = 0; j < in_vecs_len; j++) {
         // inner_codeword - кодовое слово внутреннего кода
         for (k = 0; k < in_cd_len; k++) {
            inner_codeword[k] = inner_matrix_cds[j][k];
            inner_dec_codeword[k] = INF;
         }
         for (k = 0; k < in_cd_len; k++) {
            inner_codeword[k] /= dd->sg22;
         }
         bcjr2(inner_codeword, sqrt(2.0 / dd->sg22), in_cd_len, inner_matrix_size, inner_matrix, inner_dec_codeword, inner_dec_codeword_int);
         out_dec_input[i*in_vecs_len + j] = inner_dec_codeword[0];
         if (inner_dec_codeword_int[0] == -1) {
            for (k = 0; k < in_cd_len; k++) {
               if (inner_matrix[k]) {
                  inner_matrix_cds[j][k] = -inner_matrix_cds[j][k];
               }
            }
         }
      }
      if (cdc_opts[i] == 0) {
         printf("%d\n", cur_pos_dec);
         MLDec(out_dec_input + i * in_vecs_len, in_vecs_len, dd->dc.node_table + i*in_vecs_len,
               cdc_sizes + i*2, cdc_mats[i], x_dec + cur_pos_dec, ret_codeword);
         cur_pos_dec += cdc_sizes[2*i];
         printf("check\n");
      }
      else if (cdc_opts[i] == 1) {
         SyndrDec(out_dec_input + i * in_vecs_len, in_vecs_len, dd->dc.node_table + i*in_vecs_len,
                  cdc_sizes + i*2, cdc_mats[i], x_dec + cur_pos_dec, grid4);
      }
      // исправление ошибок с учетом декодирования внешнего кода
      for (j = 0; j < in_vecs_len; j++) {
         if ((1 - 2 * ret_codeword[j]) * out_dec_input[i * in_vecs_len + j] < 0) {
            for (k = 0; k < in_cd_len; k++) {
               if (inner_matrix[k]) {
                  inner_matrix_cds[j][k] = -inner_matrix_cds[j][k];
               }
            }
         }
      }
   }

   free(in_dec_input);
   free(out_dec_input);
   free(inner_codeword);
   free(inner_dec_codeword);
   free(inner_dec_codeword_int);
   free(inner_matrix_size);
   free(inner_matrix);
   for (i = 0; i < in_vecs_len; i++) {
      free(inner_matrix_cds[i]);
   }
   free(inner_matrix_cds);
   free(ret_codeword);
   return 0;
}
#endif