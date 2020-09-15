//=============================================================================
// Polar SCL.
// Can return the list or the best candidate.
// Internal format can be changed.
// Derived from dtrm_glp
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
// Configuration switches.

//#define DBG

//-----------------------------------------------------------------------------
// Includes.

#ifdef DBG
#include <stdio.h>
#endif


#include <stdio.h> // DELETE, ONLY FOR DEBUG
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../common/typedefs.h"
#include "../common//std_defs.h"
#include "polar_scl_inner.h"


//-----------------------------------------------------------------------------
// Internal defines.

#define NEAR_ZERO (1e-300)
#define INF (1e10)

// List access defines.
#define GETX(i) (lind2xl[i]->x)
#define PUSHX(v, i) { next_xle_ptr->x = v; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }
#define PUSHZEROX(i) { next_xle_ptr->x = 0; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }
#define PUSHONEX(i) { next_xle_ptr->x = 1; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }

#ifdef DBG2
#define YLIST_POP(i) (printf("pop i: %d, yitems[i]: %p, yfrindp[i]: %d\n", i, yitems[i], yfrindp[i]), yitems[i] + (1 << (i)) * (yfrindp[i]++));
#else // DBG2
#define YLIST_POP(i) (yitems[i] + (1 << (i)) * (yfrindp[i]++))
#endif // DBG2
#define YLIST_RESET(i) { yfrindp[i] = 0; }

#define YLISTP(i, j) ylist[(i) * dd->c_m + (j)]
#define YLISTPP(i, j) (ylist + (i) * dd->c_m + (j))


//-----------------------------------------------------------------------------
// Internal typedefs.

//-----------------------------------------------------------------------------
// Global data.

uint32 *node_table; // The table of list sizes at the border nodes.
int node_counter;

slitem *slist;
xlist_item *xlist;
int *plist; // Permutation index list.
int *parent; // index of the parent list element.
xlist_item **lind2xl; // list index to xlist

ylitem **yitems; // peak_lsiz * 2 * c_n
ylitem **ylist; // current buffer - YLISTP(list_index, m)
int yfrindp[32]; // yfrindp[i] points to the next available row for ylist in yitems[i].

int *lorder; // list ordering
int cur_lsiz; // Current list size.

int *frind; // Stack containing indexes of free list cells.
int frindp; // frind pointer. Points to the next available element.

xlist_item *next_xle_ptr; // pointer to the next available xlist element

xlitem *xtmp;

#ifdef FLIPPING
lv *lv_array_cur;
int fbit;
#endif // FLIPPING

#ifdef LISTFLIPPING
double *Malpha_cur;
int fstep;
slitem *slist_alpha;
double alpha_cur;
#endif // LISTFLIPPING

#if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
int fstep;
#endif

#if defined LISTFLIPPINGFAST
int is_enough;
double min_abs_llr;
int inv_subcodes_counter;
int *subcodes_first_symb;
double *subcodes_min_abs_llr;
#endif

//-----------------------------------------------------------------------------
// Functions.

void vgetx(int listInd, xlitem *res_x, int n) {
  int i;
  xlist_item *xle = lind2xl[listInd];
  for (i = 0; i < n; i++) {
    res_x[n - i - 1] = xle->x;
    xle = xle->p;
  }
}

#ifdef LISTFLIPPINGFAST
void polar_transform_mat(int *word, int n) {
  if (n == 1) return;
  int *lpart, *rpart;
  int i, n2;
  n2 = n / 2;
  lpart = (int *)malloc(sizeof(int) * n2);
  rpart = (int *)malloc(sizeof(int) * n2);
  for (i = 0; i < n2; i++) {
    lpart[i] = word[i] ^ word[i + n2];
    rpart[i] = word[i + n2];
  }
  polar_transform_mat(lpart, n2);
  polar_transform_mat(rpart, n2);
  for (i = 0; i < n2; i++) {
    word[i] = lpart[i];
  }
  for (i = n2; i < n; i++) {
    word[i] = rpart[i-n2];
  }
  free(lpart);
  free(rpart);
}

void polar_transf_x(int listInd, xlitem *res_x, int n) {
  int i;
  xlist_item *xle = lind2xl[listInd];
  polar_transform_mat(res_x, n);
  for (i = 0; i < n; i++) {
    xle->x = res_x[n - i - 1];
    xle = xle->p;
  }
}

void polar_transf_x_spc(int listInd, xlitem *res_x, int n) {
  int i;
  polar_transform_mat(res_x, n);
  lind2xl[listInd] = lind2xl[listInd]->p;
  xlist_item *xle = lind2xl[listInd];
  for (i = 0; i < n - 1; i++) {
    xle->x = res_x[n - i - 1];
    xle = xle->p;
  }
}
#endif

int branch(int i0, xlist_item *xle0, int ret_list) {
  int i1 = frind[frindp++];
  parent[i1] = i0;
  lind2xl[i1] = xle0;
  lorder[cur_lsiz++] = i1;
  slist[i1] = slist[i0];
  plist[i1] = plist[i0];
#ifdef LISTFLIPPING
  slist_alpha[i1] = slist_alpha[i0];
#endif // LISTFLIPPING
#ifdef FLIPPING 
  memcpy(lv_array_cur[i1].llrs, lv_array_cur[i0].llrs, node_counter * sizeof(double));
#endif // FLIPPING
  return i1;
}

#ifdef DBG
int vgetx_all(int listInd, xlitem *res_x) {
  int i, n;
  xlitem x;
  xlist_item *xle = lind2xl[listInd];

  i = 0;
  while (xle->p != NULL) {
    res_x[i++] = xle->x;
    xle = xle->p;
  }
  n = i;
  for (i = 0; i < n / 2; i++) {
    x = res_x[i];
    res_x[i] = res_x[n - i - 1];
    res_x[n - i - 1] = x;
  }
  return n;
}
void print_list(char rem[]) {
  int i, j, n;
  xlitem tmp[4096];
  if (rem != NULL) {
    printf("%s: ", rem);
  }
  for (i = 0; i < cur_lsiz; i++) {
    printf("%3.1f:", slist[lorder[i]]);
    n = vgetx_all(lorder[i], tmp);
    for (j = 0; j < n; j++) {
      printf("%1d", tmp[j]);
    }
    printf(" ");
  }
  printf("\n");
#ifdef DBG2
  printf("Y buffers: ");
  for (i = 0; i < 4; i++) {
    printf("%d:%p:%d, ", i, yitems[i], yfrindp[i]);
  }
  printf("\n");
#endif // DBG2
}
#endif // DBG

#define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }

void qpartition(int *v, int len, int k) {
  int i, st, tmp;

  for (st = i = 0; i < len - 1; i++) {
    if (slist[v[i]] < slist[v[len - 1]]) continue;
    SWAP(i, st);
    st++;
  }

  SWAP(len - 1, st);

  if (k == st) return;
  if (st > k) qpartition(v, st, k);
  else qpartition(v + st, len - st, k - st);
}

// For qsort().
static int lst_comp(const void *i, const void *j)
{
   if (slist[*(int *)i] < slist[*(int *)j]) return 1;
   return -1;
}

void polar0_branch(
   decoder_type *dd, // Decoder instance data.
   int peak_lsiz
) {
  int cur_lsiz_old = cur_lsiz;
  int cur_ind0, cur_ind1;
  ylitem y;
  xlitem x;
  ylitem *yp;
  int i;

#ifdef DBG
  printf("rm00 branch.\n");
  print_list("a");
#endif

  for (i = 0; i < cur_lsiz_old; i++) {

    cur_ind0 = lorder[i];
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0], dd->ret_list);

    yp = YLISTP(cur_ind0, 0);
    if (CHECK_EST0(yp[0])) {
      y = yp[0];
      x = 0;
    }
    else {
      y = INV_EST(yp[0]);
      x = 1;
    }

    PUSHX(x, cur_ind0);

    #if !defined LISTFLIPPINGFAST
    slist[cur_ind0] += EST0_TO_LNP0(y);
    #endif

    PUSHX(1 - x, cur_ind1);

    #if defined LISTFLIPPINGFAST
    slist[cur_ind1] -= y;
    #endif

    #if !defined LISTFLIPPINGFAST
    slist[cur_ind1] += EST0_TO_LNP1(y);
    #endif

    #ifdef LISTFLIPPING
    slist_alpha[cur_ind0] += EST0_TO_LNP0(alpha_cur * y);
    slist_alpha[cur_ind1] += EST0_TO_LNP1(alpha_cur * y);
    #endif // LISTFLIPPING

    
#ifdef FLIPPING
    // TODO: DIFFERENT FORMATS
    lv_array_cur[cur_ind0].llrs[node_counter] = -yp[0];
    lv_array_cur[cur_ind1].llrs[node_counter] = -yp[0];
#endif // FLIPPING
  }

  #ifdef LISTFLIPPING
  Malpha_cur[node_counter] = -INF;
  #endif // LISTFLIPPING

#ifdef DBG
  print_list("b");
#endif
  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) { 
    #ifdef LISTFLIPPING
    qpartition(lorder, cur_lsiz, peak_lsiz);
    ALPHA_CALC(slist_alpha, lorder, peak_lsiz, cur_lsiz, Malpha_cur, node_counter);
    if (node_counter != fstep) {
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    }
    else {
      for (i = 0; i < peak_lsiz; ++i) frind[--frindp] = lorder[i];
      for (i = peak_lsiz; i < cur_lsiz; ++i) lorder[i - peak_lsiz] = lorder[i];
      cur_lsiz = peak_lsiz;
    }
    #endif // LISTFLIPPING

    #if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
    qpartition(lorder, cur_lsiz, peak_lsiz);
    if (node_counter != fstep) {
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    }
    else {
      for (i = 0; i < peak_lsiz; ++i) frind[--frindp] = lorder[i];
      for (i = peak_lsiz; i < cur_lsiz; ++i) lorder[i - peak_lsiz] = lorder[i];
      cur_lsiz = peak_lsiz;
    }
    #endif

    #if !defined LISTFLIPPING && !defined LISTFLIPPINGPRECALC && !defined LISTFLIPPINGFAST
    qpartition(lorder, cur_lsiz, peak_lsiz);
    while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    #endif // LISTFLIPPING
  }

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, 0), YLISTPP(cur_ind0, 0), dd->c_m * sizeof(ylitem *));
    }
    x = GETX(cur_ind1);
    yp = YLISTP(cur_ind1, 0) = YLIST_POP(0);
    yp[0] = x ? YLDEC1 : YLDEC0;
#ifdef FLIPPING
    if (fbit == node_counter) {
      if (yp[0] == YLDEC1)
        yp[0] = YLDEC0;
      else
        yp[0] = YLDEC1;
      GETX(cur_ind1) = 1 - x;
    }
#endif // FLIPPING
  }

#ifdef DBG
  print_list("c");
#endif

}

void polar0_skip(
  decoder_type *dd // Decoder instance data.
) {
  int cur_ind0;
  ylitem *yp;
  int i;

#ifdef DBG
  printf("polar0 skip.\n");
  print_list("a");
#endif

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    yp = YLISTP(cur_ind0, 0);
    
    #if defined LISTFLIPPINGFAST
    slist[cur_ind0] += (yp[0] < 0) * yp[0];
    #endif

    #if !defined LISTFLIPPINGFAST
    slist[cur_ind0] += EST_TO_LNP0(yp[0]);
    #endif
    
#ifdef FLIPPING
    lv_array_cur[cur_ind0].llrs[node_counter] = -yp[0];
#endif // FLIPPING

#ifdef LISTFLIPPING
    Malpha_cur[node_counter] = -INF;
#endif // LISTFLIPPING
    yp = YLISTP(cur_ind0, 0) = YLIST_POP(0);
    yp[0] = YLDEC0;
  }

#ifdef DBG
  print_list("b");
#endif
}

#ifdef LISTFLIPPINGFAST

#define L_extra 4

void DecRate0(decoder_type *dd, int m, int n) {
  int i, j, cur_ind0;
  double sum_rate;
  ylitem *yp1;

  sum_rate = 0;
  
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    yp1 = YLISTP(cur_ind0, m);
    sum_rate = 0;
    for (j = 0; j < n; j++) {
      if (yp1[j] < 0) sum_rate += yp1[j];
    }
    slist[cur_ind0] += sum_rate;
    for (j = 0; j < n; j++) {
      yp1[j] = YLDEC0;
    }
  }
}

void DecRate1(decoder_type *dd, int m, int n, int peak_lsiz) {
  int cur_lsiz_old = cur_lsiz;
  int i, j, k, cur_ind0, cur_ind1, cur_ind2, cur_ind3, min_ind1, min_ind2;
  ylitem *yp1, *yp;
  int **words;
  xlitem *x_tmp;
  double logP, min1, min2, max_slist;

  words = malloc(sizeof(int *) * L_extra);
  for (i = 0; i < L_extra; i++) {
    words[i] = (int *)malloc(sizeof(int) * n);
  }
  x_tmp = (xlitem *)malloc(sizeof(xlitem) * n);

  max_slist = -INF;
  for (i = 0; i < cur_lsiz_old; i++) {
    cur_ind0 = lorder[i];
    yp1 = YLISTP(cur_ind0, m);
    logP = 0;
    for (j = 0; j < n; j++) {
      logP += -log(1 + exp(-fabs(yp1[j])));
    }
    slist[cur_ind0] += logP;
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0], dd->ret_list);
    cur_ind2 = branch(cur_ind0, lind2xl[cur_ind0], dd->ret_list);
    cur_ind3 = branch(cur_ind0, lind2xl[cur_ind0], dd->ret_list);
    for (j = 0; j < n; j++) {
      for (k = 0; k < L_extra; k++) {
        if (yp1[j] < 0) words[k][j] = 1;
        else words[k][j] = 0;
      }
    }
    min_ind1 = 0, min_ind2 = 1;
    min1 = fabs(yp1[0]), min2 = fabs(yp1[1]);
    if (min1 > min2) {
      min_ind1 = 1, min_ind2 = 0;
      min1 = fabs(yp1[1]), min2 = fabs(yp1[0]);
    }
    for (j = 2; j < n; j++) {
      if (fabs(yp1[j]) < min1) {
        min2 = min1;
        min_ind2 = min_ind1;
        min1 = fabs(yp1[j]);
        min_ind1 = j;
      }
      else if (fabs(yp1[j]) < min2) {
        min2 = fabs(yp1[j]);
        min_ind2 = j;
      }
    }
    words[1][min_ind1] ^= 1;
    words[3][min_ind1] ^= 1;
    words[2][min_ind2] ^= 1;
    words[3][min_ind2] ^= 1;
    slist[cur_ind1] -= min1;
    slist[cur_ind2] -= min2;
    slist[cur_ind3] -= (min1 + min2);
    for (k = 0; k < n; k++) {
      PUSHX(words[0][k], cur_ind0);
    }
    for (k = 0; k < n; k++) {
      PUSHX(words[1][k], cur_ind1);
    }
    for (k = 0; k < n; k++) {
      PUSHX(words[2][k], cur_ind2);
    }
    for (k = 0; k < n; k++) {
      PUSHX(words[3][k], cur_ind3);
    }
    if (max_slist < slist[cur_ind0]) {
      min_abs_llr = INF;
      for (j = 0; j < n; j++) {
        min_abs_llr = fmin(min_abs_llr, fabs(yp1[j]));
      }
      max_slist = slist[cur_ind0];
    }
  }

  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) {
    is_enough = 1; 

    qpartition(lorder, cur_lsiz, peak_lsiz);
    if (node_counter > fstep || node_counter + n <= fstep) {
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    }
    else {
      for (i = 0; i < peak_lsiz; ++i) frind[--frindp] = lorder[i];
      for (i = peak_lsiz; i < cur_lsiz; ++i) lorder[i - peak_lsiz] = lorder[i];
      cur_lsiz -= peak_lsiz;
      if (cur_lsiz > peak_lsiz) qpartition(lorder, cur_lsiz, peak_lsiz);
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
      cur_lsiz = peak_lsiz;
    }
  }

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, m), YLISTPP(cur_ind0, m), dd->c_m * sizeof(ylitem *));
    }
    vgetx(cur_ind1, x_tmp, n);
    yp = YLISTP(cur_ind1, m) = YLIST_POP(m);
    for (j = 0; j < n; j++) {
      if (x_tmp[j]) yp[j] = YLDEC1;
      else yp[j] = YLDEC0;
    }
    polar_transf_x(cur_ind1, x_tmp, n);
  }
  for (i = 0; i < L_extra; i++) {
    free(words[i]);
  }
  free(words);
  free(x_tmp);
}

void DecRepCode(decoder_type *dd, int m, int n, int peak_lsiz) {
  int cur_lsiz_old = cur_lsiz;
  int i, j, cur_ind0, cur_ind1;
  ylitem *yp1, *yp;
  xlitem *x_tmp;
  double sum_rate, max_slist;

  x_tmp = (xlitem *)malloc(sizeof(xlitem));

  max_slist = -INF;
  for (i = 0; i < cur_lsiz_old; i++) {
    cur_ind0 = lorder[i];
    yp1 = YLISTP(cur_ind0, m);
    sum_rate = 0;
    for (j = 0; j < n; j++) {
      if (yp1[j] < 0) sum_rate += yp1[j];
    }
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0], dd->ret_list);
    slist[cur_ind0] += sum_rate;
    sum_rate = 0;
    for (j = 0; j < n; j++) {
      if (yp1[j] > 0) sum_rate -= yp1[j];
    }
    slist[cur_ind1] += sum_rate;
    PUSHX(0, cur_ind0);
    PUSHX(1, cur_ind1);
    if (max_slist < slist[cur_ind0]) {
      min_abs_llr = INF;
      for (j = 0; j < n; j++) {
        min_abs_llr = fmin(min_abs_llr, fabs(yp1[j]));
      }
      max_slist = slist[cur_ind0];
    }
    if (max_slist < slist[cur_ind1]) {
      min_abs_llr = INF;
      for (j = 0; j < n; j++) {
        min_abs_llr = fmin(min_abs_llr, fabs(yp1[j]));
      }
      max_slist = slist[cur_ind1];
    }
  }

  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) { 
    is_enough = 1;

    qpartition(lorder, cur_lsiz, peak_lsiz);
    if (node_counter > fstep || node_counter + n <= fstep) {
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    }
    else {
      for (i = 0; i < peak_lsiz; ++i) frind[--frindp] = lorder[i];
      for (i = peak_lsiz; i < cur_lsiz; ++i) lorder[i - peak_lsiz] = lorder[i];
      cur_lsiz -= peak_lsiz;
    }
  }
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, m), YLISTPP(cur_ind0, m), dd->c_m * sizeof(ylitem *));
    }
    vgetx(cur_ind1, x_tmp, 1);
    yp = YLISTP(cur_ind1, m) = YLIST_POP(m);
    for (j = 0; j < n; j++) {
      if (x_tmp[0]) yp[j] = YLDEC1;
      else yp[j] = YLDEC0;
    }
  }
  free(x_tmp);
}

void DecSPCCode(decoder_type *dd, int m, int n, int peak_lsiz) {
  int cur_lsiz_old = cur_lsiz;
  int *cur_inds, *min_inds;
  int i, j, k, cur_ind0, cur_ind1, main_inv;
  ylitem *yp1, *yp;
  int **words;
  xlitem *x_tmp;
  double logP, *mins, max_slist;

  cur_inds = (int *)malloc(sizeof(int) * 8);
  min_inds = (int *)malloc(sizeof(int) * 4);
  mins = (double *)malloc(sizeof(double) * 4);
  words = malloc(sizeof(int *) * 8);
  for (i = 0; i < 8; i++) {
    words[i] = (int *)malloc(sizeof(int) * n);
  }
  x_tmp = (xlitem *)malloc(sizeof(xlitem) * n);

  max_slist = -INF;
  for (i = 0; i < cur_lsiz_old; i++) {
    cur_inds[0] = lorder[i];
    yp1 = YLISTP(cur_inds[0], m);

    logP = 0;
    for (j = 0; j < n; j++) {
      logP += -log(1 + exp(-fabs(yp1[j])));
    }
    slist[cur_inds[0]] += logP;
    parent[cur_inds[0]] = -1;
    for (j = 1; j < 8; j++) cur_inds[j] = branch(cur_inds[0], lind2xl[cur_inds[0]], dd->ret_list);
    main_inv = 0;
    for (j = 0; j < n; j++) {
      for (k = 0; k < 8; k++) {
        if (yp1[j] < 0) {
          words[k][j] = 1;
        }
        else words[k][j] = 0;
      }
      if (yp1[j] < 0) main_inv ^= 1;
    }
    for (j = 0; j < 4; j++) mins[j] = INF;
    for (j = 0; j < n; j++) {
      if (fabs(yp1[j]) < mins[0]) {
        mins[3] = mins[2];
        min_inds[3] = min_inds[2];
        mins[2] = mins[1];
        min_inds[2] = min_inds[1];
        mins[1] = mins[0];
        min_inds[1] = min_inds[0];
        mins[0] = fabs(yp1[j]);
        min_inds[0] = j;
      }
      else if (fabs(yp1[j]) < mins[1]) {
        mins[3] = mins[2];
        min_inds[3] = min_inds[2];
        mins[2] = mins[1];
        min_inds[2] = min_inds[1];
        mins[1] = fabs(yp1[j]);
        min_inds[1] = j;
      }
      else if (fabs(yp1[j]) < mins[2]) {
        mins[3] = mins[2];
        min_inds[3] = min_inds[2];
        mins[2] = fabs(yp1[j]);
        min_inds[2] = j;
      }
      else if (fabs(yp1[j]) < mins[3]) {
        mins[3] = fabs(yp1[j]);
        min_inds[3] = j;
      }
    }
    if (main_inv) {
      for (j = 0; j < 8; j++) {
        words[j][min_inds[0]] ^= 1;
        slist[cur_inds[j]] -= mins[0];
      }
      mins[0] = -mins[0];
    }
    // (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
    words[1][min_inds[0]] ^= 1; words[1][min_inds[1]] ^= 1;
    slist[cur_inds[1]] -= (mins[0] + mins[1]);

    words[2][min_inds[0]] ^= 1; words[2][min_inds[2]] ^= 1;
    slist[cur_inds[2]] -= (mins[0] + mins[2]);

    words[3][min_inds[0]] ^= 1; words[3][min_inds[3]] ^= 1;
    slist[cur_inds[3]] -= (mins[0] + mins[3]);

    words[4][min_inds[1]] ^= 1; words[4][min_inds[2]] ^= 1;
    slist[cur_inds[4]] -= (mins[1] + mins[2]);

    words[5][min_inds[1]] ^= 1; words[5][min_inds[3]] ^= 1;
    slist[cur_inds[5]] -= (mins[1] + mins[3]);

    words[6][min_inds[2]] ^= 1; words[6][min_inds[3]] ^= 1;
    slist[cur_inds[6]] -= (mins[2] + mins[3]);

    words[7][min_inds[0]] ^= 1; words[7][min_inds[1]] ^= 1;
    words[7][min_inds[2]] ^= 1; words[7][min_inds[3]] ^= 1;
    slist[cur_inds[7]] -= (mins[0] + mins[1] + mins[2] + mins[3]);

    for (j = 0; j < 8; j++) {
      for (k = 0; k < n; k++) {
        PUSHX(words[j][k], cur_inds[j]);
      }
    }
    if (max_slist < slist[cur_inds[0]]) {
      min_abs_llr = INF;
      for (j = 0; j < n; j++) {
        min_abs_llr = fmin(min_abs_llr, fabs(yp1[j]));
      }
      max_slist = slist[cur_inds[0]];
    }
  }

  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) { 
    is_enough = 1;

    qpartition(lorder, cur_lsiz, peak_lsiz);
    if (node_counter > fstep || node_counter + n <= fstep) {
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
    }
    else {
      for (i = 0; i < peak_lsiz; ++i) frind[--frindp] = lorder[i];
      for (i = peak_lsiz; i < cur_lsiz; ++i) lorder[i - peak_lsiz] = lorder[i];
      cur_lsiz -= peak_lsiz;
      if (cur_lsiz > peak_lsiz) qpartition(lorder, cur_lsiz, peak_lsiz);
      while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
      cur_lsiz = peak_lsiz;
    }
  }

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, m), YLISTPP(cur_ind0, m), dd->c_m * sizeof(ylitem *));
    }
    vgetx(cur_ind1, x_tmp, n);
    yp = YLISTP(cur_ind1, m) = YLIST_POP(m);
    for (j = 0; j < n; j++) {
      if (x_tmp[j]) yp[j] = YLDEC1;
      else yp[j] = YLDEC0;
    }
    polar_transf_x_spc(cur_ind1, x_tmp, n);
    vgetx(cur_ind1, x_tmp, n-1);
  }

  free(cur_inds);
  free(min_inds);
  free(mins);
  for (i = 0; i < 8; i++) {
    free(words[i]);
  }
  free(words);
  free(x_tmp);
}

#endif 

void polar_dec_inner(
  decoder_type *dd, // Decoder instance data.
  int m
) {

  int n = 1 << m;
  int n2 = n / 2;
  int cur_ind0;
  ylitem *yp1, *yp2, *vp, *up;
  ylitem y1;
  int i, j;
  #ifdef LISTFLIPPINGFAST
  int is_rate_0, is_rate_1, is_rep_code, is_spc_code;
  double sum_rate;
  #endif 

  if (m == 0) {
    if (node_table[node_counter] == 0) {
      polar0_skip(dd);
    }
    else {
      polar0_branch(dd, node_table[node_counter]);
    }
    node_counter++;
    return;
  }

  #ifdef LISTFLIPPINGFAST
  is_rate_0 = 1;
  is_rep_code = 1;
  for (i = node_counter; i < node_counter + n - 1; i++) {
    if (node_table[i] != 0) {
      is_rate_0 = 0;
      is_rep_code = 0;
    }
  }
  if (node_table[node_counter + n - 1] == 0) {
    is_rep_code = 0;
  }
  else {
    is_rate_0 = 0;
  }
  is_rate_1 = 1;
  is_spc_code = 1;
  if (node_table[node_counter] == 0) {
    is_rate_1 = 0;
  }
  else {
    is_spc_code = 0;
  }
  for (i = node_counter + 1; i < node_counter + n; i++) {
    if (node_table[i] == 0) {
      is_rate_1 = 0;
      is_spc_code = 0;
    }
  }
  if (is_rate_0) {
    DecRate0(dd, m, n);
    node_counter += n;
    return;
  }
  if (is_rate_1) {
    DecRate1(dd, m, n, node_table[node_counter]);
    if (is_enough) subcodes_first_symb[inv_subcodes_counter] = node_counter;
    else subcodes_first_symb[inv_subcodes_counter] = -1;
    subcodes_min_abs_llr[inv_subcodes_counter] = min_abs_llr;
    node_counter += n;
    inv_subcodes_counter++;
    return;
  }
  if (is_rep_code) {
    DecRepCode(dd, m, n, node_table[node_counter + n - 1]);
    if (is_enough) subcodes_first_symb[inv_subcodes_counter] = node_counter + n - 1;
    else subcodes_first_symb[inv_subcodes_counter] = -1;
    subcodes_min_abs_llr[inv_subcodes_counter] = min_abs_llr;
    node_counter += n;
    inv_subcodes_counter++;
    return;
  }
  if (is_spc_code) {
    DecSPCCode(dd, m, n, node_table[node_counter + n - 1]);
    if (is_enough) subcodes_first_symb[inv_subcodes_counter] = node_counter + 1;
    else subcodes_first_symb[inv_subcodes_counter] = -1;
    subcodes_min_abs_llr[inv_subcodes_counter] = min_abs_llr;
    node_counter += n;
    inv_subcodes_counter++;
    return;
  }
  #endif

#ifdef DBG2
  printf("innr a: m=%d, n=%3d.\n", m, n);
  print_list("innr a");
#endif

  // Calculate y_v = y_1 xor y_2.
  for (i = 0; i < cur_lsiz; i++) {
     cur_ind0 = lorder[i];
     yp1 = YLISTP(cur_ind0, m);
     yp2 = yp1 + n2;
     vp = YLISTP(cur_ind0, m - 1) = YLIST_POP(m - 1);
     VXOR_EST(yp1, yp2, vp, n2);
  }

  polar_dec_inner(dd, m - 1);

#ifdef DBG2
  printf("innr b: m=%d, n=%3d.\n", m, n);
  print_list("innr b");
#endif

  // Calculate y_u = y_1 xor v + y_2.
  // Also save ref to v in y.
  for (i = 0; i < cur_lsiz; i++) {
     cur_ind0 = lorder[i];
     yp1 = YLISTP(cur_ind0, m);
     yp2 = yp1 + n2;
     vp = YLISTP(cur_ind0, m) = YLISTP(cur_ind0, m - 1); // y <-- v
     up = YLISTP(cur_ind0, m - 1) = YLIST_POP(m - 1);
     for (j = 0; j < n2; j++) {
        y1 = EST_XOR_YLDEC(yp1[j], vp[j]); // y1 <-- y_1[j] xor v[j].
        ADD_EST(y1, yp2[j], up[j]); // y_u <-- y1 + y_2[j];
     }
  }

#ifdef DBG2
  printf("innr c: m=%d, n=%3d.\n", m, n);
#endif

  polar_dec_inner(dd, m - 1);

#ifdef DBG2
  printf("innr d: m=%d, n=%3d.\n", m, n);
#endif

  // y_dec <-- (u xor v | u).
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    vp = YLISTP(cur_ind0, m);
    yp1 = YLISTP(cur_ind0, m) = YLIST_POP(m);
    yp2 = yp1 + n2;
    up = YLISTP(cur_ind0, m - 1);
    for (j = 0; j < n2; j++) {
      yp1[j] = XOR_YLDEC(vp[j], up[j]);
      yp2[j] = up[j];
    }
  }

  YLIST_RESET(m - 1);
}

// Decoding procedure.
int
polar_dec(
   decoder_type *dd, // Decoder instance data.
   double *y_input, // Decoder input.
   int *x_dec, // Decoded information sequence.
   double *s_dec, // Metric of x_dec (set NULL, if not needed).
   lv *lv_array, // Comparator value + llrs.
   int Ti, // Index of the flipping bit.
   int Si, // Index of the flipping step.
   double alpha, // SCLFlip constant.
   double *Malpha, // List of alphas.
   int *sfs,
   double *smal
)
{
  int i, j, i1;
  int c_n = dd->c_n;
  int c_k = dd->c_k;
  int c_m = dd->c_m;
  int peak_lsiz; // Peak list size.
  int flsiz; // Size of allocated list.
  ylitem *y_in;
  int n2 = dd->c_n / 2; // Half of the current length (n).
  int cur_ind0;
  ylitem *vp, *up;
  int *p1; // For permutations.
  slitem s1;
  ylitem y1;
  int *xp;
  uint8 *mem_buf_ptr;
  int mem_buf_size;

  // ----- Shortcuts assignment.

  node_table = dd->node_table;
  peak_lsiz = dd->peak_lsiz;
  flsiz = MAX(dd->peak_lsiz, dd->p_num) * FLSIZ_MULT;

  // ----- Memory allocation.

  mem_buf_ptr = dd->mem_buf;

  y_in = (ylitem *)mem_buf_ptr;
  mem_buf_ptr += c_n * sizeof(ylitem);

  frind = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  lorder = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  xlist = (xlist_item *)mem_buf_ptr;
  mem_buf_ptr += (c_k * flsiz + 1) * sizeof(xlist_item);

  lind2xl = (xlist_item **)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(xlist_item *);

  ylist = (ylitem **)mem_buf_ptr;
  mem_buf_ptr += dd->c_m * flsiz * sizeof(ylitem *);

  yitems = (ylitem **)mem_buf_ptr;
  mem_buf_ptr += dd->c_m * sizeof(ylitem *);
#ifdef DBG2
  printf("\nStart buffers: %ld\n", mem_buf_ptr - dd->mem_buf);
#endif
  for (i = 0; i < dd->c_m; i++) {
    yitems[i] = (ylitem *)mem_buf_ptr;
    mem_buf_ptr += (1 << i) * MAX(dd->peak_lsiz, dd->p_num) * 4 * sizeof(ylitem);
#ifdef DBG2
    printf("Buffer %d: yitems[i]=%p, %ld\n", i, yitems[i], mem_buf_ptr - dd->mem_buf);
#endif
  }

  slist = (slitem *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(slitem);

  plist = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  parent = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  xtmp = (xlitem *)mem_buf_ptr;
  mem_buf_ptr += c_n * sizeof(xlitem);

#ifdef FLIPPING
  lv_array_cur = (lv *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(lv);
  for (i = 0; i < flsiz; i++) {
      lv_array_cur[i].llrs = (double *)mem_buf_ptr;
      mem_buf_ptr += c_n * sizeof(double);
  }
  fbit = Ti;
#endif // FLIPPING

#ifdef LISTFLIPPING
  Malpha_cur = (double *)mem_buf_ptr;
  mem_buf_ptr += c_n * sizeof(double);
  fstep = Si;
  alpha_cur = alpha;
  slist_alpha = (slitem *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(slitem);
#endif // LISTFLIPPING

#if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
  fstep = Si;
#endif

#if defined LISTFLIPPINGFAST
is_enough = 0;
inv_subcodes_counter = 0;
subcodes_first_symb = sfs;
subcodes_min_abs_llr = smal;
#endif

#ifdef DBG2
  printf("\nAssigned buffer size: %ld\n", mem_buf_ptr - dd->mem_buf);
#endif

  // ----- Initial assignments.

  VY_TO_FORMAT(y_input, y_in, c_n);

  for (i = 0; i < flsiz; i++) frind[i] = i;
  for (i = 0; i < 32; i++) YLIST_RESET(i);

  node_counter = 0;

  //---------- Root node.

  // Calculate y_v = y_1 xor y_2.
  // Branch with permutations at (c_m, c_m).
  // List size is going to be dd->p_num.
  // All path metrics are 0 so far.
  
  for (i = 0; i < dd->p_num; i++) {
    slist[i] = 0.0;
    #ifdef LISTFLIPPING
    slist_alpha[i] = 0.0;
    #endif // LISTFLIPPING
    lorder[i] = i;
    lind2xl[i] = xlist;
    p1 = dd->pyarr + c_n * i;
    vp = YLISTP(i, c_m - 1) = YLIST_POP(c_m - 1);
    // TODO: Permutations are not usable with general Polar, but let's leave it as it is for now.
    for (j = 0; j < n2; j++) {
      vp[j] = XOR_EST(y_in[p1[j]], y_in[p1[j + n2]]);
    }
    plist[i] = i;
  }

  cur_lsiz = frindp = dd->p_num;
  next_xle_ptr = xlist + 1;
  xlist[0].x = 0;
  xlist[0].p = NULL;

#ifdef DBG2
  printf("root 1: m=%d, n=%3d.\n", c_m, c_n);
  print_list("root a");
#endif

  // Decode v.
  polar_dec_inner(dd, c_m - 1);

#ifdef DBG2
  print_list("root b");
#endif

  // Calculate y_u = y_1 xor v + y_2.
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    vp = YLISTP(cur_ind0, c_m - 1);
    up = YLISTP(cur_ind0, c_m - 1) = YLIST_POP(c_m - 1);
    // TODO: Permutations.
    p1 = dd->pyarr + c_n * plist[cur_ind0];
    for (j = 0; j < n2; j++) {
      y1 = EST_XOR_YLDEC(y_in[p1[j]], vp[j]); // y1 <-- y_1[j] xor v[j].
      ADD_EST(y1, y_in[p1[j + n2]], up[j]); // y_u <-- y1 + y_2[j];
    }
  }

  // Decode u.
  polar_dec_inner(dd, c_m - 1);

#ifdef DBG
  print_list("final");
#endif

#ifdef LISTFLIPPINGFAST
subcodes_first_symb[inv_subcodes_counter] = -2;
#endif

  if (dd->ret_list) {

    // Sort.
    qsort(lorder, cur_lsiz, sizeof(int), lst_comp);
#ifdef DBG
    print_list("sorted");
#endif
    for (i = 0; i < cur_lsiz; i++) {
      vgetx(lorder[i], xtmp, c_k);
      // TODO: Permutations.
      p1 = dd->pxarr + dd->c_k * plist[lorder[i]];
      xp = x_dec + i * c_k;
      for (j = 0; j < c_k; j++) xp[p1[j]] = xtmp[j];
#ifdef FLIPPING
      // TODO: Permutations => update lv_array_cur.llrs
      memcpy(lv_array[i].llrs, lv_array_cur[lorder[i]].llrs, c_n * sizeof(double));
      lv_array[i].l = slist[lorder[i]];
#endif // FLIPPING
#ifdef LISTFLIPPING
      memcpy(Malpha, Malpha_cur, c_n * sizeof(double));
#endif // LISTFLIPPING
    }

  }
  else {

    // Find the best.
    i1 = 0;
    for (i = 1; i < cur_lsiz; i++) {
      if (slist[lorder[i]] > slist[lorder[i1]]) {
        i1 = i;
      }
    }

#ifdef DBG
    printf("Best: %3.1f", slist[lorder[i1]]);
#endif

    vgetx(lorder[i1], xtmp, c_k);
    // TODO: Permutations.
    p1 = dd->pxarr + dd->c_k * plist[lorder[i1]];
    for (j = 0; j < c_k; j++) x_dec[p1[j]] = xtmp[j];

    if (s_dec != NULL) {
      (*s_dec) = slist[lorder[i1]];
    }
  }
  return 0;
}
