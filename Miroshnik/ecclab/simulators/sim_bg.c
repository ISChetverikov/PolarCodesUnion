//=============================================================================
// Simulation in the binary Gaussian channel
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

#define SIM_BG_C

//-----------------------------------------------------------------------------
// Configuration switches.

// Use the stdlib rand() function as the uniform random numbers generator.
// #define USE_STDLIB_RND

//-----------------------------------------------------------------------------
// Includes.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/spf_par.h"
#include "../interfaces/simul.h"
#include "../interfaces/codec.h"
#include "../interfaces/ui_utils.h"

//-----------------------------------------------------------------------------
// Internal defines.

// Max # of SNR values.
#define SNR_NUM_MAX        50
// Max length of file names.
#define FN_LEN_MAX         1000

// Defines for pseudorandom numbers generator.
#ifdef USE_STDLIB_RND
#define RND ((double)rand() / RAND_MAX)
#else // USE_STDLIB_RND
// #define rnd_A     16807UL
// #define rnd_M     2147483647UL
#define rnd_A     16807.0
#define rnd_M     2147483647.0 // (2^31 - 1) - a prime number.
#define rnd_D     4.656612875e-10
// #define RND ((double)(rnd_seed = (rnd_A * rnd_seed) % rnd_M) * rnd_D)
#define RND ((rnd_seed = fmod(rnd_A * rnd_seed, rnd_M)) * rnd_D)
#endif // else USE_STDLIB_RND

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
   void *dc_inst; // Codec instance.
   char rf_name[FN_LEN_MAX]; // Results file name.
   int ret_int; // Return interval.
   int code_n;
   int code_k;
   int use_rndcw; // If 1 use random codeword.
   double fixedR; // If set, simulate as if the code has this rate (error rate vs 1/R * Es/N0).
   int snr_num;
   int csnrn;
   double snr_db[SNR_NUM_MAX];
   int trn_req[SNR_NUM_MAX];
   int trn[SNR_NUM_MAX];
   int trn_saved[SNR_NUM_MAX];
   int en_bit[SNR_NUM_MAX];
   int en_bit_saved[SNR_NUM_MAX];
   int en_bl[SNR_NUM_MAX];
   int en_bl_saved[SNR_NUM_MAX];
   int do_ml; // If 1 evaluate ML LB.
   int do_ml_hd; // If 1 evaluate ML LB for hard dec. decoder.
   int enml_bl[SNR_NUM_MAX];
   int enml_bl_saved[SNR_NUM_MAX];
   #ifdef LISTFLIPPING
   int flip_num;
   double flip_alpha;
   #endif // LISTFLIPPING
   #if defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
   int flip_num;
   #endif
} sim_bg_inst;

//-----------------------------------------------------------------------------
// Global data.

// For pseudorandom numbers generator.
// unsigned long rnd_seed = 1;
double rnd_seed = 1.0;

// Simulation parameters file keywords.
char res_file_token[] = "res_file";
char snr_val_trn_token[] = "SNR_val_trn";
char ml_lb_token[] = "ml_lb";
char random_codeword_token[] = "random_codeword";
char fixedR_token[] = "fixed_R";

//-----------------------------------------------------------------------------
// Functions.

static double db2val(double x) {
   return exp(log(10.0) * x / 10.0);
}

static void set_rnd_seed(long s) {
#ifdef USE_STDLIB_RND
   srand((unsigned)s);
#else
   rnd_seed = (double) s;
#endif
}

/* Pauses for a specified number of milliseconds. */
static void sleep( clock_t wait )
{
   clock_t goal;
   goal = wait + clock();
   while( goal > clock() )
      ;
}

// Copy and add Gaussian noise using Marsaglia polar method.
static void copy_add_noise(
   double const ys[],
   int n, // Length of y[] (must be even!).
   double sg, // Standard deviation (sigma).
   double yd[]
)
{
   double v1, v2, r;
   int i = 0;

   while (i < n) {
      do {
         v1 = 2.0 * RND - 1.0;
         v2 = 2.0 * RND - 1.0;
         r = v1 * v1 + v2 * v2;
      } while (r >= 1.0);
      r = sqrt((-2.0 * log(r)) / r);
      yd[i] = ys[i] + v1 * r * sg;
      i++;
      yd[i] = ys[i] + v2 * r * sg;
      i++;
   }
}

// Allocate and init a simulation instance.
int sim_init(
   sim_init_params *sp, // Simulation parameters.
   void **inst // Simulation instance.
)
{
   char *sp_str; // Simulation parameters string.
   sim_bg_inst *sim; // Simulator instance.
   // Temporary variables.
   char *token, *str1;

   // Read and preparse simulation parameters.
   if (spf_read_preparse(sp->spf_name, &sp_str)) return 1;

   // Allocate simulator instance.
   sim = (sim_bg_inst *)malloc(sizeof(sim_bg_inst));
   if (sim == NULL) {
      err_msg("sim_init(): short of memory!");
      return 1;
   }
   memset(sim, 0, sizeof(sim_bg_inst));

   // Allocate temporary string for parsing.
   str1 = (char *)malloc(SP_STR_MAX);
   if (str1 == NULL) {
      err_msg("sim_init(): short of memory!");
      return 1;
   }
   str1[0] = 0;

   // Parse simulation parameters.
   strcpy(str1, sp_str);
   token = strtok(str1, tk_seps_prepared);
   while (token != NULL) {

      if (strcmp(token, res_file_token) == 0) {
         token = strtok(NULL, tk_seps_prepared);
         if (strlen(token) >= FN_LEN_MAX) {
            err_msg("sim_init(): Output file name is too long.");
            return 1;
         }
         strcpy(sim->rf_name, token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, snr_val_trn_token) == 0) {
         int i = 0;
         token = strtok(NULL, tk_seps_prepared);
         if (token[0] == GROUP_CHAR_BEGIN) {
            token = strtok(NULL, tk_seps_prepared);
            while ((token[0] != GROUP_CHAR_END) && (i < SNR_NUM_MAX)) {
               sim->snr_db[i] = atof(token);
               token = strtok(NULL, tk_seps_prepared);
               sim->trn_req[i++] = atoi(token);
               token = strtok(NULL, tk_seps_prepared);
            }
         }
         else {
            err_msg("sim_init(): Error in parameters file: can't read SNR_val_trn.");
            return 1;
         }
         sim->snr_num = i;
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      TRYGET_ONOFF_TOKEN(token, ml_lb_token, sim->do_ml);
      TRYGET_ONOFF_TOKEN(token, random_codeword_token, sim->use_rndcw);
      TRYGET_ONOFF_TOKEN(token, "ml_lb_hard", sim->do_ml_hd);
      TRYGET_FLOAT_TOKEN(token, fixedR_token, sim->fixedR);

      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   // Init codec.
   strcpy(str1, sp_str);
   if (cdc_init(str1, &(sim->dc_inst))) {
      err_msg("Cannot init codec.");
      return 1;
   }

   sim->code_n = cdc_get_n(sim->dc_inst);
   sim->code_k = cdc_get_k(sim->dc_inst);
   #if defined LISTFLIPPING || defined LISTFLIPPINGPRECALC || defined LISTFLIPPINGFAST
   sim->flip_num = cdc_get_flip_num(sim->dc_inst);
   #endif
   #ifdef LISTFLIPPING
   sim->flip_alpha = cdc_get_flip_alpha(sim->dc_inst);
   #endif

   // Check if all necessary parameters are specified.
   if ((sim->code_n == 0) || (sim->code_k == 0)) {
      err_msg("sim_init: code_n or code_k are not specified.");
      return 1;
   }

   // Assign other settings and init counters in the simulation instance.
   sim->csnrn = 0;
   sim->ret_int = sp->ret_int;
   if (sp->dont_randomize == 0) set_rnd_seed((int)time(NULL));
   if (sim->do_ml_hd) sim->do_ml = 1;

   free(sp_str);
   free(str1);

   (*inst) = sim;

   return 0;
}

// Delete simulation instance.
void sim_close(
   void *inst // Simulation instance.
) {
   sim_bg_inst *sim;

   sim = (sim_bg_inst *)inst;
   cdc_close(sim->dc_inst);
   free(sim);
}

#if !defined GCCDEC
// Main simulation cycle routine.
// Return: 0 - completed, >0 - not completed, <0 - error.
int sim_run(
   void *inst // Simulation instance.
) 
{
   sim_bg_inst *sim;
   int c_n; // Code length.
   int c_k; // Code dimension.
   double c_R; // Code rate.
   int csnrn; // Current SNR value number.
   // int trn; // (Current) trials number.
   double noise_sg; // Current value of sigma.
   int *x; // Information vector.
   int *x_dec; // Decoded information vector.
   double *c_in; // Channel input.
   double *c_out; // Channel output.
   time_t start_time;
   int en = 0;
   int i;
   #ifdef LISTFLIPPINGFAST
   int *Malpha_maxs, freq;
   FILE* precalc_bits;
   #endif

   sim = (sim_bg_inst *)inst;
   
   c_n = sim->code_n;
   c_k = sim->code_k;
   c_R = sim->fixedR ? sim->fixedR : ((double)c_k / (double)c_n);
   csnrn = sim->csnrn; // Current SNR value number.

   time(&start_time); // Remember start time.

   // Allocate buffers.
   x = (int *)malloc(c_n * sizeof(int));
   x_dec = (int *)malloc(c_n * sizeof(int));
   c_in = (double *)malloc(c_n * sizeof(double));
   c_out = (double *)malloc(c_n * sizeof(double));

   #if defined LISTFLIPPINGTHRESHOLD
   int cur_th;
   int max_threshold = 15;
   int* cnt_errors_threshold = malloc(sizeof(int) * max_threshold);
   int* cnt_T_threshold = malloc(sizeof(int) * max_threshold);
   memset(cnt_errors_threshold, 0, sizeof(int) * max_threshold);
   memset(cnt_T_threshold, 0, sizeof(int) * max_threshold);
   #endif

   /*#ifdef LISTFLIPPINGFAST
   Malpha_maxs = (int *)malloc(sim->flip_num * sizeof(int));
   precalc_bits = fopen("./precalc_bits/precalc_polar7_ca6_k64_L8_T10.spf", "r");
   i = 0;
   while (!feof(precalc_bits)) {
      if (i == sim->flip_num) break;
      fscanf(precalc_bits, "%d %d", &Malpha_maxs[i], &freq);
      i++;
   }
   fclose(precalc_bits);
   #endif*/
   if ((x == NULL) || (x_dec == NULL) || (c_in == NULL) || (c_out == NULL)) {
      err_msg("sim_run(): Short of memory.");
      return -1;
   }

   // Main simulation loop.
   while (csnrn < sim->snr_num) {

      noise_sg = 1 / sqrt(2 * c_R * db2val(sim->snr_db[csnrn]));

#ifdef DEC_NEEDS_CSNRN
      // Set current SNR # in the codec.
      cdc_set_csnrn(sim->dc_inst, csnrn);
#endif // DEC_NEEDS_CSNRN

#ifdef DEC_NEEDS_SIGMA
      // Set sg in the codec.
      cdc_set_sg(sim->dc_inst, noise_sg);
#endif // DEC_NEEDS_SIGMA

      int *var1 = malloc(sizeof(int) * c_n), *var2 = malloc(sizeof(int) * c_n), *var3 = malloc(sizeof(int) * c_n);
      for (i = 0; i < c_n; i++) var1[i] = var2[i] = var3[i] = 0;

      while (sim->trn[csnrn] < sim->trn_req[csnrn]) {
         
         sim->trn[csnrn]++;

         // Generate random or zero code word.
         if (sim->use_rndcw) {
            for (i = 0; i < c_k; i++) x[i] = (RND > 0.5) ? 1 : 0;
            enc_bpsk(sim->dc_inst, x, c_in);
         }
         else {
            for (i = 0; i < c_k; i++) x[i] = 0;
            for (i = 0; i < c_n; i++) c_in[i] = -1.0;
         }

         // Generate channel output for the simulated codeword.
         // (Channel simulation.)
         // c_out <-- c_in + noise.
         copy_add_noise(c_in, c_n, noise_sg, c_out);

         // Decode.

         #ifdef FLIPPING
         if (dec_bpsk_flipping(sim->dc_inst, c_out, x_dec, 32)) {
            err_msg("sim_run(): Error while decoding.");
            return -1;
         }
         #endif // FLIPPING

         #ifdef LISTFLIPPING
         #if !defined LISTFLIPPINGOPT
            if (dec_bpsk_list_flipping(sim->dc_inst, c_out, x_dec, sim->flip_num, sim->flip_alpha, var1, var2, var3)) {
               err_msg("sim_run(): Error while decoding.");
               return -1;
            }
         #endif

         #if defined LISTFLIPPINGOPT
         
         #if !defined LISTFLIPPINGTHRESHOLD
            if (dec_bpsk_list_flipping_opt(sim->dc_inst, c_out, x_dec, sim->flip_num, sim->flip_alpha, var1, var2, var3, x)) {
               err_msg("sim_run(): Error while decoding.");
               return -1;
            }
         #endif

         #if defined LISTFLIPPINGTHRESHOLD
            for (cur_th = 0; cur_th < max_threshold; cur_th++) {
               var3[0] = 0;
               if (dec_bpsk_list_flipping_threshold(sim->dc_inst, c_out, x_dec, c_n, sim->flip_alpha, var1, var2, var3, x, (double)cur_th)) {
                  err_msg("sim_run(): Error while decoding.");
                  return -1;
               }
               en = 0;
               for (i = 0; i < c_k; i++) if (x_dec[i] != x[i]) en++;
               if (en) {
                  cnt_errors_threshold[cur_th]++;
                  if (var3[0]) var2[0]++;
               }
               cnt_T_threshold[cur_th] += var1[0];
            }
         #endif

         #endif
         #endif // LISTFLIPPING

         #ifdef LISTFLIPPINGPRECALC
         if (dec_bpsk_list_flipping_precalc(sim->dc_inst, c_out, x_dec, sim->flip_num)) {
            err_msg("sim_run(): Error while decoding.");
            return -1;
         }
         #endif // LISTFLIPPINGPRECALC

         #ifdef LISTFLIPPINGFAST
         if (dec_bpsk_list_flipping_fast(sim->dc_inst, c_out, x_dec, sim->flip_num)) {
            err_msg("sim_run(): Error while decoding.");
            return -1;
         }
         #endif // LISTFLIPPINGFAST

         #if !defined FLIPPING && !defined LISTFLIPPING && !defined LISTFLIPPINGPRECALC && !defined LISTFLIPPINGFAST
         if (dec_bpsk(sim->dc_inst, c_out, x_dec)) {
            err_msg("sim_run(): Error while decoding.");
            return -1;
         }
         #endif

         // Count the number of incorrect information bits.
         en = 0;
         for (i = 0; i < c_k; i++) if (x_dec[i] != x[i]) en++;
         // Adjust error counters.
         if (en) {
            sim->en_bit[csnrn] += en;
            sim->en_bl[csnrn]++;
         }

         // Check if the ML decoding would fail too.
         if (sim->do_ml) {
            double d1 = 0.0, d2 = 0.0;
            if (sim->do_ml_hd)
               for (i = 0; i < c_n; i++) c_out[i] = (c_out[i] > 0.0) ? 1.0 : -1.0;
            // d1 <-- dist(c_in, c_out).
            for (i = 0; i < c_n; i++) d1 += (c_in[i] - c_out[i]) * (c_in[i] - c_out[i]);
            // d2 <-- dist[encode(x_dec), c_out].
            enc_bpsk(sim->dc_inst, x_dec, c_in);
            for (i = 0; i < c_n; i++) d2 += (c_in[i] - c_out[i]) * (c_in[i] - c_out[i]);
            if (d1 > d2) sim->enml_bl[csnrn]++;
         }

         // Check elapsed time.
         if (time(NULL) - start_time > sim->ret_int) {
            free(c_out);
            free(c_in);
            free(x_dec);
            free(x);
            return 1;
         }
      }

      #if defined LISTFLIPPINGTHRESHOLD
      printf("Number of CRC errors: %d\n", var2[0]);
      printf("threshold FER T\n");
      for (cur_th = 0; cur_th < max_threshold; cur_th++) {
         printf("%d %f %f\n", cur_th, (double) cnt_errors_threshold[cur_th] / sim->trn_req[csnrn], (double) cnt_T_threshold[cur_th] / sim->trn_req[csnrn]);
      }
      #endif

      sim->csnrn = ++csnrn;

      // precalc for SCLFlip
      /*for (i = 0; i < c_n; i++) {
         printf("%d %d\n", i, var1[i]);
      }
      printf("\n");
      for (i = 0; i < c_n; i++) {
         printf("%d %d\n", i, var2[i]);
      }
      printf("\n");
      for (i = 0; i < c_n; i++) {
         printf("%d %d\n", i, var3[i]);
      }
      printf("\n");*/

   }

   free(c_out);
   free(c_in);
   free(x_dec);
   free(x);

   return 0;
}
#endif // !GCCDEC

#ifdef GCCDEC
int sim_run(
   void *inst // Simulation instance.
) 
{
   sim_bg_inst *sim;
   int c_n; // Code length.
   int c_k; // Code dimension.
   double c_R; // Code rate.
   int csnrn; // Current SNR value number.
   // int trn; // (Current) trials number.
   double noise_sg; // Current value of sigma.
   int *x; // Information vector.
   int *x_dec; // Decoded information vector.
   double *c_in; // Channel input.
   double *c_out; // Channel output.
   time_t start_time;
   int en = 0;
   int i, j, k, el, *outer_matrix, *inner_matrix, *inner_matrix_inv, in_cd_len, out_cd_len;
   int *cdc_opts, *cdc_sizes, **cdc_mats, ***grid4;
   int err_nums1, err_nums2, err_nums3, err_nums4;

   sim = (sim_bg_inst *)inst;
   
   c_n = sim->code_n;
   c_k = sim->code_k;
   c_R = sim->fixedR ? sim->fixedR : ((double)c_k / (double)c_n);
   csnrn = sim->csnrn; // Current SNR value number.
   in_cd_len = 4, out_cd_len = 128;

   time(&start_time); // Remember start time.

   // Allocate buffers.
   x = (int *)malloc(c_n * sizeof(int));
   x_dec = (int *)malloc(c_n * sizeof(int));
   c_in = (double *)malloc(c_n * sizeof(double));
   c_out = (double *)malloc(c_n * sizeof(double));
   if ((x == NULL) || (x_dec == NULL) || (c_in == NULL) || (c_out == NULL)) {
      err_msg("sim_run(): Short of memory.");
      return -1;
   }
   outer_matrix = (int *)malloc(sizeof(int) * out_cd_len * out_cd_len);
   inner_matrix = (int *)malloc(sizeof(int) * in_cd_len * in_cd_len);
   inner_matrix_inv = (int *)malloc(sizeof(int) * in_cd_len * in_cd_len);
   cdc_opts = (int *)malloc(sizeof(int) * in_cd_len);
   cdc_sizes = (int *)malloc(sizeof(int) * 2 * in_cd_len);
   grid4 = (int ***)malloc(sizeof(int**) * 4);
   for (i = 0; i < 4; i++) {
      grid4[i] = (int **)malloc(sizeof(int*) * 33);
      for (j = 0; j < 33; j++) {
         grid4[i][j] = (int *)malloc(sizeof(int) * 5);
         for (k = 0; k < 5; k++) {
            grid4[i][j][k] = 0;
         }
      }
   }

   cdc_sizes[0] = 2, cdc_sizes[1] = 32;
   cdc_sizes[2] = 15, cdc_sizes[3] = 32;
   cdc_sizes[4] = 17, cdc_sizes[5] = 32;
   cdc_sizes[6] = 2, cdc_sizes[7] = 32;


   cdc_mats = (int **)malloc(sizeof(int *) * in_cd_len);
   cdc_mats[0] = (int *)malloc(sizeof(int) * 2 * 32);
   cdc_mats[1] = (int *)malloc(sizeof(int) * 15 * 32);
   cdc_mats[2] = (int *)malloc(sizeof(int) * 17 * 32);
   cdc_mats[3] = (int *)malloc(sizeof(int) * 2 * 32);

   cdc_opts[0] = 0, cdc_opts[1] = 0, cdc_opts[2] = 0, cdc_opts[3] = 1;

   //TODO: перенести в init

   FILE* outer_code = fopen("./GCCmats/sylvester7.mat", "r");

   for (i = 0; i < out_cd_len; i++) {
      for (j = 0; j < out_cd_len; j++) {
         fscanf(outer_code, "%d", &el); 
         outer_matrix[i*out_cd_len + j] = el;
      }
   }
   fclose(outer_code);

   FILE* inner_code = fopen("./GCCmats/sylvester2.mat", "r");

   for (i = 0; i < in_cd_len; i++) {
      for (j = 0; j < in_cd_len; j++) {
         fscanf(inner_code, "%d", &el);
         inner_matrix[i * in_cd_len + j] = el;
      }
   }
   fclose(inner_code);

   FILE* inner_code_inv = fopen("./GCCmats/sylvester2.mat", "r");

   for (i = 0; i < in_cd_len; i++) {
      for (j = 0; j < in_cd_len; j++) {
         fscanf(inner_code_inv, "%d", &el);
         inner_matrix_inv[i * in_cd_len + j] = el;
      }
   }
   fclose(inner_code_inv);

   // следующие 4 матрицы - порождающие матрицы внешних кодов

   FILE* outer_code_inv1 = fopen("./GCCmats/sylvester7sub1.mat", "r");

   for (i = 0; i < 2; i++) {
      for (j = 0; j < 32; j++) {
         fscanf(outer_code_inv1, "%d", &el);
         cdc_mats[0][i * 32 + j] = el;
      }
   }
   fclose(outer_code_inv1);

   FILE* outer_code_inv2 = fopen("./GCCmats/sylvester7sub2.mat", "r");

   for (i = 0; i < 15; i++) {
      for (j = 0; j < 32; j++) {
         fscanf(outer_code_inv2, "%d", &el);
         cdc_mats[1][i * 32 + j] = el;
      }
   }
   fclose(outer_code_inv2);

   FILE* outer_code_inv3 = fopen("./GCCmats/sylvester7sub3.mat", "r");

   for (i = 0; i < 17; i++) {
      for (j = 0; j < 32; j++) {
         fscanf(outer_code_inv3, "%d", &el);
         cdc_mats[2][i * 32 + j] = el;
      }
   }
   fclose(outer_code_inv3);

   FILE* outer_code_inv4 = fopen("./GCCmats/sylvester7sub4.mat", "r");

   for (i = 0; i < 2; i++) {
      for (j = 0; j < 32; j++) {
         fscanf(outer_code_inv4, "%d", &el);
         cdc_mats[3][i * 32 + j] = el;
      }
   }
   fclose(outer_code_inv4);

   grid4[0][0][4] = 1;

   for (j = 0; j < 32; j++) {
      for (i = 0; i < 4; i++) {
         el = (2*cdc_mats[3][2*j] + cdc_mats[3][2*j + 1]) ^ i;
         if (grid4[i][j][4]) {
            grid4[i][j+1][4] = 1;
            grid4[i][j+1][i] = 1;
            grid4[el][j+1][4] = 1;
            grid4[el][j+1][i] = 1;
         }
      }
   }

   // Main simulation loop.
   while (csnrn < sim->snr_num) {

      err_nums1 = 0, err_nums2 = 0, err_nums3 = 0, err_nums4 = 0;

      noise_sg = 1 / sqrt(2 * c_R * db2val(sim->snr_db[csnrn]));

#ifdef DEC_NEEDS_CSNRN
      // Set current SNR # in the codec.
      cdc_set_csnrn(sim->dc_inst, csnrn);
#endif // DEC_NEEDS_CSNRN

#ifdef DEC_NEEDS_SIGMA
      // Set sg in the codec.
      cdc_set_sg(sim->dc_inst, noise_sg);
#endif // DEC_NEEDS_SIGMA

      while (sim->trn[csnrn] < sim->trn_req[csnrn]) {
         
         sim->trn[csnrn]++;

         // Generate random or zero code word.
         // coding
         if (sim->use_rndcw) {
            for (i = 0; i < c_k; i++) x[i] = (RND > 0.5) ? 1 : 0;
            enc_bpsk_gcc(sim->dc_inst, x, c_in, outer_matrix, inner_matrix, out_cd_len, in_cd_len);
         }
         else {
            for (i = 0; i < c_k; i++) x[i] = 0;
            for (i = 0; i < c_n; i++) c_in[i] = -1.0;
         }

         // Generate channel output for the simulated codeword.
         // (Channel simulation.)
         // c_out <-- c_in + noise.
         copy_add_noise(c_in, c_n, noise_sg, c_out);

         // Decode.
         if (dec_bpsk_gcc(sim->dc_inst, c_out, x_dec, inner_matrix_inv, in_cd_len, out_cd_len, cdc_opts,
                          cdc_sizes, cdc_mats, grid4, c_in + 96)) {
            err_msg("sim_run(): Error while decoding.");
            return -1;
         }

         // Count the number of incorrect information bits.
         en = 0;
         for (i = 0; i < c_k-30; i++) {
            if (x_dec[i] != x[i]) en++;
         }
         for (i = 0; i < 32; i++) {
            if (1 - 2 * x_dec[c_k - 30 + i] != c_in[c_n - 32 + i]) en++;
         }

         for (i = 0; i < 2; i++) {
            if (x_dec[i] != x[i]) {
               err_nums1++;
               break;
            }
         }

         for (i = 2; i < 17; i++) {
            if (x_dec[i] != x[i]) {
               err_nums2++;
               break;
            }
         }

         for (i = 17; i < 34; i++) {
            if (x_dec[i] != x[i]) {
               err_nums3++;
               break;
            }
         }


         for (i = 0; i < 32; i++) {
            if (1 - 2 * x_dec[c_k - 30 + i] != c_in[c_n - 32 + i]) {
               err_nums4++;
               break;
            }
         }

         // Adjust error counters.

         if (en) {
            sim->en_bit[csnrn] += en;
            sim->en_bl[csnrn]++;
         }

         // Check elapsed time.
         if (time(NULL) - start_time > sim->ret_int) {
            free(c_out);
            free(c_in);
            free(x_dec);
            free(x);
            free(outer_matrix);
            free(inner_matrix);
            free(inner_matrix_inv);
            free(cdc_opts);
            free(cdc_sizes);
            for (i = 0; i < 4; i++) {
               for (j = 0; j < 33; j++) {
                  free(grid4[i][j]);
               }
               free(grid4[i]);
            }
            free(grid4);
            for (i = 0; i < 4; i++) {
               free(cdc_mats[i]);
            }
            free(cdc_mats);
            return 1;
         }
      }

      printf("%f, %f, %f, %f, %f \n", sim->snr_db[csnrn], (double) err_nums1 / sim->trn_req[csnrn],
                                     (double) err_nums2 / sim->trn_req[csnrn], (double) err_nums3 / sim->trn_req[csnrn],
                                     (double) err_nums4 / sim->trn_req[csnrn]);
      //printf("%d %d %d %d %d\n", err_nums1, err_nums2, err_nums3, err_nums4, sim->en_bl[csnrn]);

      sim->csnrn = ++csnrn;
   }

   free(c_out);
   free(c_in);
   free(x_dec);
   free(x);
   free(outer_matrix);
   free(inner_matrix);
   free(inner_matrix_inv);
   free(cdc_opts);
   free(cdc_sizes);
   for (i = 0; i < 4; i++) {
      for (j = 0; j < 33; j++) {
         free(grid4[i][j]);
      }
      free(grid4[i]);
   }
   free(grid4);
   for (i = 0; i < 4; i++) {
      free(cdc_mats[i]);
   }
   free(cdc_mats);
   return 0;
}
#endif // GCCDEC

// Return short simulation status string.
void
sim_state_str(
   void *inst, // Simulation instance.
   char ds[] // Destination string.
)
{
   sim_bg_inst *sim;
   int n;

   sim = (sim_bg_inst *)inst;
   n = sim->csnrn;
   if (sim->do_ml) {
      sprintf(ds,
         "SNR: %3.2f (%d / %d), trn: %d / %d, en: %d, ML en: %d.",
         sim->snr_db[n], n + 1, sim->snr_num,
         sim->trn[n], sim->trn_req[n],
         sim->en_bl[n], sim->enml_bl[n]
      );
   }
   else {
      sprintf(ds,
         "SNR: %3.2f (%d / %d), trn: %d / %d, en: %d.",
         sim->snr_db[n], n + 1, sim->snr_num,
         sim->trn[n], sim->trn_req[n],
         sim->en_bl[n]
      );
   }
}

// Return current simulation results string.
void
sim_cur_res_str(
   void *inst, // Simulation instance.
   char ds[] // Destination string.
)
{
   sim_bg_inst *sim;
   int n, snrs_to_show;
   char str1[1024];

   sim = (sim_bg_inst *)inst;
   snrs_to_show = (sim->csnrn < sim->snr_num) ? sim->csnrn + 1 : sim->snr_num;
   if (sim->do_ml) {
      sprintf(ds, "SNR\t ep_bit\t\t ep_bl\t\t epml_bl");
      for (n = 0; n < snrs_to_show; n++) {
         sprintf(str1,
            "\n%3.2f\t %.3e\t %.3e\t %.3e",
            sim->snr_db[n],
            (double)(sim->en_bit[n]) / (sim->trn[n] * sim->code_k),
            (double)(sim->en_bl[n]) / sim->trn[n],
            (double)(sim->enml_bl[n]) / sim->trn[n]
         );
         strcat(ds, str1);
      }
   }
   else {
      sprintf(ds, "SNR\t ep_bit\t\t ep_bl");
      for (n = 0; n < snrs_to_show; n++) {
         sprintf(str1,
            "\n%3.2f\t %.3e\t %.3e",
            sim->snr_db[n],
            (double)(sim->en_bit[n]) / (sim->trn[n] * sim->code_k),
            (double)(sim->en_bl[n]) / sim->trn[n]
         );
         strcat(ds, str1);
      }
   }
}

// Save simulation results.
int
sim_save_res(
   void *inst // Simulation instance.
)
{
   sim_bg_inst *sim;
   FILE *fp;
   char bsyfn[FN_LEN_MAX]; // Busy flag file name.
   char *str1, *token; // For old res file parsing.
   int snr_num = 0;
   double snr_db[SNR_NUM_MAX];
   int trn[SNR_NUM_MAX];
   int ml_trn[SNR_NUM_MAX];
   int en_bit[SNR_NUM_MAX];
   int en_bl[SNR_NUM_MAX];
   int enml_bl[SNR_NUM_MAX];
   int trn_n = 0;
   int ml_trn_n = 0;
   int en_bit_n = 0;
   int en_bl_n = 0;
   int enml_bl_n = 0;
   int i, j, n, snrs_to_save;

   sim = (sim_bg_inst *)inst;

   // Check if the file is busy.
   strcpy(bsyfn, sim->rf_name);
   strcat(bsyfn, ".bsy");
   n = 0;
   while ((fp = fopen(bsyfn, "r")) != NULL) {
      fclose(fp);
      if (n++ > 30) {
         err_msg("sim_save_res(): Output file is busy too long.");
         return 1;
      }
      sleep(1000);
   }

   // Create the busy flag.
   if ((fp = fopen(bsyfn, "w")) == NULL) {
      err_msg("sim_save_res(): Cannot create busy flag file.");
      return 1;
   }
   fprintf(fp, "Busy");
   fclose(fp);

   // Try to read and preparse old res file.
   if (spf_tryread_preparse(sim->rf_name, &str1) == 0) {
      
      // Get old data.
      token = strtok(str1, tk_seps_prepared);
      while (token != NULL) {
         // TRYGET_INT_TOKEN(token, code_n_token, code_n);
         // TRYGET_INT_TOKEN(token, code_k_token, code_k);
         TRYGET_GRDOUBLE_TOKEN(
            token, SNR_token, snr_db, snr_num, SNR_NUM_MAX,
            "Error in old res file.");
         TRYGET_GRINT_TOKEN(
            token, tr_num_token, trn, trn_n, SNR_NUM_MAX,
            "Error in old res file.");
         TRYGET_GRINT_TOKEN(
            token, en_bit_token, en_bit, en_bit_n, SNR_NUM_MAX,
            "Error in old res file.");
         TRYGET_GRINT_TOKEN(
            token, en_bl_token, en_bl, en_bl_n, SNR_NUM_MAX,
            "Error in old res file.");
         TRYGET_GRINT_TOKEN(
            token, ml_tr_num_token, ml_trn, ml_trn_n, SNR_NUM_MAX,
            "Error in old res file.");
         TRYGET_GRINT_TOKEN(
            token, enml_bl_token, enml_bl, enml_bl_n, SNR_NUM_MAX,
            "Error in old res file.");
         SPF_SKIP_UNKNOWN_PARAMETER(token);
      }
      free(str1);

      // Check validity of the data.
      if ((snr_num != trn_n) || (en_bit_n != trn_n) || (en_bl_n != trn_n)
         || (enml_bl_n != trn_n) || (ml_trn_n != trn_n)) {
         err_msg("sim_save_res(): Invalid old res file.");
         return 1;
      }

      // Update the data.
      for (i = 0; (i <= sim->csnrn) && (i < sim->snr_num); i++) {
         if (sim->trn[i] == sim->trn_saved[i]) continue;
         n = 0;
         while ((snr_db[n] < sim->snr_db[i]) && (n < snr_num)) n++;
         if ((snr_db[n] > sim->snr_db[i]) || (n >= snr_num)) {
            // Make room for the new at n if required.
            for (j = snr_num; j > n; j--) {
               snr_db[j] = snr_db[j - 1];
               trn[j] = trn[j - 1];
               en_bit[j] = en_bit[j - 1];
               en_bl[j] = en_bl[j - 1];
               ml_trn[j] = ml_trn[j - 1];
               enml_bl[j] = enml_bl[j - 1];
            }
            snr_num++;
            // Add new at n.
            snr_db[n] = sim->snr_db[i];
            trn[n] = 0;
            en_bit[n] = 0;
            en_bl[n] = 0;
            ml_trn[n] = 0;
            enml_bl[n] = 0;
         }
         // Update at n.
         trn[n] += sim->trn[i] - sim->trn_saved[i];
         en_bit[n] += sim->en_bit[i] - sim->en_bit_saved[i];
         en_bl[n] += sim->en_bl[i] - sim->en_bl_saved[i];
         if (sim->do_ml) {
            ml_trn[n] += sim->trn[i] - sim->trn_saved[i];
            enml_bl[n] += sim->enml_bl[i] - sim->enml_bl_saved[i];
         }
      }
      snrs_to_save = snr_num;
   }
   else { 
      
      // Assume there was no old res file.
      snrs_to_save = (sim->csnrn < sim->snr_num) ? sim->csnrn + 1 : sim->snr_num;
      for (n = 0; n < snrs_to_save; n++) {
         snr_db[n] = sim->snr_db[n];
         trn[n] = sim->trn[n] - sim->trn_saved[n];
         sim->trn_saved[n] = sim->trn[n];
         en_bit[n] = sim->en_bit[n] - sim->en_bit_saved[n];
         sim->en_bit_saved[n] = sim->en_bit[n];
         en_bl[n] = sim->en_bl[n] - sim->en_bl_saved[n];
         sim->en_bl_saved[n] = sim->en_bl[n];
         if (sim->do_ml) {
            ml_trn[n] = trn[n];
            enml_bl[n] = sim->enml_bl[n] - sim->enml_bl_saved[n];
            sim->enml_bl_saved[n] = sim->enml_bl[n];
         }
         else {
            ml_trn[n] = 0;
            enml_bl[n] = 0;
         }
      }
   }
   for (i = 0; (i <= sim->csnrn) && (i < sim->snr_num); i++) {
      sim->trn_saved[i] = sim->trn[i];
      sim->en_bit_saved[i] = sim->en_bit[i];
      sim->en_bl_saved[i] = sim->en_bl[i];
      sim->enml_bl_saved[i] = sim->enml_bl[i];
   }

   // Save updated data.

   fp = fopen(sim->rf_name, "wt");
   if (fp == NULL) {
      err_msg("sim_save_res(): Cannot open file to save results.");
      return 1;
   }

   fprintf(fp, "code_n %d\ncode_k %d\n\n", sim->code_n, sim->code_k);

   fprintf(fp, "SNR { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%g ", snr_db[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "tr_num { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", trn[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "en_bit { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", en_bit[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "en_bl { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", en_bl[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "ml_tr_num { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", ml_trn[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "enml_bl { ");
   for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", enml_bl[n]);
   fprintf(fp, "}\n");

   fprintf(fp, "\n");
   fprintf(fp, "%% SNR   ep_bit      ep_bl       epml_bl");
   for (n = 0; n < snrs_to_save; n++) {
      fprintf(fp,
         "\n%% %3.2f  %.3e  %.3e  %.3e",
         snr_db[n],
         (double)(en_bit[n]) / ((double)trn[n] * sim->code_k),
         (double)(en_bl[n]) / trn[n],
         (ml_trn[n]) ? (double)(enml_bl[n]) / ml_trn[n] : 0.0
      );
   }

   // Delete busy flag file.
   remove(bsyfn);

   fclose(fp);
   return 0;
}

// Control simulation parameters "on the fly".
void sim_control(
   void *inst, // Simulation instance.
   int ctrl_code
)
{
   sim_bg_inst *sim;
   int csnrn; // Current SNR value number.

   sim = (sim_bg_inst *)inst;
   csnrn = sim->csnrn;

   switch (ctrl_code) {
   case SIM_CTRL_CUR_MORE :
      sim->trn_req[csnrn] += sim->trn_req[csnrn] / 2;
      break;
   case SIM_CTRL_CUR_LESS :
      sim->trn_req[csnrn] -= sim->trn_req[csnrn] / 3;
      if (sim->trn_req[csnrn] < sim->trn[csnrn]) sim->trn_req[csnrn] = sim->trn[csnrn];
      break;
   case SIM_CTRL_CUR_NEXT : sim->trn_req[csnrn] = sim->trn[csnrn];
   }
}