clear;
addpath('Decoding_Index/')
addpath('GA/')
addpath('GetCriticalSet/')
n = 6;
N = 2^n;
K = 38;
max_iter = 50;
max_err = 20;
max_runs = 1e7;
resolution = 1e5;
ebno_vec = 0:0.2:3;
[bler, ber] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K);


