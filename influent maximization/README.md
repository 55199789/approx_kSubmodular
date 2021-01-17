Compiled with g++ 7.3.0 under Ubuntu Linux release 16.04.10

# To compile:
g++ -std=c++11 -O3 approx_subm.cpp jlog.cpp -o approx_subm

# Examples:
## ADR 
./approx_subm -algo=lazy -n=55 -k=3 -tick=65535 -budget=13 -0=./dat/temps.tsv -1=./dat/hums.tsv -2=./dat/lights.tsv -repeats=1 -epsilon=1 -seed=532973 
## AS
./approx_subm -algo=lazy_H -n=55 -k=3 -tick=65535 -budget=13 -0=./dat/temps.tsv -1=./dat/hums.tsv -2=./dat/lights.tsv -repeats=1 -epsilon=1 -seed=532973 

# Parameters
+ algo is the algorithm: 
  - ADR: run greedy algorithm on ADR functions
  - AS: run greedy algorithm on AS functions
+ n: the number of locations (i.e., size of set |V|): fixed to 55 in our experiments
+ budget: B_1=B_2=B_3
+ k: the number of dimensions (1, 2, 3)
+ tick: the number of observations in data
+ 0 1 2: the datasets path (sensor measurements for sensors of types 0, 1, 2)
+ repats: repeats is a number, how many times we will execute lazy_greedy. You can always give 1. If you give a number A, it runs it A times and computes the average H.
+ epsilon: is the parameter epsilon
+ seed: random seed




