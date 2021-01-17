Compiled with g++ 7.3.0 under Ubuntu Linux release 16.04.10

# To compile
g++ -std=c++11 -O3 greedy.cpp jlog.cpp -o greedy

# Examples:
## ADR
./greedy -algo=ADR -edge=./dat/digg_graph_100.tsv -prob=./dat/digg_prob_100_1000_pow10.tsv -k=5 -budget=100 -epsilon=1 -beta=100 -seed=1435837 -method=max
## AS 
./greedy -algo=AS -edge=./dat/digg_graph_100.tsv -prob=./dat/digg_prob_100_1000_pow10.tsv -k=5 -budget=100 -epsilon=1 -beta=100 -seed=1435837 -method=sum

# Parameters
+ algo is the algorithm:
  - ADR: run greedy algorithm on ADR functions
  - AS: run greedy algorithm on AS functions
+ edge: edge file path
+ prob: prob file path
+ k: 1, 2, 3, ..., or 10, the number of topics to run
+ budget: total size constraint $B$
+ epsilon: approximation ratio
+ beta: number of simulations
+ seed: random seed
+ method:
    - sum: $\frac{\sum_{x\in \pmb{x}}\xi(x)}{|supp(\pmb{x})|}$
    - max: $\max_{x\in\pmb{x}}\xi(x)$
    - new: randomly generated $\xi(\pmb{x})$
