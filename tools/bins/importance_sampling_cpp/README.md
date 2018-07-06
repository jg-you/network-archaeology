Compilation
-----------
In the random_sampling_cpp directory:

    cmake .;
    make

Execution
---------
note that the executable are in the parent dir (execs). General parameters are :
1. edge list file path (string).
2. model parameter gamma (double).
3. model parameter b (double).
4. Total size (int) of the desired sample.
5. Bias exponent (double). Let us denote alpha for the exponent and k1,k2 the degree of the endpoints of an edge. This edge is chosen proportionnally to (k1k2)^alpha as the seed. alpha = 0 corresponds to choosing an edge randomly among the edge list.
6. (optional) seed (int) for the random number generator.
