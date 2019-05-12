# -*- coding: utf-8 -*-
"""
Example of script.

This script sweeps values of gamma at a fixed value of b.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
if __name__ == '__main__':
    # Relative paths
    import sys
    from os import path, remove
    sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)),
                    "tools/"))
    # Import
    import consistency_check as cc
    import sys
    import time

    # Global settings
    method_params = dict()
    model_params = dict()

    outfile = "test.txt"
    method = "SMC"  # options : OD, degree, SMC
    b = 1  # any float in [0, 1]
    T = 50  # 50 edges
    if method == "SMC":
        method_params["num_samples"] = 1000
        method_params["min_ESS"] = 500
        method_params["use_truncated"] = True
    seed = 42

    # Parameter grid
    method_params["b"] = b
    model_params["b"] = b
    gamma_range = [-10, -7.5, -5, -3,
                   -2, -1, -0.5, -0.25, -0.05,
                   0, 0.05, 0.25, 0.50, 0.75, 0.95,
                   1, 1.05, 1.25, 1.50, 1.75, 1.95,
                   2]

    # Output header
    with open(outfile, 'a') as f:
        print("#model_params\tT\tmethod\tmethod_params\tscore", file=f)
    
    # Sweep gammas.
    for gamma in gamma_range:
        model_params["gamma"] = gamma
        method_params["gamma"] = gamma
        tmp_path = "/tmp/varying_kernell_" +\
                   str(int(time.time())) + ".txt"
        try:
            scores = cc.run("generalized_gn",
                            model_params,
                            T,
                            method,
                            method_params,
                            1,  # num_iter,
                            tmp_path=tmp_path,
                            seed=seed,
                            verbose=False,
                            simplified_interface=True)
            for s in scores:
                with open(outfile, 'a') as f:
                    print(model_params,
                          T,
                          method,
                          method_params,
                          s,
                          sep='\t',
                          file=f)
            print("gamma =", gamma, "|", "correlation:", scores[0])
        except Exception as e:
            print("#", str(e))  # log exceptions
            remove(tmp_path)
            pass
