# -*- coding: utf-8 -*-
"""
Wrapper around C++ implementation of the generative models.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import subprocess
import numpy as np
import time
from os import path


def _unpack(model, model_params):
    if model == "gn":
        return [model_params["gamma"], 1]
    if model == "generalized_gn":  # we use the full generalized GN with a constant birth probability.
        return [0, 0, model_params["b"],
                1, 1, abs(model_params["gamma"]), 0, int((np.sign(-model_params["gamma"]) + 1) / 2)]


def run(model, model_params, T, seed=None, verbose=False,
        path_to_generator_exec="../generators/growth",
        simplified_interface=True):
    """Generate a random instance of a generative model.

    Parameters
    ----------
    model : str
      Name of the generative model.
      Implemented models are: gn, generalized_gn.
    model_params : dict
      Parameters of the generative model.
    T : int
      Number of generative step.
    seed : int (optional)
      Specific seed of the generative model.
    verbose : bool
      Output logs to stdout.
    path_to_generator_exec : str / path
      Relative path to compiled executable of the generative models.
    simplified_interface : bool
      Assume that the generator is compiled with a simplified interface (i.e., not Boost).

    Returns
    -------
    X : list of tuples
      History vector. Tuples represent the edges of a multigraph (with a tag).
    """
    # Construct c++ call
    absolute_path = path.dirname(path.abspath(__file__))

    if not simplified_interface:
        call = [path.join(absolute_path, path_to_generator_exec)]
        call += ["-m", model,
                 "-p", *[str(p) for p in _unpack(model, model_params)],
                 "-t", str(T)]
        if seed is not None:
            call += ['-d', str(seed)]
    else:
        if seed is None:
            raise Exception("Random seed must be specified when the " +
                            "generator uses the simplified interface.")
        call = [path.join(absolute_path, path_to_generator_exec)]
        call += [model, str(T), str(seed), *[str(p) for p in _unpack(model, model_params)]]
    if verbose:
        print(" ".join(call))
        proc = subprocess.run(call, stdout=subprocess.PIPE)
    else:
        proc = subprocess.run(call,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
    # Construct and return captured output
    X = []
    edge_tags = dict()
    for edge in proc.stdout.decode().strip().split("\n"):
        e = edge.strip().split()
        e = (int(e[0]), int(e[1]))
        if edge_tags.get(e) is None:
            edge_tags[e] = 0
        else:
            edge_tags[e] += 1
        X.append((e[0], e[1], edge_tags[e]))
    return X
