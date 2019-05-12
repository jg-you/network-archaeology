# -*- coding: utf-8 -*-
"""
Wrapper around C++/external implementation of the inference methods.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import subprocess
from os import path


def _build_call(history_file_path, method, method_params):
    absolute_path = path.dirname(path.abspath(__file__))
    if method == "OD":
        call = [path.join(absolute_path, "bins", "OD")]
        call += [history_file_path]
    if method == "degree":
        call = [path.join(absolute_path, "bins", "degree.py")]
        call += [history_file_path]
    if method == "SMC":
        call = [path.join(absolute_path, "bins", "SMC.py")]
        call += [history_file_path]
        call += [str(method_params["gamma"])]
        call += [str(method_params["b"])]
        call += [str(method_params["num_samples"])]
        call += [str(method_params["min_ESS"])]
        if method_params.get('use_truncated') is not None:
            if method_params.get('use_truncated') is True:
                call += ["1"]
            else:
                call += ["0"]
        else:
            call += ["0"]
    return call


def run(history_file_path, method, method_params=None, verbose=False):
    """Infer history from an instance of a generative model.

    Parameters
    ----------
    history_file_path : str / path
      Path of the original history file.
    method : str
      Name of the inference method.
      Implemented methods are: degree, OD, random_expand,
                               snowball_sampling, biased_snowball_sampling.
    method_params : dict (optional)
      Parameters of the inference method, if any.
    verbose : bool
      Output logs to stdout.

    Returns
    -------
    X : list of tuples
      Inferred history vector.
      The first entry of the pair is an edge (with tags).
      The second entry of the pair is the estimated arrival time of the edge.

    Notes
    -----
    The history is not passed to the inference method as an object, but instead
    through a file (previously written to disk). The rationale is that this
    simplifies communication among processes and avoid passing the same
    information multiple times.
    """
    # Construct call
    call = _build_call(history_file_path, method, method_params)
    if verbose:
        print(" ".join(call))
        proc = subprocess.run(call, stdout=subprocess.PIPE)
    else:
        proc = subprocess.run(call,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
    # Construct and return captured output
    X = []
    for line in proc.stdout.decode().strip().split("\n"):
        data = line.strip().split()
        X.append(((int(data[0]), int(data[1]), int(data[2])), float(data[3])))
    return X
