# -*- coding: utf-8 -*-
"""
Complete consitency-check pipeline.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import generative_models as gn
import inference_methods as im
import comparison as cp
from obfuscation import obfuscate_history, deobfuscate_history
from os import remove


available_models = {"gn", "generalized_gn"}
available_methods = {"degree", "OD",
                     "snowball_sampling"}


def _write_history(X, path):
    with open(path, 'w') as f:
        for edge in X:
            print(str(edge[0]), str(edge[1]), str(edge[2]), file=f)


def run(model, model_params, T,
        method, method_params, num_iter,
        tmp_path="/tmp/consistency_check.txt",
        seed=None,
        verbose=False,
        simplified_interface=True):
    """
    Wrapper around the full consistency check pipeline.

    Parameters
    ----------
    model : str
      Name of the generative model.
      Implemented models are: gn, generalized_gn.
    model_params : dict
      Parameters of the generative model.
    T : int
      Number of generative step.
    method : str
      Name of the inference method.
      Implemented methods are: degree, OD, random_expand,
                               snowball_sampling, biased_snowball_sampling.
    method_params : dict
      Parameters of the inference method.
    num_iter : int
      Number of repetition of the inference method.
      Note that all repetitions run on the same model instance.
    tmp_path : str
      Location where temporary files will be written.
    verbose : bool
      Output logs to stdout.
    simplified_interface : bool
      Assume that the generator is compiled with a simplified interface (i.e., not Boost).

    Returns
    -------
    scores : list of dict
        A list of scores (one per repetition).
        Each entry of the list corresponds to a repetition of the method.
        An entry in the list is a dictionary, whose key is the name of
        the comparison measure.

    Warning
    -------
    This function has side-effects. It writes and read from a temporary
    location (defaulted to /tmp/) to communicate with pre-compiled modules.
    If multiple instances run at the same time, make sure to pass different
    temporary paths to each instances.
    """
    # Tests
    if {model} & available_models == set():
        raise NotImplementedError("Model '" + str(model) +
                                  "' not implemented.")
    if {method} & available_methods == set():
        raise NotImplementedError("Method '" + str(method) +
                                  "' not implemented.")
    # Generate history
    generated_history = gn.run(model, model_params, T,
                               verbose=verbose, seed=seed,
                               simplified_interface=simplified_interface)
    encoded_history, encoding, tag_encoding = obfuscate_history(generated_history, seed=seed)
    _write_history(encoded_history, tmp_path)
    # Infer and compute similarity
    scores = []
    for i in range(num_iter):
        output = im.run(tmp_path, method, method_params, verbose=verbose)
        if len(generated_history) != len(output):
            RuntimeError("Length of generated and inferred data don't match.")
        inferred = deobfuscate_history([x[0] for x in output], encoding, tag_encoding)
        res = cp.corr(generated_history,
                      [(e, _[1]) for e, _ in zip(inferred, output)])
        scores.append(res)
    # Garbage collection
    remove(tmp_path)
    return scores
