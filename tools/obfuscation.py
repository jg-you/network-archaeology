# -*- coding: utf-8 -*-
"""
Obfuscation tools for avoiding accidental biases in inference methods.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import random


def obfuscate_history(X, max_id=None, seed=None):
    """Take an history vector and return an obfuscated version.

    The obfuscated version has the same structure but nodes are randomly
    mapped to new nodes, as to avoid any accidental bias. Edges are also
    randomly re-ordered and re-tagged.

    Parameters
    ----------
    X : list of tuples
      History vector. Tuples represent the edges of a multigraph (with a tag).
    max_id: int
      If set, skip the calculation of the set of nodes and
      assumes that nodes are 0-indexed and contiguous up to max_id.

    Return
    ------
    Y : list of tuples
      Obfuscated history vector.
    encoding: dict
      Information required to reverse the node obfuscation.
      Maps nodes identifiers to their original.
    tag_encoding: dict of dict
      Information required to reverse the edge tag obfuscation.
      Maps nodes edge tags to their original.
      The outer dict key are un-tagged edge.
      The inner dict key are tags.
    """
    # Construct array of nodes, and create a random ordering.
    if max_id is not None:
        nodes = list(range(max_id + 1))
    if max_id is None:
        nodes = list({v for edge in X for v in edge})  # lazy, could optimize
    if seed is not None:
        random.seed(seed)
    random.shuffle(nodes)
    # Construct encoding dictionary
    encoding = {nodes[idx]: idx for idx in range(len(nodes))}
    # Construct obfuscated history
    partial = [0 for _ in range(len(X))]
    for t, e in enumerate(X):
        e0 = nodes[e[0]]
        e1 = nodes[e[1]]
        partial[t] = (e0, e1, e[2]) if e0 < e1 else (e1, e0, e[2])
    random.shuffle(partial)
    # Obfuscate tags
    Y = [0 for _ in range(len(X))]
    tag_encoding = dict()
    for t, e in enumerate(partial):
        if tag_encoding.get((e[0], e[1])) is None:
            tag_encoding[(e[0], e[1])] = {e[2]: 0}
        else:
            tag_encoding[(e[0], e[1])][e[2]] = len(tag_encoding.get((e[0], e[1])))
        Y[t] = (e[0], e[1], tag_encoding[(e[0], e[1])][e[2]])
    return Y, encoding, tag_encoding


def deobfuscate_history(Y, encoding, tag_encoding):
    """De-obfuscate an history, given the encoding.

    Parameters
    ----------
    Y : list of tuples
      Obfuscated history vector.
    encoding: dict
      Information required to reverse the node obfuscation.
      Maps nodes identifiers to their original.
    tag_encoding: dict of dict
      Information required to reverse the edge tag obfuscation.
      Maps nodes edge tags to their original.
      The outer dict key are un-tagged edge.
      The inner dict key are tags.

    Return
    ------
    X : list of tuples
      De-obfuscated history vector.
    """
    X = [0 for _ in range(len(Y))]
    for t, e in enumerate(Y):
        e0 = encoding[e[0]]
        e1 = encoding[e[1]]
        e2 = tag_encoding[(e[0], e[1])][e[2]]
        X[t] = (e1, e0, e2) if e0 > e1 else (e0, e1, e2)
    return X
