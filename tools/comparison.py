# -*- coding: utf-8 -*-
"""
Comparison tools of history vectors.

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import scipy as sp


def corr(X, Y):
    """Compare two histories event by event and give a similarity score.

    Warning
    -------
    Note the asymmetry of X and Y; the latter is inferred and can therefore
    contain ties. We add an additional variable to denote the time
    of birth of an edge.

    Parameters
    ----------
    X : list of tuples
      Reference history vector; tuples represent edges.
      Position corresponds to time.
    Y : list of pairs
      Inferred history vector with ranking information.
      The first entry of the pair contains an edge (pair)
      The second entry contains the rank of the edge (float).

    Return
    ------
    score : float
      Correlation of of generated and infered history.
    """
    # Augment reference history with arrival times
    X = [(_, t) for t, _ in enumerate(X)]
    # Sort based on edges
    X = sorted(X, key=lambda x: x[0])
    Y = sorted(Y, key=lambda x: x[0])
    corr = sp.corrcoef([x[1] for x in X], [y[1] for y in Y])[0, 1]
    return corr
