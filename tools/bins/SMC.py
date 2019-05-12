#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SMC sampling for generalized PA with densification.

Author: Jean-Gabriel Young <jgyou@umich.edu>
"""
import networkx as nx
import numpy as np
import copy
from SamplableSet import EdgeSamplableSet

# Global status dict
[UNEXPLORED, BOUNDARY, EXPLORED] = range(3)
[TENDRIL, CLOSURE] = range(2)


def order(edge):
    """Canonical ordering of an edge's representation."""
    if edge[0] > edge[1]:
        return (edge[1], edge[0], edge[2])
    else:
        return edge


class History(object):
    """Edge history of a MultiGraphs."""

    def __init__(self, g, gamma, b, seed, use_snowball=True):
        """History constructor.

        Parameters
        ----------
        g : nx.MultiGraph
            Graph whose history we will sample from.
        gamma : float
            Exponent of the attachment Kernel.
        b : float in [0, 1]
            Closure probability.
        seed : int
            RNG seed
        use_snowball : bool
            Use snowball proposal distribution.
        """
        self.t = 0
        self.use_snowball = use_snowball
        self.gamma = gamma
        self.b = b
        self.X = [(None, None, None)] * g.number_of_edges()
        self.deg = [0] * g.number_of_nodes()
        self.Z = 0
        self.k_max = max(list(nx.degree(g)), key=lambda x: x[1])[1]
        if use_snowball:
            self.step = self._step_snowball
            self.closures = EdgeSamplableSet(1, 1, seed)
            self.tendrils = EdgeSamplableSet(1, 1, seed)
        else:
            self.step = self._step_weighted
            if gamma >= 0:
                self.closures = EdgeSamplableSet(1, self.k_max ** (2 * gamma), seed)  # noqa
                self.tendrils = EdgeSamplableSet(1, self.k_max ** gamma, seed)
            else:
                self.closures = EdgeSamplableSet(self.k_max ** (2 * gamma), 1, seed)  # noqa
                self.tendrils = EdgeSamplableSet(self.k_max ** gamma, 1, seed)
        self.status = {order(e): UNEXPLORED for e in g.edges}

    # Private members
    def _sample_boundary_unif(self):
        """Sample an edge in the boundary."""
        w_tendrils = self.tendrils.size()
        w_closures = self.closures.size()
        if np.random.rand() < w_tendrils / (w_tendrils + w_closures):
            # sample from tendrils
            (e, _) = self.tendrils.sample()
            return (e,
                    w_tendrils + w_closures,
                    TENDRIL)
        else:
            # sample from closures
            (e, _) = self.closures.sample()
            return (e,
                    w_tendrils + w_closures,
                    CLOSURE)

    def _sample_boundary_weighted(self):
        """Sample an edge in the boundary."""
        w_tendrils = self.b * self.tendrils.total_weight()
        w_closures = (1 - self.b) * self.closures.total_weight() / self.Z
        if np.random.rand() < w_tendrils / (w_tendrils + w_closures):
            # sample from tendrils
            (e, _) = self.tendrils.sample()
            return (e,
                    (w_tendrils + w_closures) / self.Z,
                    TENDRIL)
        else:
            # sample from closures
            (e, _) = self.closures.sample()
            return (e,
                    (w_tendrils + w_closures) / self.Z,
                    CLOSURE)

    def _update_degrees(self, edge):
        """Update degree vector and norm upon adding `edge`."""
        v0 = edge[0]
        v1 = edge[1]
        if self.deg[v0] == 0 and self.deg[v1] == 0:
            pass
        elif self.deg[v0] == 0:    # new node is v0
            self.Z -= self.deg[v1] ** self.gamma
            self.Z += (self.deg[v1] + 1) ** self.gamma + 1
        elif self.deg[v1] == 0:  # new node is v1
            self.Z -= self.deg[v0] ** self.gamma
            self.Z += (self.deg[v0] + 1) ** self.gamma + 1
        else:
            if v1 == v0:  # count self-loops
                self.Z -= self.deg[v0] ** self.gamma
                self.Z += (self.deg[v0] + 2) ** self.gamma
            else:
                self.Z -= self.deg[v0] ** self.gamma +\
                    self.deg[v1] ** self.gamma
                self.Z += (self.deg[v0] + 1) ** self.gamma +\
                    (self.deg[v1] + 1) ** self.gamma
        self.deg[v0] += 1
        self.deg[v1] += 1

    def _get_neighboring_edges(self, g, edge):
        """Get edges neighboring `edge`, excluding itself."""
        x = [(edge[0], v, _) for v in nx.neighbors(g, edge[0]) for _
             in range(g.number_of_edges(edge[0], v)) if v != edge[1]] +\
            [(edge[0], edge[1], tag) for tag in
             range(g.number_of_edges(edge[0], edge[1])) if tag != edge[2]]
        if edge[0] != edge[1]:  # not a self-loop
            x += [(edge[1], v, _) for v in nx.neighbors(g, edge[1]) for _
                  in range(g.number_of_edges(edge[1], v)) if v != edge[0]]
        return set([order(e) for e in x])

    def _is_closure(self, e):
        return self.deg[e[0]] > 0 and self.deg[e[1]] > 0

    def _is_tendril(self, e):
        return self.deg[e[0]] == 0 or self.deg[e[1]] == 0

    # Public methods
    def init(self, init_edge):
        """Initialize history at init_edge."""
        self.t = 0
        self.tendrils.insert(init_edge, 1)
        self.status[init_edge] = BOUNDARY
        self.Z = 2

    def _step_snowball(self, g):
        """Step the sampler forward (snowball proposal)."""
        # Sample boundary
        (edge, step_weight, e_type) = self._sample_boundary_unif()
        if e_type == TENDRIL:
            w = 1
            if self.deg[edge[0]] > 0:
                w = self.deg[edge[0]] ** self.gamma
            elif self.deg[edge[1]] > 0:
                w = self.deg[edge[1]] ** self.gamma
            step_weight *= w / self.Z * self.b
        else:
            w = self.deg[edge[0]] ** self.gamma * self.deg[edge[1]] ** self.gamma  # noqa
            step_weight *= w / (self.Z ** 2) * (1 - self.b)

        # Update history's internal state
        self.X[self.t] = edge
        self.t += 1
        self._update_degrees(edge)

        # Update the status of `edge`
        self.status[edge] = EXPLORED
        if e_type == TENDRIL:
            self.tendrils.erase(edge)
        else:
            self.closures.erase(edge)

        # Update affected edges' status and weights
        for e in self._get_neighboring_edges(g, edge):
            status = self.status[e]
            if status == EXPLORED:
                continue
            # The edge is either UNEXPLORED or in the BOUNDARY
            if self._is_tendril(e):
                if status == UNEXPLORED:
                    self.status[e] = BOUNDARY
                    self.tendrils.insert(e, 1)
                else:
                    self.tendrils.set_weight(e, 1)
            else:
                if status == UNEXPLORED:
                    self.status[e] = BOUNDARY
                    self.closures.insert(e, 1)
                else:
                    # Was it in the correct boundary?
                    if self.tendrils.count(e) > 0:
                        self.tendrils.erase(e)
                        self.closures.insert(e, 1)
                    else:
                        self.closures.set_weight(e, 1)
        return step_weight

    def _step_weighted(self, g):
        """Step the sampler forward (truncated posterior proposal)."""
        # Sample boundary
        (edge, step_weight, e_type) = self._sample_boundary_weighted()

        # Update history's internal state
        self.X[self.t] = edge
        self.t += 1
        self._update_degrees(edge)

        # Update the status of `edge`
        self.status[edge] = EXPLORED
        if e_type == TENDRIL:
            self.tendrils.erase(edge)
        else:
            self.closures.erase(edge)

        # Update affected edges' status and weights
        for e in self._get_neighboring_edges(g, edge):
            status = self.status[e]
            if status == EXPLORED:
                continue
            # The edge is either UNEXPLORED or in the BOUNDARY
            if self._is_tendril(e):
                if self.deg[e[0]] > 0:
                    w = self.deg[e[0]] ** self.gamma
                else:
                    w = self.deg[e[1]] ** self.gamma
                if status == UNEXPLORED:
                    self.status[e] = BOUNDARY
                    self.tendrils.insert(e, w)
                else:
                    self.tendrils.set_weight(e, w)
            else:
                w = self.deg[e[0]] ** self.gamma *\
                    self.deg[e[1]] ** self.gamma
                if status == UNEXPLORED:
                    self.status[e] = BOUNDARY
                    self.closures.insert(e, w)
                else:
                    # Was it in the correct boundary?
                    if self.tendrils.count(e) > 0:
                        self.tendrils.erase(e)
                        self.closures.insert(e, w)
                    else:
                        self.closures.set_weight(e, w)
        return step_weight

    def duplicate(self, seed):
        """Return a copy of the History with new RNG seed."""
        cls = self.__class__
        new_history = cls.__new__(cls)
        new_history.use_snowball = self.use_snowball
        if self.use_snowball:
            new_history.step = new_history._step_snowball
        else:
            new_history.step = new_history._step_weighted
        new_history.t = self.t
        new_history.gamma = self.gamma
        new_history.b = self.b
        new_history.X = copy.copy(self.X)
        new_history.deg = copy.copy(self.deg)
        new_history.Z = self.Z
        new_history.k_max = self.k_max
        new_history.closures = EdgeSamplableSet(self.closures, seed)
        new_history.tendrils = EdgeSamplableSet(self.tendrils, seed)
        new_history.status = copy.copy(self.status)
        return new_history

    def tau(self, e):
        """Postition of edge 'e' in the history."""
        if len(e) == 2:
            e = (e[0], e[1], 0)  # allow for incorrect formatting
        return self.X.index(order(e))


class BridgeSampler(object):
    """Implements Bloem-Reddy & Orbanz Bridge sampling algorithm."""

    def __init__(self):
        """Empty constructor."""
        pass

    def draw(self, g, gamma, b, n, min_ESS, use_snowball=True):
        """Generate a batch of samples with adaptive SMC.

        Parameters
        ----------
        g : nx.MultiGraph
            Graph whose history we will sample from.
        gamma : float
            Exponent of the attachment Kernel.
        b : float in [0, 1]
            Closure probability.
        n : int
            Number of samples
        min_ESS : float
            Lower bound on the ESS that will not trigger resampling.
        use_snowball : bool
            Use snowball proposal distribution.
        """
        # First step
        T = g.number_of_edges()
        weight_list = np.ones(n)
        history_list = [History(g, gamma, b, i, use_snowball) for i in range(n)]  # noqa
        for i in range(n):
            # selection somewhat complicated by the presence of self-loops
            id_init = np.random.randint(0, T - g.number_of_selfloops())
            idx = -1
            for e in g.edges:
                if e[0] != e[1]:
                    idx += 1
                if idx == id_init:
                    init_edge = order(e)
                    break
            history_list[i].init(order(init_edge))
            history_list[i].step(g)  # mark first edge as explored

        # Bulk of the steps
        for t in range(1, T):
            for i in range(n):
                weight_list[i] *= history_list[i].step(g)
            # check for ESS
            ESS = np.sum(weight_list) ** 2 / np.sum(weight_list ** 2)

            if ESS < min_ESS:
                # Re-sample indices
                weight_list /= np.sum(weight_list)
                indices = np.random.choice(range(n),
                                           size=n,
                                           replace=True,
                                           p=weight_list)
                tmp_list = [None for i in range(n)]
                for i in range(n):
                    new_seed = np.random.randint(2147483647)   # max int32
                    tmp_list[i] = history_list[indices[i]].duplicate(new_seed)
                    weight_list[i] = 1
                history_list = tmp_list
                del tmp_list
        return history_list, weight_list


if __name__ == '__main__':
    # Read graph
    import sys
    g = nx.read_edgelist(sys.argv[1],
                         nodetype=int,
                         create_using=nx.MultiGraph(),
                         data=[("key", int)])
    gamma = float(sys.argv[2])
    b = float(sys.argv[3])
    n = int(sys.argv[4])
    min_ESS = float(sys.argv[5])
    if int(sys.argv[6]) == 1:
        use_snowball = False
    else:
        use_snowball = True

    sampler = BridgeSampler()
    hl, wl = sampler.draw(g, gamma, b, n, min_ESS, use_snowball)

    # Construct raw arrival estimators
    edge_list = g.edges(keys=True)
    tau = [np.mean([h.tau(order(e)) * w
           for h, w in zip(hl, wl)]) for e in edge_list]
    T = g.number_of_edges()
    norm = (2 * np.sum(tau)) / (T * (T - 1))

    # Pool multi edges
    unique_edges = {(order(e)[0], order(e)[1]): [] for e in edge_list}
    for e, t in zip(edge_list, tau):
        ep = order(e)
        unique_edges[(ep[0], ep[1])].append((e[2], t / norm))

    # Out
    for base_edge in unique_edges:
        avg = np.mean([t for (v, t) in unique_edges[base_edge]])
        for (v, _) in unique_edges[base_edge]:
            print(base_edge[0], base_edge[1], v, avg)
