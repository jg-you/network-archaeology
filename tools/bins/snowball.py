#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Importance sampling with the snowball proposition distribution in O(k_max x T).

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import networkx as nx
import numpy as np


# Global status dict
status = {"unexplored": 0,
          "boundary": 1,
          "explored": 2}


class History(object):
    """Edge history vectors in MultiGraphs."""

    def __init__(self, X=None):
        """Constructor with automagic."""
        if X is not None:
            self.X = X
            for e in self.X:
                self.X = self.order_edge(e)
        else:
            self.X = []
        self.idx = 0

    def __str__(self):
        return self.X.__str__()

    def __repr__(self):
        return self.X.__repr__()

    def __iter__(self):
        return self

    def __next__(self):
        self.idx += 1
        try:
            return self.X[self.idx - 1]
        except IndexError:
            self.idx = 0
            raise StopIteration  # Done iterating.
    next = __next__  # python2.x compatibility.

    def order_edge(self, edge):
        """Canonical ordering of an edge's representation."""
        if edge[0] > edge[1]:
            return (edge[1], edge[0], edge[2])
        else:
            return edge

    # Public methods
    def load_history(self, filepath, is_simple=False):
        """Load history from file.

        Note
        ----
        The file must contain one edge per line, represented by a pair of
        integers, with an additionnal tag to count multiedges if is a
        multigraphs. The ordering is assumed to represent the history's order.

        Parameters
        ----------
        filepath : str
            Path to history file.
        is_simple : bool
            Assume that the graph is simple and ignore edge tags.
        """
        with open(filepath, 'r') as f:
            self.X = []
            for line in f:
                e = tuple(int(v) for v in line.strip().split())
                if is_simple:
                    self.X.append(self.order_edge((e[0], e[1], 0)))
                else:
                    self.X.append(self.order_edge(e))
            return self.X

    def save_history(self, filepath):
        """Save history to file.

        Note
        ----
        The file must contain one edge per line, represented by a pair of
        integers, with an additionnal tag to count multiedges if is a
        multigraphs. The ordering is assumed to represent the history's order.

        Parameters
        ----------
        filepath : str
            Path to history file.
        is_simple : bool
            Assume that the graph is simple and ignore edge tags.
        """
        with open(filepath, 'r') as f:
            for e in self.X:
                print(list(e).join(" "), file=f)

    def set_size(self, T):
        """Set total size of history when used in `set` mode."""
        self.X = [(0, 0, 0)] * T

    def set(self, t, e):
        """Set edge at time t."""
        self.X[t] = self.order_edge(e)

    def clear(self):
        """Clear container for reuse."""
        self.X.clear()

    def tau(self, e):
        """Postition of edge 'e' in the history."""
        if len(e) == 2:
            e = (e[0], e[1], 0)
        return self.X.index(self.order_edge(e))


class Sampler(object):
    """Baes Sampler class."""

    def __init__(self, g):
        """Construct sampler for a graph."""
        self.N = g.number_of_nodes()
        self.g = g
        self.reset()
        # super(Sampler, self).__init__()

    def reset(self):
        """Reset all sample dependents values."""
        self.deg = [0] * self.N
        self.Z = 0

    def loglikelihood_update(self, edge, gamma, b):
        """Compute the change in log-likelihood due to adding `edge`."""
        v0 = edge[0]
        v1 = edge[1]
        if self.deg[v0] == 0 and self.deg[v1] == 0:
            return 0
            # raise Exception("Impossible event!")
        if self.deg[v0] == 0:    # new node is v0
            delta = np.log(b * (self.deg[v1] ** gamma / self.Z))
        elif self.deg[v1] == 0:  # new node is v1
            delta = np.log(b * (self.deg[v0] ** gamma / self.Z))
        else:
            delta = np.log((1 - b)) +\
                np.log((self.deg[v0] * self.deg[v1]) ** gamma / (self.Z ** 2))
        return delta

    def order_edge(self, edge):
        """Canonical ordering of an edge's representation."""
        if edge[0] > edge[1]:
            return (edge[1], edge[0], edge[2])
        else:
            return edge

    def update_degrees(self, edge, gamma, b):
        """Update degree vector and norm upon adding `edge`."""
        v0 = edge[0]
        v1 = edge[1]
        if self.deg[v0] == 0 and self.deg[v1] == 0:
            pass
        elif self.deg[v0] == 0:    # new node is v0
            self.Z -= self.deg[v1] ** gamma
            self.Z += (self.deg[v1] + 1) ** gamma + 1
        elif self.deg[v1] == 0:  # new node is v1
            self.Z -= self.deg[v0] ** gamma
            self.Z += (self.deg[v0] + 1) ** gamma + 1
        else:
            if v1 == v0:  # count self-loops
                self.Z -= self.deg[v0] ** gamma
                self.Z += (self.deg[v0] + 2) ** gamma
            else:
                self.Z -= self.deg[v0] ** gamma +\
                    self.deg[v1] ** gamma
                self.Z += (self.deg[v0] + 1) ** gamma +\
                    (self.deg[v1] + 1) ** gamma
        self.deg[edge[0]] += 1
        self.deg[edge[1]] += 1

    def get_neighboring_edges(self, edge):
        """Get edges neighboring `edge`, excluding itself."""
        x = [(edge[0], v, _) for v in nx.neighbors(self.g, edge[0]) for _
             in range(self.g.number_of_edges(edge[0], v)) if v != edge[1]] +\
            [(edge[0], edge[1], tag) for tag in
             range(self.g.number_of_edges(edge[0], edge[1])) if tag != edge[2]]
        if edge[0] != edge[1]:  # not a self-loop
            x += [(edge[1], v, _) for v in nx.neighbors(self.g, edge[1]) for _
                  in range(self.g.number_of_edges(edge[1], v)) if v != edge[0]]
        return set([self.order_edge(e) for e in x])

    def get(self, gamma, b, alpha=0):
        """Get a sample and its weight."""
        # Reset everything
        self.reset()
        X = History()
        X.set_size(self.g.number_of_edges())
        nx.set_edge_attributes(self.g, status["unexplored"], 'status')
        boundary = [0] * self.g.number_of_edges()
        self.Z = 2
        log_Q = -np.log(self.g.number_of_edges())
        log_P = 0
        t = 0

        # Draw initial edge
        if np.isclose(alpha, 0):
            weights = np.zeros(self.g.number_of_edges() -
                               self.g.number_of_selfloops())
            idx = -1
            for e in self.g.edges:
                if e[0] != e[1]:
                    idx += 1
                    weights[idx] = (self.g.degree(e[0]) * self.g.degree(e[1])) ** alpha
            id_init = np.random.choice(range(len(weights)), p=weights/np.sum(weights))
        else:
            id_init = np.random.randint(0, self.g.number_of_edges() -
                                        self.g.number_of_selfloops())
        idx = -1
        for e in self.g.edges:
            if e[0] != e[1]:
                idx += 1
            if idx == id_init:
                init_edge = self.order_edge(e)
                break

        # X.(0, edge)
        boundary_size = 1
        boundary[0] = init_edge
        self.g[init_edge[0]][init_edge[1]][init_edge[2]]["status"] = status["boundary"]

        # Iterate
        while boundary_size > 0:
            # draw an edge from the boundary [O(1)]
            id_edge = np.random.randint(0, boundary_size)
            edge = boundary[id_edge]
            # add to history and compute the change in likelihoods
            X.set(t, edge)
            log_P += self.loglikelihood_update(edge, gamma, b)
            log_Q -= np.log(boundary_size)
            self.update_degrees(edge, gamma, b)
            # delete from boundary
            self.g[edge[0]][edge[1]][edge[2]]["status"] = status["explored"]
            boundary[id_edge] = boundary[boundary_size - 1]
            boundary[boundary_size - 1] = 0
            boundary_size -= 1
            # update boundary  [O(k_max)]
            neighbors = self.get_neighboring_edges(edge)
            for e in neighbors:
                if self.g[e[0]][e[1]][e[2]]["status"] == status["unexplored"]:
                    self.g[e[0]][e[1]][e[2]]["status"] = status["boundary"]
                    boundary[boundary_size] = e
                    boundary_size += 1
            t += 1
        return (X, log_P, log_Q)



def update(g, X, delta_P, delta_Q):
    """Update estimators on the edges."""
    for t, e in enumerate(X):
        g[e[0]][e[1]][e[2]]["tau"] += t * np.exp(delta_P - delta_Q)


if __name__ == '__main__':
    if not float(nx.__version__) >= 2.0:
        print("Warning! Designed for networkx 2.0")
    import argparse as ap
    prs = ap.ArgumentParser(prog="snowball",
                            description="Snowball sampler for generalized PA.")
    prs.add_argument('--gamma', '-g', type=float, default=0,
                     help='Gamma, exponent of the attachment kernel.')
    prs.add_argument('--b', '-b', type=float, default=1,
                     help='b, node creation probability.')
    prs.add_argument('--alpha', '-a', type=float, default=0,
                     help='Seed bias.')
    prs.add_argument('--num_samples', '-n', type=float, default=100,
                     help='Number of samples.')
    prs.add_argument("--dist_root", '-r', action='store_true',
                     help="Output root probability distribution.")
    prs.add_argument("--find_rood", '-f', action='store_true',
                     help="Perform root-finding.")
    prs.add_argument("--complete", '-c', type=float, nargs='*',
                     help="Perform full inference with outputs at the specified times.")
    prs.add_argument("--order", '-o', type=int, nargs=2,
                     help="Find order of the edges in the position N and M \
                           of the edge list, starting from 1.")
    prs.add_argument("edge_list", type=str, help='Path to an edge list.')
    args = prs.parse_args()

    # Read graph
    g = nx.read_edgelist(args.edge_list,
                         nodetype=int,
                         create_using=nx.MultiGraph(),
                         data=[("key", int)])
    edge_to_idx = dict()
    idx_to_edge = dict()
    for idx, e in enumerate(g.edges):
        edge_to_idx[e] = idx
        idx_to_edge[idx] = e

    # Declare sampler
    sampler = Sampler(g)

    # Calibration step
    ns = min((1000, args.num_samples))
    log_P = np.zeros(ns)
    log_Q = np.zeros(ns)
    samples = [0] * ns
    for i in range(ns):
        (samples[i], log_P[i], log_Q[i]) = sampler.get(args.gamma, args.b)
    avg_log_P = np.mean(log_P)
    avg_log_Q = np.mean(log_Q)

    if args.complete:
        # init
        nx.set_edge_attributes(g, 0, "tau")
        # Re-use calibration histories
        for i in range(ns):
            update(g, samples[i], log_P[i] - avg_log_P, log_Q[i] - avg_log_Q)

        # Sample from the process
        for i, k in enumerate([ns] + args.complete + [args.num_samples]):
            if k == args.num_samples:
                break
            elif k != args.complete[-1]:
                kp1 = args.complete[i]
            else:
                kp1 = args.num_samples
            for i in range(int(k), int(kp1)):
                (sample, log_P, log_Q) = sampler.get(args.gamma, args.b)
                update(g, sample, log_P - avg_log_P, log_Q - avg_log_Q)
            # Normalize and output
            actual_norm = sum(nx.get_edge_attributes(g, "tau").values())
            target_norm = g.number_of_edges() * (g.number_of_edges() - 1) / 2
            # Y  = [0] * g.number_of_edges()
            print("# t", int(kp1))
            for t, e in enumerate(nx.get_edge_attributes(g, "tau")):
                print(e[0], e[1], g[e[0]][e[1]][e[2]]["tau"] / actual_norm * target_norm, flush=True)


    if args.find_rood or args.dist_root:
        # init
        root_distribution = np.zeros(g.number_of_edges())
        # Re-use calibration histories
        for i in range(ns):
            root_distribution[edge_to_idx[samples[i].X[0]]] += np.exp(log_P[i] - avg_log_P - log_Q[i] + avg_log_Q)

        # Sample from the process
        for i in range(ns, int(args.num_samples)):
            (sample, log_P, log_Q) = sampler.get(args.gamma, args.b)
            root_distribution[edge_to_idx[sample.X[0]]] += np.exp(log_P - avg_log_P - log_Q + avg_log_Q)
        
        # output
        root_distribution = root_distribution / np.sum(root_distribution)
        if args.dist_root:
            for e in edge_to_idx:
                print(e, root_distribution[edge_to_idx[e]])
        if args.find_rood:
            max_id = np.argmax(root_distribution)
            print(idx_to_edge[max_id], root_distribution[max_id])
