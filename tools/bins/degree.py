#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classify edges by degree.

Give the age of the oldest nodes to each edge (assuming high degree
represents age).

Author: Jean-Gabriel Young <info@jgyoung.ca>
"""
import sys
import networkx as nx

if __name__ == '__main__':
    if not float(nx.__version__) >= 2.0:
        print("Warning! Designed for networkx 2.0")

    # Read graph as a MultiGraph.
    g = nx.read_edgelist(sys.argv[1],
                         nodetype=int,
                         create_using=nx.MultiGraph(),
                         data=[("key", int)])

    # Combine edges with the reverse sorted pair of the degree of their nodes,
    # then sort based on degrees pairs, breaking ties with the smallest degree.
    data = sorted([(e, (g.degree(e[0]), g.degree(e[1])))
                    if g.degree(e[0]) > g.degree(e[1])
                    else (e, (g.degree(e[1]),g.degree(e[0])))
                    for e in g.edges],
                  key=lambda x: x[1], reverse=True)

    T = nx.number_of_edges(g)
    avg = (T - 1) / 2

    # Output ordering, with an average rank for tied edges.
    class_content = []
    class_type = data[0][1]  # degree-degree pair
    class_base_time = 0
    for t, x in enumerate(data):
        if x[1] != class_type:
            for c in class_content:
                print(*c, class_base_time + (len(class_content) + 1) / 2 - 1)
            class_content.clear()
            class_type = x[1]
            class_base_time = t
        class_content.append((x[0][1], x[0][0], x[0][2])
                             if x[0][1] < x[0][0] else
                             (x[0][0], x[0][1], x[0][2]))

    # Output last class
    for c in class_content:
        print(*c, class_base_time + (len(class_content) + 1) / 2 - 1)
