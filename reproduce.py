#!/usr/bin/env python3
"""
Reproduce all tables and figures from:
  "Semiprime Bottlenecks in Arithmetic Graphs"
  by Tobias Canavesi

Usage:  python reproduce.py
Output: Tables 1-6 printed to stdout; three PDF figures saved to cwd.
Requires: networkx, sympy, numpy, matplotlib, scipy
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange, factorint, isprime
from scipy.stats import pearsonr


# ── Graph construction ──────────────────────────────────────────────

def build_graph(N):
    """Build arithmetic graph G_N."""
    G = nx.Graph()
    G.add_nodes_from(range(1, N + 1))
    for i in range(1, N):
        G.add_edge(i, i + 1)
    for p in primerange(2, N + 1):
        for k in range(1, N // p):
            G.add_edge(k * p, (k + 1) * p)
    return G


def node_omega(v):
    """Number of distinct prime factors of v."""
    return len(factorint(v)) if v > 1 else 0


def node_type(v):
    """Classify node: '2p', '3p', 'prime', or 'other'."""
    if v <= 1:
        return "other"
    f = factorint(v)
    keys = sorted(f.keys())
    if len(keys) == 1 and f[keys[0]] == 1:
        return "prime"
    if len(keys) == 2 and all(f[k] == 1 for k in keys):
        if keys[0] == 2:
            return "2p"
        if keys[0] == 3:
            return "3p"
    return "other"


# ── Table 1: Degree verification (Table 1 in paper) ────────────────

def table_degree(G, N):
    print("\n=== Table 1: Degree formula verification ===")
    nodes = [106, 173, 210, 331, 2310]
    print(f"{'v':>6}  {'Factorization':<20} {'w(v)':>4} {'deg':>4} {'Predicted':<16}")
    for v in nodes:
        if v > N:
            continue
        w = node_omega(v)
        d = G.degree[v]
        pred = f"2+2*{w}={2+2*w}" if d == 2 + 2 * w else f"boundary"
        f = factorint(v) if v > 1 else {}
        fstr = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
        if isprime(v):
            fstr = "prime"
        print(f"{v:>6}  {fstr:<20} {w:>4} {d:>4} {pred:<16}")


# ── Table 2: Graph metrics (Table 2 in paper) ──────────────────────

def table_scaling():
    print("\n=== Table 2: Graph metrics at different scales ===")
    print(f"{'N':>6} {'Edges':>8} {'Density':>8} {'Diam':>6} {'AvgPath':>8}")
    for N in [100, 500, 1000, 2000, 5000]:
        G = build_graph(N)
        e = G.number_of_edges()
        d = nx.density(G)
        diam = nx.diameter(G)
        apl = nx.average_shortest_path_length(G)
        print(f"{N:>6} {e:>8} {d:>8.4f} {diam:>6} {apl:>8.2f}")


# ── Tables 3-4: Betweenness analysis ───────────────────────────────

def betweenness_analysis(scales):
    """Compute BC at multiple scales; return data for figures."""
    results = []
    full_data = {}

    for N in scales:
        G = build_graph(N)
        if N <= 5000:
            bc = nx.betweenness_centrality(G)
        else:
            bc = nx.betweenness_centrality(G, k=1000)

        degrees = [G.degree[v] for v in G.nodes()]
        bcs = [bc[v] for v in G.nodes()]
        rho = pearsonr(degrees, bcs)[0]

        top50 = sorted(bc, key=bc.get, reverse=True)[:50]
        semi_count = sum(1 for v in top50 if node_type(v) in ("2p", "3p", "prime"))
        top15 = top50[:15]
        n2p = sum(1 for v in top15 if node_type(v) == "2p")
        n3p = sum(1 for v in top15 if node_type(v) == "3p")
        method = "exact" if N <= 5000 else "sampled"

        results.append((N, rho, semi_count, n2p, n3p, method))
        full_data[N] = (G, bc, degrees, bcs)

    print("\n=== Table 3: Semiprime dominance across scales ===")
    print(f"{'N':>6} {'rho':>7} {'Top50%':>7} {'2p/3p':>6} {'Method':<8}")
    for N, rho, sc, n2, n3, m in results:
        print(f"{N:>6} {rho:>7.3f} {sc*2:>6}% {n2:>2}/{n3:<2}  {m:<8}")

    # Table 4: Top-10 at largest scale
    N_max = max(scales)
    G, bc, _, _ = full_data[N_max]
    top10 = sorted(bc, key=bc.get, reverse=True)[:10]
    print(f"\n=== Table 4: Top 10 BC nodes at N={N_max} ===")
    print(f"{'v':>6} {'Factorization':<16} {'w':>3} {'deg':>4} {'BC':>8}")
    for v in top10:
        f = factorint(v)
        fstr = " * ".join(str(p) for p in sorted(f.keys()))
        print(f"{v:>6} {fstr:<16} {node_omega(v):>3} {G.degree[v]:>4} {bc[v]:>8.5f}")

    return results, full_data


# ── Table 5: Highway BC distribution ───────────────────────────────

def table_highway(G, bc, N):
    print(f"\n=== Table 5: Highway BC distribution (N={N}) ===")
    primes = [p for p in primerange(N // 4 + 1, N // 3 + 1) if 3 * p <= N][:5]
    print(f"{'p':>6} {'BC(p)':>8} {'BC(2p)':>8} {'BC(3p)':>8} {'BC(2p)/S':>8}")
    for p in primes:
        bcp = bc.get(p, 0)
        bc2p = bc.get(2 * p, 0)
        bc3p = bc.get(3 * p, 0)
        total = bcp + bc2p + bc3p
        share = bc2p / total * 100 if total > 0 else 0
        print(f"{p:>6} {bcp:>8.5f} {bc2p:>8.5f} {bc3p:>8.5f} {share:>7.1f}%")


# ── Table 6: Routing cost U-shape ──────────────────────────────────

def table_ushape(G, N, n_samples=500):
    print(f"\n=== Table 6: Routing cost U-shape (N={N}) ===")
    primes = [p for p in primerange(N // 4 + 1, N // 3 + 1) if 3 * p <= N][:5]
    nodes = list(G.nodes())
    sample = np.random.choice(nodes, min(n_samples, len(nodes)), replace=False)

    print(f"{'p':>6} {'E[d(.,p)]':>10} {'E[d(.,2p)]':>11} {'E[d(.,3p)]':>11}")
    for p in primes:
        dists = {1: [], 2: [], 3: []}
        for s in sample:
            for k in [1, 2, 3]:
                t = k * p
                if t <= N and t != s:
                    dists[k].append(nx.shortest_path_length(G, s, t))
        d1 = np.mean(dists[1]) if dists[1] else float("nan")
        d2 = np.mean(dists[2]) if dists[2] else float("nan")
        d3 = np.mean(dists[3]) if dists[3] else float("nan")
        print(f"{p:>6} {d1:>10.2f} {d2:>11.2f} {d3:>11.2f}")


# ── Figure 1: Degree vs BC scatter ─────────────────────────────────

def fig_degree_vs_bc(G, bc, N, rho):
    nodes = list(G.nodes())
    degs = [G.degree[v] for v in nodes]
    bcs = [bc[v] for v in nodes]
    types = [node_type(v) for v in nodes]

    colors = {"2p": "red", "3p": "orange", "prime": "royalblue", "other": "lightgrey"}
    zorder = {"other": 1, "prime": 2, "3p": 3, "2p": 4}
    labels = {"other": "Other", "prime": r"Prime $p$", "3p": r"Semiprime $3p$", "2p": r"Semiprime $2p$"}

    fig, ax = plt.subplots(figsize=(7, 5))
    for t in ["other", "prime", "3p", "2p"]:
        idx = [i for i, tp in enumerate(types) if tp == t]
        ax.scatter([degs[i] for i in idx], [bcs[i] for i in idx],
                   c=colors[t], s=15, alpha=0.6, label=labels[t], zorder=zorder[t])
    ax.set_xlabel("Degree")
    ax.set_ylabel("Betweenness Centrality")
    ax.set_title(rf"$G_{{{N}}}$: Degree vs Betweenness Centrality ($\rho = {rho:.3f}$)")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig("fig_degree_vs_bc.pdf", dpi=300)
    plt.close(fig)
    print("Saved fig_degree_vs_bc.pdf")


# ── Figure 2: Correlation decay ─────────────────────────────────────

def fig_correlation_decay(results):
    Ns = [r[0] for r in results]
    rhos = [r[1] for r in results]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(Ns, rhos, "ko-", markersize=7)
    for n, r in zip(Ns, rhos):
        ax.annotate(f"{r:.3f}", (n, r), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$\rho(\mathrm{deg}, \mathrm{BC})$")
    ax.set_title("Degree\u2013Betweenness Correlation Decay")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("fig_correlation_decay.pdf", dpi=300)
    plt.close(fig)
    print("Saved fig_correlation_decay.pdf")


# ── Figure 3: Routing cost U-shape ──────────────────────────────────

def fig_ushape(G, N, n_samples=500):
    primes = [p for p in primerange(N // 4 + 1, N // 3 + 1) if 4 * p <= N][:3]
    nodes = list(G.nodes())
    sample = np.random.choice(nodes, min(n_samples, len(nodes)), replace=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    for p in primes:
        ks = [k for k in range(1, 5) if k * p <= N]
        avg_d = []
        for k in ks:
            t = k * p
            dists = [nx.shortest_path_length(G, int(s), t) for s in sample if s != t]
            avg_d.append(np.mean(dists))
        ax.plot(ks, avg_d, "o-", label=f"$p = {p}$", markersize=6)

    ax.set_xlabel("Multiplier $k$")
    ax.set_ylabel("Average distance from random node")
    ax.set_title(f"Routing Cost to Highway Node $kp$ ($N = {N}$)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("fig_ushape.pdf", dpi=300)
    plt.close(fig)
    print("Saved fig_ushape.pdf")


# ── Edge count verification ─────────────────────────────────────────

def verify_edge_count():
    print("\n=== Edge count formula verification ===")
    for N in [100, 500, 1000, 2000, 5000]:
        G = build_graph(N)
        actual = G.number_of_edges()
        formula = (N - 1) + sum(N // p - 1 for p in primerange(2, N + 1))
        match = "OK" if actual == formula else "MISMATCH"
        print(f"  N={N:>5}: actual={actual}, formula={formula}  [{match}]")


# ── Main ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    np.random.seed(42)

    print("=" * 60)
    print("Reproducing: Semiprime Bottlenecks in Arithmetic Graphs")
    print("=" * 60)

    # Edge count verification
    verify_edge_count()

    # Table 1: Degree verification (needs N >= 5000 for node 2310)
    G5k = build_graph(5000)
    table_degree(G5k, 5000)

    # Table 2: Graph metrics
    table_scaling()

    # Tables 3-4: Betweenness at multiple scales
    scales = [500, 1000, 2000, 5000, 10000]
    print("\nComputing betweenness centrality (this takes a few minutes)...")
    results, full_data = betweenness_analysis(scales)

    # Table 5: Highway BC (at N=2000)
    G2k, bc2k, _, _ = full_data[2000]
    table_highway(G2k, bc2k, 2000)

    # Table 6: Routing cost U-shape (at N=2000)
    table_ushape(G2k, 2000)

    # Figures
    G5k, bc5k, _, _ = full_data[5000]
    rho5k = results[3][1]  # rho at N=5000
    fig_degree_vs_bc(G5k, bc5k, 5000, rho5k)
    fig_correlation_decay(results)
    fig_ushape(G2k, 2000)

    print("\nDone. All tables printed; figures saved as PDF.")
