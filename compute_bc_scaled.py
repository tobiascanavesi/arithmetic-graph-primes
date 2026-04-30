#!/usr/bin/env python3
"""
Scaled-up betweenness centrality computation for the arithmetic graph G_N.

Uses NetworKit (C++ backend) for exact Brandes BC at N up to 100k.
Replaces the sampled-BC at N=10k in the original paper.

Outputs:
  bc_results.json  -- per-scale summary stats
  bc_topnodes.json -- top-50 BC nodes at each scale (with classification)
  bc_highway.json  -- highway BC distribution at each scale
"""

import json
import time
from pathlib import Path

import networkit as nk
import numpy as np
from scipy.stats import pearsonr
from sympy import factorint, isprime, primerange


# ── Graph construction (NetworKit) ──────────────────────────────────

def build_graph_nk(N):
    """Build arithmetic graph G_N using NetworKit."""
    # NetworKit graphs are 0-indexed internally; we use nodes 0..N
    # but skip node 0 (so vertex i corresponds to integer i)
    G = nk.Graph(N + 1, weighted=False, directed=False)
    # additive edges
    for i in range(1, N):
        G.addEdge(i, i + 1)
    # multiplicative edges
    for p in primerange(2, N + 1):
        max_k = N // p
        for k in range(1, max_k):
            G.addEdge(k * p, (k + 1) * p)
    return G


# ── Node classification ─────────────────────────────────────────────

def node_omega(v):
    return len(factorint(v)) if v > 1 else 0


def node_type(v):
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


# ── Per-scale analysis ──────────────────────────────────────────────

def analyze(N):
    print(f"\n=== N = {N} ===", flush=True)
    t0 = time.time()
    G = build_graph_nk(N)
    n_edges = G.numberOfEdges()
    print(f"  built graph: {N} nodes, {n_edges} edges in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    # Brandes algorithm - exact, O(VE)
    bc_alg = nk.centrality.Betweenness(G, normalized=True)
    bc_alg.run()
    bc_scores = bc_alg.scores()  # list, indexed by node id
    print(f"  computed exact BC in {time.time()-t0:.1f}s", flush=True)

    # Skip node 0 (placeholder); use nodes 1..N
    nodes = list(range(1, N + 1))
    degrees = [G.degree(v) for v in nodes]
    bcs = [bc_scores[v] for v in nodes]

    rho, _ = pearsonr(degrees, bcs)

    # Top-50 BC
    bc_pairs = sorted(zip(nodes, bcs), key=lambda x: -x[1])
    top50 = bc_pairs[:50]
    top15 = bc_pairs[:15]
    top10 = bc_pairs[:10]

    n_semi = sum(1 for v, _ in top50 if node_type(v) in ("2p", "3p", "prime"))
    n_2p_3p_50 = sum(1 for v, _ in top50 if node_type(v) in ("2p", "3p"))
    n_2p = sum(1 for v, _ in top15 if node_type(v) == "2p")
    n_3p = sum(1 for v, _ in top15 if node_type(v) == "3p")

    # Top-10 detail
    top10_detail = []
    for v, b in top10:
        f = factorint(v)
        fstr = " * ".join(str(p) if e == 1 else f"{p}^{e}" for p, e in sorted(f.items()))
        top10_detail.append({
            "v": int(v),
            "factorization": fstr,
            "omega": node_omega(v),
            "deg": int(G.degree(v)),
            "bc": float(b),
        })

    # Highway BC at large primes (p in (N/4, N/3])
    highway = []
    for p in primerange(N // 4 + 1, N // 3 + 1):
        if 3 * p > N:
            continue
        bcp = float(bc_scores[p])
        bc2p = float(bc_scores[2 * p])
        bc3p = float(bc_scores[3 * p])
        total = bcp + bc2p + bc3p
        share2p = bc2p / total if total > 0 else 0.0
        highway.append({
            "p": int(p),
            "bc_p": bcp,
            "bc_2p": bc2p,
            "bc_3p": bc3p,
            "share_2p": share2p,
        })
        if len(highway) >= 8:
            break

    summary = {
        "N": N,
        "edges": int(n_edges),
        "rho": float(rho),
        "top50_semiprime_or_prime_pct": int(n_semi * 2),
        "top50_2p_3p_count": int(n_2p_3p_50),
        "top15_2p_count": int(n_2p),
        "top15_3p_count": int(n_3p),
        "max_degree": int(max(degrees)),
        "max_bc": float(max(bcs)),
        "top10": top10_detail,
        "highway_sample": highway,
    }

    print(f"  rho(deg, BC) = {rho:.4f}", flush=True)
    print(f"  top-50 semi/prime: {n_semi*2}%", flush=True)
    print(f"  top-15 (2p, 3p): ({n_2p}, {n_3p})", flush=True)
    print(f"  highway 2p mean share: {np.mean([h['share_2p'] for h in highway])*100:.1f}%", flush=True)
    print(f"  top-3 BC: {[f'{v}={b:.5f}' for v, b in top10[:3]]}", flush=True)

    return summary


def main():
    out_dir = Path(__file__).parent
    scales = [500, 1000, 2000, 5000, 10000, 20000, 50000]
    # Add 100000 if time allows; user can re-run with extended scales.

    all_results = []
    for N in scales:
        res = analyze(N)
        all_results.append(res)
        # incremental save
        with open(out_dir / "bc_results.json", "w") as f:
            json.dump(all_results, f, indent=2)

    print("\n=== Summary table ===")
    print(f"{'N':>8} {'edges':>9} {'rho':>8} {'top50%':>7} {'15:2p/3p':>10}")
    for r in all_results:
        print(f"{r['N']:>8} {r['edges']:>9} {r['rho']:>8.4f} "
              f"{r['top50_semiprime_or_prime_pct']:>6}% "
              f"{r['top15_2p_count']:>3}/{r['top15_3p_count']:<2}")


if __name__ == "__main__":
    main()
