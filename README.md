# Semiprime Bottlenecks in Arithmetic Graphs

Code and data for the paper *"Semiprime Bottlenecks in Arithmetic Graphs: Betweenness Centrality from Prime Highways"* by Tobias Canavesi.

## What is an arithmetic graph?

The arithmetic graph G_N has vertices {1, ..., N} with two kinds of edges:
- **Additive**: n → n+1 (backbone chain)
- **Multiplicative**: kp → (k+1)p for each prime p (prime highways)

The main finding: semiprimes 2p and 3p (p a large prime) dominate betweenness centrality, not high-degree primorials.

## Reproduce

```bash
pip install networkx sympy numpy matplotlib scipy
python reproduce.py
```

This prints all paper tables to stdout and saves three PDF figures:
- `fig_degree_vs_bc.pdf` — degree vs betweenness scatter
- `fig_correlation_decay.pdf` — correlation decay with N
- `fig_ushape.pdf` — routing cost U-shape

Runtime: ~10 minutes (mostly betweenness at N=10,000).

## Files

| File | Description |
|------|-------------|
| `paper.tex` | LaTeX source |
| `reproduce.py` | Single script reproducing all results |
| `fig_*.pdf` | Generated figures |

## License

MIT
