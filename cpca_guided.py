
#!/usr/bin/env python3
"""
cpca_guided.py — Command-line Guided CPCA for complex fMRI time series (parcellated 2D text).

Usage (mirrors your existing script style):

    python cpca_guided.py \
        -i subject_txt_list.txt \
        -n 3 \
        --templates group_templates.txt \
        --target_sim 0.6 \
        --out subj1_guided_cpca

Where:
- subject_txt_list.txt is a text file where each line is a path to a 2D txt (T x V) or .npy file.
- --templates is a (V x K) text or .npy file with group components as columns.
- -n is the number of components to return (usually K).

Outputs:
- <out>.npz with components (V x K, complex), timecourses (T_total x K, complex), sims, lambdas, evals.
- Optional: exports to plain text via --write_txt.

Author: ChatGPT
"""

import argparse
import numpy as np
import os
from numpy.linalg import norm
from scipy.signal import hilbert

# ---- Guided CPCA core (lightweight) ----

def _center_time(X):
    return X - X.mean(axis=0, keepdims=True)

def _unit_norm_cols(M, eps=1e-12):
    norms = np.linalg.norm(M, axis=0)
    scale = np.where(norms > eps, 1.0 / norms, 0.0)
    return M * scale

def _subject_covariance(Xc):
    T = Xc.shape[0]
    return (Xc.conj().T @ Xc) / max(T - 1, 1)

def _largest_eigval_hermitian(C):
    w = np.linalg.eigvalsh(C)
    return float(w[-1].real)

def _guided_covariance(C, Tmpl, lambdas):
    return C + (Tmpl * lambdas[np.newaxis, :]) @ Tmpl.conj().T

def _top_k_eig(C, k):
    w, V = np.linalg.eigh(C)  # ascending
    idx = np.argsort(w)[::-1][:k]
    return w[idx].real, V[:, idx]

def _match_components(U, Tmpl):
    K = Tmpl.shape[1]
    S = np.abs(Tmpl.conj().T @ U)
    order = np.full(K, -1, dtype=int)
    sims = np.zeros(K, dtype=float)
    taken = np.zeros(U.shape[1], dtype=bool)
    for k in range(K):
        s = S[k].copy()
        s[taken] = -np.inf
        j = int(np.argmax(s))
        order[k] = j
        sims[k] = float(S[k, j])
        taken[j] = True
    return order, sims

def _phase_align(U, Tmpl, order):
    U_aligned = U.copy()
    for k, j in enumerate(order):
        z = Tmpl[:, k].conj().T @ U[:, j]
        phi = -np.angle(z)
        U_aligned[:, j] = U[:, j] * np.exp(1j * phi)
    return U_aligned

def guided_cpca(X, templates, n_components=None,
                lambdas=None, target_sim=None,
                lambda_grid=(0.0, 0.1, 0.5, 1.0, 2.0),
                scale_by_leading=True, ridge=0.0):
    X = np.asarray(X)
    Tmpl = np.asarray(templates)
    T, V = X.shape
    Vt, K = Tmpl.shape
    if n_components is None:
        n_components = K
    assert V == Vt, f"Feature mismatch: X has V={V}, templates have V={Vt}"
    assert 1 <= n_components <= V

    # center and covariance
    Xc = _center_time(X)
    C = _subject_covariance(Xc).astype(np.complex128, copy=False)

    # normalize templates
    Tmpl = _unit_norm_cols(Tmpl.astype(np.complex128, copy=False))

    if lambdas is not None:
        lambdas = np.asarray(lambdas, dtype=float)
        assert lambdas.shape == (K,)
        lambda_used = lambdas.copy()
    else:
        L1 = _largest_eigval_hermitian(C) if scale_by_leading else 1.0
        best_s = None
        achieved = None
        for base in lambda_grid:
            s = float(base) * float(L1)
            lambdas_try = np.full(K, s, dtype=float)
            Ct = _guided_covariance(C, Tmpl, lambdas_try)
            if ridge > 0:
                Ct = Ct + (ridge * np.eye(V))
            evals, U = _top_k_eig(Ct, n_components)
            U = U / np.maximum(np.linalg.norm(U, axis=0, keepdims=True), 1e-12)
            order, sims = _match_components(U, Tmpl[:, :n_components])
            med = np.median(sims)
            if target_sim is None:
                best_s = s; achieved = med; break
            if med >= target_sim and best_s is None:
                best_s = s; achieved = med; break
            achieved = med
        if best_s is None:
            best_s = float(lambda_grid[-1]) * float(L1)
        lambda_used = np.full(K, best_s, dtype=float)

    Ct = _guided_covariance(C, Tmpl, lambda_used)
    if ridge > 0:
        Ct = Ct + (ridge * np.eye(V))
    evals, U = _top_k_eig(Ct, n_components)
    U = U / np.maximum(np.linalg.norm(U, axis=0, keepdims=True), 1e-12)
    order, sims = _match_components(U, Tmpl[:, :n_components])
    U = _phase_align(U, Tmpl[:, :n_components], order)
    U = U[:, order]
    S = Xc @ U
    return {
        "components_": U,
        "timecourses_": S,
        "evals_": evals,
        "sims_": sims,
        "lambdas_": lambda_used
    }

# ---- I/O helpers for list-of-files workflow ----

def _load_matrix(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npy":
        A = np.load(path)
    else:
        # default: text file with whitespace-separated numbers
        A = np.loadtxt(path)
    if A.ndim != 2:
        raise ValueError(f"{path}: expected 2D matrix, got shape {A.shape}")
    return A

def load_listfile(listfile):
    with open(listfile, "r") as f:
        paths = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    if len(paths) == 0:
        raise ValueError("Input list is empty.")
    mats = [_load_matrix(p) for p in paths]
    # ensure all have same number of columns (V)
    Vset = {M.shape[1] for M in mats}
    if len(Vset) != 1:
        raise ValueError(f"All matrices must have same V (columns). Got {Vset}")
    X = np.vstack(mats)  # concatenate along time
    return X, paths

def load_templates(path):
    ext = os.path.splitext(path)[1].lower()
    Tmpl = np.load(path) if ext == ".npy" else np.loadtxt(path)
    Tmpl = np.asarray(Tmpl)
    if Tmpl.ndim == 1:
        Tmpl = Tmpl[:, None]      # make (V,) -> (V,1)
    if Tmpl.ndim != 2:
        raise ValueError(f"{path}: expected 2D (V x K), got {Tmpl.shape}")
    return Tmpl


def maybe_hilbert_complex(X, use_complex):
    if not use_complex:
        return X
    # Compute analytic signal per column and return complex-conjugated (as in your repo)
    Xh = hilbert(X, axis=0).conj()
    return Xh

def write_txt_or_npz(base, out, write_txt=False):
    comps = out["components_"]
    tcs = out["timecourses_"]
    sims = out["sims_"]
    lambdas = out["lambdas_"]
    evals = out["evals_"]
    np.savez(base + ".npz",
             components=comps, timecourses=tcs,
             sims=sims, lambdas=lambdas, evals=evals)
    if write_txt:
        np.savetxt(base + "_components_real.txt", comps.real)
        np.savetxt(base + "_components_imag.txt", comps.imag)
        np.savetxt(base + "_timecourses_real.txt", tcs.real)
        np.savetxt(base + "_timecourses_imag.txt", tcs.imag)
        np.savetxt(base + "_sims.txt", sims[None, :])
        np.savetxt(base + "_lambdas.txt", lambdas[None, :])
        np.savetxt(base + "_evals.txt", evals[None, :])

# ---- CLI ----

def main():
    ap = argparse.ArgumentParser(description="Guided CPCA with template prior (txt or npy inputs).")
    ap.add_argument("-i", "--input", required=True,
                    help="Path to a text file listing subject matrices (.txt or .npy), one per line (rows=time, cols=parcels).")
    ap.add_argument("-n", "--n_comps", required=True, type=int,
                    help="Number of components to return (usually equals # template columns).")
    ap.add_argument("--templates", required=True,
                    help="Path to group templates (V x K) as .txt or .npy; columns are components.")
    ap.add_argument("--complex", action="store_true",
                    help="Apply Hilbert transform (analytic signal) per column and conjugate (like your CPCA).")
    ap.add_argument("--target_sim", type=float, default=0.6,
                    help="Target median |<t_k, u_k>| for auto lambda selection.")
    ap.add_argument("--lambda_grid", default="0,0.1,0.5,1,2",
                    help="Comma-separated grid (scaled by leading eigval of covariance).")
    ap.add_argument("--ridge", type=float, default=1e-8,
                    help="Optional diagonal ridge added to guided covariance.")
    ap.add_argument("--out", default=None,
                    help="Output prefix (default: derived from first input).")
    ap.add_argument("--write_txt", action="store_true",
                    help="Also export real/imag parts and metadata as .txt alongside the .npz.")
    args = ap.parse_args()

    X, paths = load_listfile(args.input)
    Tmpl = load_templates(args.templates)

    if args.n_comps > Tmpl.shape[1]:
        raise ValueError(f"n_comps={args.n_comps} but templates have only K={Tmpl.shape[1]} columns")

    Xcplx = maybe_hilbert_complex(X, args.complex)

    grid = tuple(float(x) for x in args.lambda_grid.split(","))
    out = guided_cpca(
        Xcplx, Tmpl,
        n_components=args.n_comps,
        lambdas=None,
        target_sim=args.target_sim,
        lambda_grid=grid,
        scale_by_leading=True,
        ridge=args.ridge
    )

    base = args.out if args.out else os.path.splitext(os.path.basename(paths[0]))[0] + "_guided_cpca"
    write_txt_or_npz(base, out, write_txt=args.write_txt)
    print(f"[Guided CPCA] Saved: {base}.npz")
    if args.write_txt:
        print(f"[Guided CPCA] Also wrote text exports with prefix: {base}_*")

if __name__ == "__main__":
    main()
