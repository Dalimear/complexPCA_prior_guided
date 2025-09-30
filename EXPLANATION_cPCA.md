# Complex PCA and Guided cPCA: Outputs and Interpretation (GitHub-safe math)

## Summary
Complex PCA (cPCA) is applied to analytic (complex) BOLD signals to obtain complex components (time courses and spatial loadings with phase). The reconstruction produces a matrix with dimensions (bins x parcels), where each entry represents a z-scored BOLD level at a given phase bin of one oscillatory cycle. Guided cPCA augments the subject covariance with a low-rank prior aligned to group templates, which encourages but does not force alignment.

## 1) Data and normalization
Let $X \in \mathbb{C}^{T \times P}$ denote timepoints by parcels (or voxels) after standard preprocessing. When `--normalize zscore` is used, each parcel is standardized:
$$
z_p(t) = \frac{x_p(t) - \overline{x}_p}{\operatorname{sd}(x_p)}.
$$
Analytic (complex) signals are then formed via the Hilbert transform, and cPCA is applied to $X$.

## 2) cPCA via SVD
The singular value decomposition is
$$
X = U\,\Sigma\,V^{H}.
$$
Here, $H$ denotes the conjugate transpose.

Definitions:
- Component time courses (complex scores): $s_k(t) = (U\Sigma)_{t,k}$.
- Spatial loadings (complex): $\ell_{k,p}$ from columns of $V$ (with the script’s scaling).
- Explained variance (eigenvalues of the covariance): $\lambda_k = \sigma_k^2/(T-1)$, real and non-negative.

Interpretation:
- $|\ell_{k,p}|$ indicates the participation strength of parcel $p$ in component $k$.
- $\arg(\ell_{k,p})$ indicates the phase lag or lead of parcel $p$ relative to the component time course.
- $\lambda_k$ quantifies the total variance captured by component $k$; it is not written into the reconstruction matrix directly.

## 3) Reconstruction and phase binning
For a component $k$, the contribution to parcel $p$ at time $t$ is, up to consistent scaling,
$$
\widehat{x}_k(t,p) = \Re\{\, s_k(t)\,\ell_{k,p} \,\}.
$$

The phase of the component is $\phi_k(t) = \arg s_k(t)$. The interval $[-\pi, \pi)$ is divided into $N$ bins (e.g., $N=30$), and values are averaged within each bin to form
$$
M[b,p] = \mathbb{E}\big[\,\widehat{x}_k(t,p)\,\big|\, \phi_k(t)\in \text{bin } b\big].
$$

Shapes and units:
- $M$ has shape (bins x parcels); for BNA-246 and $N=30$, this is $30 \times 246$.
- Entries of $M$ are in z-units (standard deviations relative to each parcel’s mean) because inputs were z-scored.
- Rows correspond to phase snapshots across one cycle; columns correspond to parcels.
- A parcel’s reconstructed waveform across the cycle is a column of $M$.
- These values are phase-conditioned mean amplitudes, not correlations and not eigenvalues.

## 4) Mapping bins to seconds
The $N$ bins represent one cycle of the component. If the dominant frequency is $f$ Hz, the period is $1/f$ seconds, and
$$
\text{seconds per bin} = \frac{1}{f\,N}.
$$
The dominant frequency may be estimated from the component time course (e.g., via an FFT peak). If TR is uncertain, a cycle length can first be estimated in volumes and then converted to seconds once TR is established.

## 5) Guided cPCA
Let $X_c$ denote column-centered data and
$$
C = \frac{X_c^{H} X_c}{T-1}
$$
denote the subject covariance. Let $T \in \mathbb{C}^{P \times K}$ stack unit-normalized group template maps as columns. Guided cPCA forms
$$
C_{\text{guided}} = C + T\,\mathrm{diag}(\lambda)\,T^{H},
$$
where $\lambda \in \mathbb{R}_{\geq 0}^{K}$ contains non-negative guidance weights (variance units).

Term meanings:
- $C$: subject covariance that would be diagonalized in conventional PCA/cPCA.
- $T$: group template directions.
- $\mathrm{diag}(\lambda)$: guidance weights, one per template direction.
- $C_{\text{guided}}$: covariance with a rank-$K$ prior added along the template directions.

An eigendecomposition of $C_{\text{guided}}$ yields guided components. The script selects the smallest $\lambda$ (from a grid scaled to data variance) that achieves a target template–component similarity, which encourages but does not force alignment. Guided components are matched and phase-aligned to templates, and subject time courses are obtained by projection, $S = X_c\,U$.

Notes on $\lambda$:
- $\lambda$ is not the eigenvalue of the group component; it is a guidance weight controlling the strength of the prior.
- As $\lambda \to 0$, ordinary cPCA is recovered; larger $\lambda$ increases alignment to the template.

## 6) Relation to the global signal
This workflow can be run without global signal regression. cPCA typically isolates an in-phase, spatially widespread mode as one component and a phase-opposed, propagating pattern (QPP-like) as another. Apparent “anticorrelations” in zero-lag functional connectivity are represented here as phase opposition within a spatiotemporal mode, not as a by-product of subtracting the global mean.

## 7) Common interpretation questions
- Are reconstruction values correlations? No. They are z-scored, phase-conditioned means.
- Where do eigenvalues appear? As explained variances per component; they are not used as individual cells in the reconstruction matrix.
- What sets the sign or polarity across bins? The relative phase between $\ell_{k,p}$ and $s_k(t)$; sign reversals across bins reflect phase progression.
- Can the bin count be changed? Yes, via `--n_recon_bins`, trading temporal detail for signal-to-noise.
- Is double bandpass filtering recommended? If data have already been filtered upstream, an additional bandpass step is generally avoided to prevent unintended narrowing of the passband and phase distortion.

## One-sentence summary
cPCA is applied to analytic BOLD; a component’s instantaneous contribution is $\Re\{s_k(t)\ell_{k,\cdot}\}$, which is then phase-binned into a bins x parcels, z-scored template of one cycle; guided cPCA adds a low-rank prior $T\,\mathrm{diag}(\lambda)\,T^{H}$ so components remain data-driven while being gently aligned to a group template.
