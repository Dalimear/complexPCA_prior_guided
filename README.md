# complexPCA_prior_guided

**Guided Complex PCA for Traveling-Wave Analysis in fMRI Data**

This repository provides a Python implementation of **guided complex PCA (CPCA)** for capturing traveling and standing wave patterns in fMRI data, with an extra twist:  
we guide each subject’s decomposition toward a group template so that components can be reliably compared across subjects.

---

## Background

Complex PCA (CPCA) is an extension of PCA that operates on the **analytic signal** (via the Hilbert transform) of fMRI time series.  
This allows each component to encode not only **amplitude** but also **relative phase (time lag/lead)** across brain regions, making it ideal for identifying traveling waves.

The **guided** variant implemented here introduces a prior:
- Group-level templates (complex spatial maps) are provided.
- Subject decompositions are softly biased toward those templates.
- Components are then matched and phase-aligned to the templates.

This ensures that *Component 2 in subject A* is directly comparable to *Component 2 in subject B*.

---

## Why Guided CPCA?

- Standard CPCA gives components with arbitrary **sign**, **phase**, and **order** across subjects.
- Guided CPCA:
  - Matches subject components to group templates using the **absolute complex inner product** (a similarity score that ignores arbitrary sign/phase flips).
  - Rotates subject components so their inner product with the template is **real and positive** → consistent orientation across subjects.
- Result: **apples-to-apples comparisons** across subjects for amplitude, phase, and time courses.

---

## Usage

Example:

```bash
python cpca_guided.py \
  -i subject_list.txt \
  -n 3 \
  --templates group_templates.npy \
  --complex \
  --target_sim 0.6 \
  --out subj1_guided


Guided Complex PCA (CPCA) — Usage Guide
---------------------------------------

Inputs
------
subject_list.txt
    Text file with paths to subject matrices (rows = time, cols = parcels).

--templates
    Group template file (V × K), either text or .npy.
    Columns correspond to group components.

-n
    Number of components to return (usually equals K).

--complex
    Apply Hilbert transform to build analytic signals.

--target_sim
    Target similarity threshold for lambda tuning (e.g., 0.6).

--out
    Output file prefix (default: derived from input filename).

Outputs
-------
<out>.npz containing:
    components_   Complex spatial maps
    timecourses_  Complex time courses
    sims_         Similarity scores to templates
    evals_        Eigenvalues
    lambdas_      Guidance weights

Optional:
    --write_txt   Export real/imag parts and metadata as plain text.

Requirements
------------
Python 3.x

Dependencies:
    numpy
    scipy

Install via:
    pip install -r requirements.txt

Acknowledgements
----------------
This code is directly inspired by and adapted from the work of Taylor Bolt
(https://github.com/tsb46/complex_pca), developed during his postdoc and subsequent collaboration with the Keilholz MIND lab at Emory/GT (https://sites.google.com/view/keilholz-lab/people?authuser=0).

The guided extension builds on his framework to add template matching and
alignment for subject-to-group comparison.

Notes & Caveats
---------------
- Guidance is soft: tune lambda to control how strongly subjects are nudged
  toward group templates.
- Phase values wrap at ±π — be cautious when interpreting values near the boundary.

