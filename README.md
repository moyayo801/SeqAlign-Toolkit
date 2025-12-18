# 🧬 SeqAlign-Toolkit: Advanced Biological Sequence Alignment

A Python implementation of classic dynamic programming algorithms for biological sequence analysis. This toolkit was developed to demonstrate expertise in algorithm design, evolutionary scoring matrices, and biological data processing.

## 🚀 Features

* **Needleman-Wunsch (Global Alignment):** Full sequence comparison using match/mismatch or substitution matrices.
* **Smith-Waterman (Local Alignment):** Identification of high-similarity regions with a zero-floor scoring rule.
* **Affine Gap Penalties (GOP/GEP):** Biologically accurate gap modeling using three-matrix ($M, I, D$) state transitions.
* **Top-K Local Search:** Ability to extract multiple local alignments surpassing a user-defined threshold $T$.
* **Substitution Matrices:** Built-in support for **BLOSUM62** (Proteins) and Nucleotide matrices (DNA/RNA).
* **Visual Analytics:** Heatmap generation for scoring matrices and gap states.

## 🧬 Mathematical Foundations

### Global Alignment (Needleman-Wunsch)
The engine computes the optimal global score using the recurrence:
$$
H_{i,j} = \max 
\begin{cases} 
H_{i-1,j-1} + S(a_i, b_j) & \text{(Match/Mismatch)} \\
H_{i-1,j} + d & \text{(Deletion)} \\
H_{i,j-1} + d & \text{(Insertion)}
\end{cases}
$$

### Affine Gap Model
To better model biological events, we separate gap opening ($o$) and extension ($e$):
$$I_{i,j} = \max(M_{i-1, j} - o, I_{i-1, j} - e)$$
$$D_{i,j} = \max(M_{i, j-1} - o, D_{i, j-1} - e)$$

---

## 🛠️ Usage

To run the specific implementations for the research tasks, use the provided scripts in the `usage/` directory:

1. **DNA Global Alignment (Task 1):**
```bash
python usage/dna_global_align.py

```

2. **Protein Global Alignment (Task 2):**
```bash
python usage/protein_constant_gap.py

```


3. **Protein Affine Gap Alignment (Task 3):**
```bash
python usage/protein_affine_gap.py

```


4. **Local Search & Thresholding (Task 4):**
```bash
python usage/local_search_top_k.py

```



## 📊 Visualizations

The toolkit generates heatmaps to analyze the scoring landscape. In local alignment mode, it highlights high-scoring sub-segments where Score \geq T.

## 📂 Project Structure

* `src/`: Core engine, scoring data, and visualization modules.
* `usage/`: Demonstration scripts for different alignment scenarios.
* `data/`: Sample sequence files in FASTA-like format.

```