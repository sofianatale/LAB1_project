# Structural HMM for Kunitz Domain Detection

## Overview

This project presents the development of a profile Hidden Markov Model (HMM) for detecting the Kunitz-type serine protease inhibitor domain. Unlike traditional sequence-based models, this classifier is trained on a manually curated structure-based alignment to better capture conserved 3D features typical of the Kunitz fold. The model’s robustness was assessed via a two-fold cross-validation strategy and benchmarked using standard classification metrics, including the Matthews Correlation Coefficient (MCC).

## Biological Background

The Kunitz domain (Pfam: `PF00014`) is a compact, cysteine-rich motif involved in the inhibition of serine proteases. Its characteristic fold, stabilized by three disulfide bridges, plays essential roles in blood coagulation, inflammation regulation, neuroprotection, and toxin function in various organisms. Due to high sequence variability, structural modeling provides an alternative route for more accurate domain detection.

## Tools and Software

This pipeline integrates a variety of command-line tools, Python packages, and web resources to build and evaluate a structural HMM for the Kunitz domain. Below is a complete list of software components and their respective purposes.

### System Tools

* **[HMMER 3.3.2](http://hmmer.org/):** Suite for building and querying profile Hidden Markov Models (HMMs).

  * `hmmbuild`: Trains the HMM from a multiple sequence alignment.
  * `hmmsearch`: Queries protein sequences against the HMM.
  * Options `--max` and `-Z 1000`: Used to improve sensitivity and normalize score calibration.

* **[PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/):** Used for structural alignment of PDB entries containing the Kunitz domain. The resulting `.ali` alignment file served as the input for HMM construction.

* **[CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/):** Reduces sequence redundancy by clustering protein sequences with ≥90% identity.

* **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download):** Filters out sequences highly similar to those used in training using `blastp` and `makeblastdb`.

* **[ChimeraX](https://www.cgl.ucsf.edu/chimerax/):** Visualizes and structurally superimposes predicted false positives and negatives against a reference structure (e.g., BPTI, PDB ID: 3TGI).

* **[WeBlogo / Skylign](https://weblogo.berkeley.edu/) & [Skylign](https://skylign.org/):** Generate sequence logos from the MSA used for HMM training.

### Python Dependencies

Used for data preprocessing, evaluation metrics, and visualizations:

| Library        | Purpose                                                                    |
| -------------- | -------------------------------------------------------------------------- |
| `pandas`       | Manage tabular data and read/write CSVs.                                   |
| `numpy`        | Perform numerical computations and array operations.                       |
| `matplotlib`   | Generate and save plots (e.g., ROC curves, score distributions).           |
| `seaborn`      | Visualize confusion matrices and classification metrics.                   |
| `scikit-learn` | Compute performance metrics like precision, recall, AUC, and MCC.          |
| `biopython`    | Used (via `Bio.SeqIO`) for parsing FASTA/FASTQ and manipulating sequences. |

To install all Python requirements:

```bash
pip install pandas numpy matplotlib seaborn scikit-learn biopython
```

### Web Tools and Databases

* **[Pfam – PF00014](https://www.ebi.ac.uk/interpro/entry/pfam/PF00014):** Reference for the Kunitz domain profile and seed alignments.
* **[InterPro](https://www.ebi.ac.uk/interpro/):** Used to confirm domain annotations and explore domain relationships.
* **[UniProt Downloads](https://www.uniprot.org/downloads):** Source of positive and negative sequences from Swiss-Prot in FASTA format.
* **[PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/):** Used again here for alignment-based domain extraction and comparison.
* **[Skylign / WebLogo](https://weblogo.berkeley.edu/):** Visualization of conserved positions in the HMM alignment.

### Installation via Conda

You can replicate the environment by running the following:

```bash
# Create and activate a dedicated conda environment
conda create -n hmm_kunitz python=3.10
conda activate hmm_kunitz

# Install core bioinformatics tools
conda install -c bioconda cd-hit
conda install -c bioconda hmmer
conda install -c bioconda blast-legacy

# Optional but recommended: Biopython
pip install biopython
```

## Pipeline Notebook

For a full walkthrough of the methodology with code, see the interactive notebook:
 [`pipeline.ipynb`](./pipeline.ipynb)

It contains data processing, model building, domain detection, and visualizations across validation folds.

## Repository Structure

### `/data/` — Input Data & Sequence Sets  
Contains raw and processed input datasets, including positive and negative sequence sets
- **clusters/**: clustering results and source PDB metadata  
- **filtered/**: filtered PDB FASTA and metadata files  
- **ids/**:  
  - `cv_splits/`: IDs for cross-validation folds (e.g., `pos_1.ids`, `neg_2.ids`)  
  - `filtering/`: IDs to retain or exclude (`to_keep.ids`, `to_remove.ids`)  
  - `random_sets/`: random positive/negative subsets  
  - `references/`: reference IDs (`all_kunitz.id`, `sp.id`)  
  - `structural/`: IDs of PDB sequences used in model training  
- **sequences/**:  
  - `cross_validation/`: FASTA files for cross-validation  
  - `negatives/`: negative sequence FASTA (`sp_negs.fasta`)  
  - `positives/`: human and non-human Kunitz positives  
  - `raw/`: unfiltered datasets from UniProt and Swiss-Prot  

### `/intermediate/` — Processing & Alignment Steps
Stores intermediate files generated during the pipeline, in particular BLAST results. This folder represents key transitional steps that bridge initial data preprocessing and final evaluation.
- Alignment files: `pdb_kunitz_rp.ali`, `pdb_kunitz_rp_formatted.ali`  
- Filtered datasets: `pdb_kunitz_rp.fasta`, `tmp_pdb_efold_ids.txt`  
- Cross-validation results: `pos_1.out`, `neg_2.out`, `fn_pos2.txt`  
- BLAST results: `pdb_kunitz_nr_23.blast -outfmt 7`  

### `/models/` — Profile HMM  
 Contains the structural HMM profile (`.hmm`) and related alignment files used for detection.
- `structural_model.hmm`: trained HMM built with `hmmbuild` from structure-based alignment  

### `/scripts/` — Automation & Evaluation  
Python scripts and notebooks used for sequence processing, evaluation, and plotting.
- `get_seq.py`: extracts sequences from ID lists  
- `performance.py`: computes MCC, precision, recall, AUC  
- `roc_curve.ipynb`: ROC curve visualization  
- `mcc_vs_threshold_plot.ipynb`: plots MCC vs E-value  
- `confusion_matrix.ipynb`: confusion matrix plotting  

### `/results/` — Classification Output  
Final output files such as performance tables, ROC curves, MCC plots, and confusion matrices.
- `classification/`
  - `fold1/`, `fold2/`: E-value prediction outputs for each fold (positive and negative sets).
  - `combined/`: Aggregated results across both folds.
- `evaluation/`
  - Contains performance summaries across thresholds  
    (e.g., `performance_set1_thresholds.txt`, `performance_set2_thresholds.txt`).
- `final_output/`
  - Final visual and tabular outputs, such as:  
    - `confusion_matrix_set_1.png`  
    - `results_set_2.txt`  
    - ROC and MCC plots, etc.

### `/figures/` — Visual material for report  
Contains all visual outputs generated during the analysis, also included in the final report:
- Confusion matrices and ROC curves  
- MCC vs E-value threshold plots  
- Sequence logos and structural superpositions

### `/final_report/` — Written report  
Includes the final PDF report with full methodology and discussion.
- `REPORT_LAB1_NATALE_SOFIA.pdf`: Final course project report

### `/` — Project Root  
- `pipeline.ipynb`: notebook describing and executing the main steps  
- `.gitignore`, `.gitattributes`: Git configuration files  
- `README.md`: project description (you are here)

## Conclusion

This project highlights the strength of **structure-guided profile HMMs** in domain detection, especially for **compact and divergent protein families** like Kunitz. The model generalizes well, achieves **perfect or near-perfect MCC**, and avoids false positives even under realistic class imbalance. The integration of structure-based features offers a reliable strategy for **high-throughput domain annotation** and could be extended to other fold families where sequence similarity fails.

## Supplementary Material

- Full project report: [`REPORT LAB1_NATALE_SOFIA.pdf`](./REPORT%20PROGETTO%20LAB1.pdf)
- GitHub repository: [https://github.com/sofianatale/LAB1_project](https://github.com/sofianatale/LAB1_project)

## References

See the full bibliography in the final report or access key tools used:

- [HMMER](http://hmmer.org/)
- [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/)
- [UniProt](https://www.uniprot.org/)
- [PDB](https://www.rcsb.org/)
- [Skylign](http://skylign.org/)
- [WebLogo](https://weblogo.berkeley.edu/)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)

## Project Details
**Author**: Sofia Natale  
**Degree Program**: MSc in Bioinformatics  
**University**: University of Bologna  
**Course**: Laboratory of Bioinformatics 1 - Module 2 
**Contact**: sofia.natale@studio.unibo.it 

>  _This work was carried out as part of the LAB1 course project and is not intended for production use without further validation._



