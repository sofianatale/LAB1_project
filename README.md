# Structural HMM for Kunitz Domain Detection

## Overview

This project presents the development of a profile Hidden Markov Model (HMM) for detecting the Kunitz-type serine protease inhibitor domain. Unlike traditional sequence-based models, this classifier is trained on a manually curated structure-based alignment to better capture conserved 3D features typical of the Kunitz fold. The model’s robustness was assessed via a two-fold cross-validation strategy and benchmarked using standard classification metrics, including the Matthews Correlation Coefficient (MCC).

## Biological Background

The Kunitz domain (Pfam: `PF00014`) is a compact, cysteine-rich motif involved in the inhibition of serine proteases. Its characteristic fold, stabilized by three disulfide bridges, plays essential roles in blood coagulation, inflammation regulation, neuroprotection, and toxin function in various organisms. Due to high sequence variability, structural modeling provides an alternative route for more accurate domain detection.

## Repository Structure

### `/data/` — Input Data & Sequence Sets  
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
- Alignment files: `pdb_kunitz_rp.ali`, `pdb_kunitz_rp_formatted.ali`  
- Filtered datasets: `pdb_kunitz_rp.fasta`, `tmp_pdb_efold_ids.txt`  
- Cross-validation results: `pos_1.out`, `neg_2.out`, `fn_pos2.txt`  
- BLAST results: `pdb_kunitz_nr_23.blast -outfmt 7`  

### `/models/` — Profile HMM  
- `structural_model.hmm`: trained HMM built with `hmmbuild` from structure-based alignment  

### `/scripts/` — Automation & Evaluation  
- `get_seq.py`: extracts sequences from ID lists  
- `performance.py`: computes MCC, precision, recall, AUC  
- `roc_curve.ipynb`: ROC curve visualization  
- `mcc_vs_threshold_plot.ipynb`: plots MCC vs E-value  
- `confusion_matrix.ipynb`: confusion matrix plotting  

### `/results/` — Classification Output  
- `classification/`:  
  - `fold1/`, `fold2/`: E-value predictions for positive and negative sets  
  - `combined/`: all results across both folds  
- `evaluation/`:  
  - Performance across thresholds (e.g., `performance_set1_thresholds.txt`)  
- `final output/`:  
  - `confusion_matrix_set_1.png`, `results_set_2.txt`, etc.

### `/` — Project Root  
- `pipeline.ipynb`: notebook describing and executing the main steps  
- `.gitignore`, `.gitattributes`: Git configuration files  
- `README.md`: project description (you are here)

## Tools and Software

- `HMMER` suite (v3.4) – for building and querying profile HMMs  
- `PDBeFold` – for 3D structural alignment of Kunitz-domain proteins  
- `CD-HIT` – to reduce redundancy in PDB datasets  
- `BLAST+` – to filter out similar sequences between sets  
- `UCSF ChimeraX` – for structural superimposition  
- `Python` – for custom preprocessing, subset generation, and performance evaluation  

## Methodology
For a detailed step-by-step implementation of the entire pipeline (from data collection to evaluation), see the Jupyter notebook [`pipeline.ipynb`](./pipeline.ipynb).

### 1. Dataset Preparation
- **Positive dataset**: Proteins annotated with the InterPro ID `IPR002223` (Kunitz-type domain) were downloaded from UniProt.
- Sequences were split by species:

  ```bash
  grep "Homo sapiens" all_kunitz.fasta > human_kunitz.fasta
  grep -v "Homo sapiens" all_kunitz.fasta > nothuman_kunitz.fasta
   ```
  
- PDB selection: Kunitz-domain sequences from high-resolution structures (≤3.5 Å) were retrieved using a custom RCSB PDB query and filtered:

 ```bash
awk -F ',' '{if (length($2)>0) {name=$2}; print name, $3, $4, $5}' rcsb_pdb_custom_report.csv | grep PF00014 > pdb_kunitz_customreported.fasta
 ```
### 2. Redundancy Filtering and Representative Selection
CD-HIT was used to cluster PDB sequences and reduce redundancy:

 ```bash
cd-hit -i pdb_kunitz_customreported.fasta -o pdb_kunitz_customreported_nr.fasta -c 0.9
 ```
Representative sequences were extracted based on clustering:

 ```bash
awk '$5 == 1 {print $1}' pdb_kunitz.clusters.txt > pdb_kunitz_rp.ids
 ```
Final representative FASTA file:

 ```bash
for i in $(cat pdb_kunitz_rp.ids); do
  grep -A 1 "^>$i" pdb_kunitz_customreported.fasta | head -n 2 >> pdb_kunitz_rp.fasta
done
 ```
### 3. Structural Alignment and HMM Building
- Manual structural alignment was performed using PDBeFold.
- The alignment was reformatted for HMMER:

```bash
awk '{ if (substr($1,1,1)==">") print "\\n" toupper($1); else printf "%s", toupper($1) }' pdb_kunitz_rp.ali > pdb_kunitz_rp_formatted.ali
```
The HMM profile was built using:

```bash
hmmbuild structural_model.hmm pdb_kunitz_rp_formatted.ali
```
### 4. Filtering for Model Evaluation
To avoid redundancy between training and testing sets, a BLASTP comparison was performed:

```bash
blastp -query pdb_kunitz_rp.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > to_remove.ids
```
The filtered positive set was created by removing high-similarity hits:

```bash
python scripts/get_seq.py --ids to_keep.ids --fasta all_kunitz.fasta --out ok_kunitz.fasta
```
### 5. Negative Dataset Construction
- The full Swiss-Prot dataset was filtered to remove known Kunitz-domain proteins.
- Balanced negative subsets were sampled from this filtered set:
  
```bash
python scripts/get_seq.py --ids neg_1.ids --fasta sp_negs.fasta --out neg_1.fasta
python scripts/get_seq.py --ids neg_2.ids --fasta sp_negs.fasta --out neg_2.fasta
```
### 6. Cross-Validation Setup
- The filtered positive set (ok_kunitz.fasta) was split into two folds (183 sequences each).
- Corresponding negative sets were matched by size for class balance.

```bash
python scripts/get_seq.py --ids pos_1.ids --fasta ok_kunitz.fasta --out pos_1.fasta
python scripts/get_seq.py --ids pos_2.ids --fasta ok_kunitz.fasta --out pos_2.fasta
```
### 7. Domain Search with HMM
hmmsearch was used to scan each fold:

```bash
hmmsearch --max -Z 1000 --tblout pos_1.out structural_model.hmm pos_1.fasta
hmmsearch --max -Z 1000 --tblout neg_1.out structural_model.hmm neg_1.fasta
```
Outputs were formatted into .class files with E-values. Missing hits were assigned E = 10.0.

### 8. Performance Evaluation
The script performance.py calculated MCC, precision, recall across thresholds from 1e-1 to 1e-10:

```bash
python scripts/performance.py results/classification/combined/set_1.class
python scripts/performance.py results/classification/combined/set_2.class
```
Plots were generated with the following Jupyter notebooks:

```bash
confusion_matrix.ipynb                  # Confusion matrices for both validation folds
roc_curve.ipynb                         # ROC curves and AUC computation
mcc_vs_threshold_plot.ipynb            # MCC values plotted against e-value thresholds
```
At the optimal threshold (1e-3), the model achieved perfect or near-perfect MCC, with no false positives, and an AUC of 1.000.

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



