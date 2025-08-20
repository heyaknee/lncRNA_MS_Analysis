# lncRNA_MS_Analysis
Differential analysis of long non-coding rnas from brain tissue in stagewise multiple sclerosis RNA-Seq Dataset

**Transcriptomic Analysis of lncRNAs in Multiple Sclerosis Brain Tissue**

This repository contains a complete reproducible RNA-Seq pipeline for identifying differentially expressed known and novel long non-coding RNAs (lncRNAs) across different stages of Multiple Sclerosis (MS) using human brain tissue samples.

 Objectives of my study:
- Identify differentially expressed lncRNAs across MS subtypes
- Detect stage-specific and shared lncRNA signatures
- Discover novel lncRNAs via transcript assembly
- Generate high-quality visualizations such as volcano plots for differntial analysis
- WGCNA 

| Step | Task                                  | Script                                   |
|------|---------------------------------------|-------------------------------------------|
| 1️⃣   | FASTQ Extraction & Integrity Check    | `scripts/01_fastq_dump_and_md5check.sh`   |
| 2️⃣   | STAR Alignment & Transcript Assembly | `scripts/02_star_alignment_stringtie.sh`  |
| 3️⃣   | Read Quantification (FeatureCounts)  | `scripts/03–05_featurecounts_*.sh`        |
| 4️⃣   | Transcript Merging & Novel lncRNA    | `scripts/06_gffcompare_novel_extraction.sh`|
| 5️⃣   | lncRNA Extraction from GTF           | `scripts/extract_lncrnas_from_gtf.sh`     |
| 6️⃣   | Differential Expression & Visualization | `notebooks/DESeq2_analysis.Rmd`, `notebooks/lncRNA_visualizations.ipynb` |
| 7️⃣   | WGCNA Network Analysis of lncRNAs    | `scripts/07_wgcna_setup.R`, `notebooks/WGCNA_analysis.Rmd` |





Important Notes Here.
These are the considerations that were taken in my data edit according options and code that fits your data.
NOT ONE SIZE FITS ALL. IS IT SO?

- RNA-Seq data used here is **reverse-stranded**. 
- STAR alignment is performed using **two-pass mode**
- StringTie used for **transcript assembly** and merging
- FeatureCounts run with appropriate strandedness (`-s 2`) (Tip: How to check strandness? 
- Known and novel lncRNAs are both analyzed
- Replicates are summed after checking batch effects
- Metadata and replicate mapping are generated programmatically


## 📊 Visual Outputs

Plots are located in the `results/` folder and include:

- PCA for condition separation
- Volcano & MA plots for each MS subtype
- Heatmaps of significant lncRNAs
- Venn diagrams for shared and unique lncRNAs


## 🧬 Novel lncRNA Discovery

- Merged transcript assemblies by biological replicate
- Compared against reference annotations using `gffcompare`
- Novel transcripts with class code `"u"` retained as potential novel lncRNAs



🛠 Tools & Libraries

| Tool/Library | Purpose |
|--------------|---------|
| **fastq-dump** | Raw data extraction |
| **FastQC** | Quality control |
| **STAR** | Splice-aware alignment |
| **StringTie** | Transcript assembly |
| **featureCounts** | Read quantification |
| **gffcompare** | Novel transcript comparison |
| **DESeq2**, **edgeR** | Differential analysis |
| **Python (pandas, seaborn)** | Data handling & visualization |

<img width="1949" height="923" alt="image" src="https://github.com/user-attachments/assets/346b3c70-7589-4fd7-83d2-62ca1e1ab9d0" />

 🧾 Citation

To be added after peer-reviewed publication.  
Repository currently under pre-publication status.  



## 🙋‍♀️ Contact

For research collaborations or questions, please contact:  
📧 nandinipatel545@gmail.con
