# lncRNA_MS_Analysis
Differential analysis of long non-coding rnas from brain tissue in stagewise multiple sclerosis RNA-Seq Dataset

**Transcriptomic Analysis of lncRNAs in Multiple Sclerosis Brain Tissue**

This repository contains a complete reproducible RNA-Seq pipeline for identifying differentially expressed known and novel long non-coding RNAs (lncRNAs) across different stages of Multiple Sclerosis (MS) using human brain tissue samples.

 Objectives of my study:
- Identify differentially expressed lncRNAs across MS subtypes
- Detect stage-specific and shared lncRNA signatures
- Discover novel lncRNAs via transcript assembly
- Generate high-quality visualizations such as volcano plots for differntial analysis


| Step | Task | Script |
|------|------|--------|
| 1️⃣   | FASTQ Extraction & Integrity Check | `fastq_dump_then_md5.sh` |
| 2️⃣   | STAR Alignment & Transcript Assembly | `string_star.sh` |
| 3️⃣   | Read Quantification (FeatureCounts) | `generate_featurecounts.sh` |
| 4️⃣   | Transcript Merging & Novel lncRNA Discovery and extraction | `run_gffcompare_extract_lncRNA.sh` |
| 5️⃣   | Differential Expression & Visualization | `DESeq2_sumtechreps_withvisualisation.R` |


---

## 📌 Important Notes

Important Notes Here.
These are the considerations that were taken in my data edit according options and code that fits your data.
NOT ONE SIZE FITS ALL. IS IT SO?


- RNA-Seq data is **reverse-stranded**, confirmed using STAR `--quantMode` output
- STAR alignment was run in **two-pass mode**
- StringTie used for both **assembly and merging**
- FeatureCounts executed with **`-s 2`** to handle stranded data
- Both **known and novel lncRNAs** are quantified
- **Technical replicates** are collapsed after **batch effect checks**
- Metadata and replicate mapping are **generated programmatically**

---

## 📊 Visual Outputs

All results are stored in the `results/` folder and include:

- **PCA plots** – sample clustering by condition
- **Volcano plots** – for each MS subtype comparison
- **MA plots** – log2 fold change vs mean expression
- **Heatmaps** – of significant lncRNAs
- **Venn diagrams** – shared vs unique lncRNAs


---

## 🛠 Tools Used

| Tool | Purpose |
|------|---------|
| `fastq-dump`, `FastQC` | Download and quality control |
| `STAR` | Splice-aware alignment |
| `StringTie` | Transcript assembly |
| `featureCounts` | Read quantification |
| `gffcompare` | Novel transcript classification |
| `DESeq2` | Differential analysis |
| `R`, `ggplot2`, `pheatmap`, `EnhancedVolcano` | Plotting & statistical analysis |

---



<img width="1949" height="923" alt="image" src="https://github.com/user-attachments/assets/346b3c70-7589-4fd7-83d2-62ca1e1ab9d0" />

 🧾 Citation

To be added after peer-reviewed publication.  
Repository currently under pre-publication status.  



## 🙋‍♀️ Contact

For research collaborations or questions, please contact:  
📧 nandinipatel545@gmail.con
