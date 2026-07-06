<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/enteroscope/main/enteroscope.png" alt="EnteroScope Banner" width="100%">
</p>

<div align="center">
  
# 🔬 EnteroScope

### **A Tailored Computational Workflow Enabling Rapid, User-Friendly Genotyping and Epidemiological Surveillance of the Enterobacter cloacae Complex**

#### **Complete CREC (Carbapenem‑Resistant *Enterobacter cloacae*) genomic analysis in minutes — not hours**

![Version](https://anaconda.org/bbeckley-hub/enteroscope/badges/version.svg)
![Latest Release Date](https://anaconda.org/bbeckley-hub/enteroscope/badges/latest_release_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/enteroscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/enteroscope/badges/license.svg)

[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/enteroscope)](https://hub.docker.com/r/bbeckleyhub/enteroscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/enteroscope/latest)](https://hub.docker.com/r/bbeckleyhub/enteroscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/enteroscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/enteroscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)
![Conda Downloads](https://img.shields.io/conda/dn/bioconda/enteroscope?label=Conda%20Downloads)

[![Powered by 🧠](https://img.shields.io/badge/powered%20by-science%20🔬-purple)](https://github.com/bbeckley-hub/enteroscope)
[![Coffee](https://img.shields.io/badge/built%20with-%E2%98%95%20coffee-orange)](https://github.com/bbeckley-hub/enteroscope)
[![Made with ❤️](https://img.shields.io/badge/made%20with-%E2%9D%A4%EF%B8%8F-red)](https://github.com/bbeckley-hub/enteroscope)
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)
[![Made for Research](https://img.shields.io/badge/made%20for-Research-0066cc.svg)](https://github.com/bbeckley-hub/enteroscope)

[![Documentation](https://img.shields.io/badge/docs-mkdocs-526CFE?logo=materialformkdocs)](https://bbeckley-hub.github.io/enteroscope)
[![RST Badge](https://img.shields.io/badge/documentation-RST-4CAF50.svg)](https://www.sphinx-doc.org/)
[![Last Commit](https://img.shields.io/github/last-commit/bbeckley-hub/enteroscope)](https://github.com/bbeckley-hub/enteroscope/commits)
[![Contributors](https://img.shields.io/github/contributors/bbeckley-hub/enteroscope)](https://github.com/bbeckley-hub/enteroscope/graphs/contributors)
[![Security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![CI](https://img.shields.io/github/actions/workflow/status/bbeckley-hub/enteroscope/ci.yml?branch=main&label=CI)](https://github.com/bbeckley-hub/enteroscope/actions)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](https://github.com/bbeckley-hub/enteroscope/tests)

[![Speed](https://img.shields.io/badge/Speed-12%20min%2F24%20samples-FF6D00.svg)](https://github.com/bbeckley-hub/enteroscope#performance-benchmarks)
[![CREC](https://img.shields.io/badge/Target-Carbapenem--Resistant%20E.%20cloacae-00BCD4.svg)](https://github.com/bbeckley-hub/enteroscope)
[![Lineages](https://img.shields.io/badge/Lineages-12%20complex%20species-673AB7.svg)](https://github.com/bbeckley-hub/enteroscope)
[![AI Ready](https://img.shields.io/badge/AI-Ready%20Reports-00C853.svg)](https://github.com/bbeckley-hub/enteroscope#ai-integration-guide)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/enteroscope)](https://github.com/bbeckley-hub/enteroscope/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/enteroscope)](https://github.com/bbeckley-hub/enteroscope/stargazers)
[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/enteroscope/#summary)

</div>

---

## 📋 Table of Contents

- [What's New in v1.1.0](#-whats-new-in-v110)
- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [⚡ Quick Start (CLI)](#-quick-start-cli)
- [🔧 Installation (CLI)](#-installation-cli)
- [🐳 Docker Usage](#-docker-usage)
- [🖥️ Singularity for HPC Clusters](#️-singularity-for-hpc-clusters)
- [🔗 Integrated External Tools & Dependencies](#-integrated-external-tools--dependencies)
- [🚀 Usage Guide (CLI)](#-usage-guide-cli)
- [📁 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [🔬 Validation & Accuracy](#-validation--accuracy)
- [🤖 AI Integration Guide](#-ai-integration-guide)
- [🔮 Future Development](#-future-development)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)
- [🐛 Troubleshooting](#-troubleshooting)
- [📚 Citation](#-citation)
- [🙏 Acknowledgements](#-acknowledgements)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)
- [📚 Third‑Party Tool Citations](#-third-party-tool-citations)

---

## 🆕 What's New in v1.1.0

**Release Date: 2026-07-06**

We’ve upgraded EnteroScope with a fully **temporary‑directory execution model** (HPC‑friendly), added **point mutation reporting**, introduced **dynamic filtering thresholds**, and improved the Ultimate Reporter with a dedicated **Mutations tab** and **grouping by MLST**.

### Major Enhancements

- ✅ **HPC‑Ready Orchestrator**: All modules now run inside temporary directories (`/tmp` by default) and are automatically cleaned up. This eliminates permission errors on clusters and keeps your working directory pristine.
- ✅ **Point Mutation Reporting (AMRfinderPlus)**: AMRfinderPlus now runs with `--plus` by default, and the Ultimate Reporter displays a gene‑centric **Mutation tab** with grouping by MLST – instantly see which clones carry mutations in `gyrA`, `parC`, `rpoB`, or 23S rRNA.
- ✅ **Dynamic Threshold Flags**:
  - `--amr-min-identity` and `--amr-min-coverage` – custom cutoffs for AMR hits (0..1).
  - `--skip-amr-mutations` – disable point mutation reporting if you only want acquired genes.
  - `--amr-force-update` – force a fresh AMR database download before analysis.
  - `--abricate-min-identity` and `--abricate-min-coverage` – custom cutoffs for ABRicate (default 80%).
- ✅ **Better Logging**: Each module’s stdout/stderr is captured in `enteroscope_run.log` for full traceability.
- ✅ **Cleaner Output**: Mutation summary files (`mutation_summary.html`, `mutation_summary.tsv`, `mutation_master_summary.json`) are now placed inside `enteroscope_amrfinder_results/` and automatically copied to the summary module.
- ✅ **Ultimate Reporter Updates**:
  - **Mutations tab** is always shown; if no mutations are found, a friendly message appears.
  - **Grouping by MLST** is now available on the AMR, Virulence, BACMET, Plasmids, and Mutations tabs – click a button to reorganise genome lists by sequence type.
  - **Sortable tables** on every tab.
  - **Scientific quotes** rotate after each module and at the end for inspiration.
- ✅ **Docker & Singularity**: The container now runs without permission errors – use `--writable-tmpfs` for Singularity and all outputs are owned by your user.

---

## 🎯 Overview

**EnteroScope** is an automated, locally‑executable pipeline designed specifically for comprehensive genomic surveillance of the **Enterobacter cloacae** – an opportunistic pathogen that increasingly important in healthcare‑associated infections and carbapenem resistance (CREC). It integrates **six essential genotyping methods** into a single, cohesive workflow.

### 🌍 The Problem
- **Emerging Threat**: *Enterobacter cloacae* is a leading cause of neonatal sepsis and hospital outbreaks, with CREC mortality rates exceeding 40%.
- **Fragmented Workflows**: Traditional analysis requires 5+ separate tools (MLST, AMR, virulence, plasmids) with conflicting dependencies.
- **Data Overload**: Raw outputs lack epidemiological context and actionable insights.
- **Resource Barriers**: Web services raise privacy concerns; local solutions are hard to install.

### 💡 Our Solution
EnteroScope delivers:
- ✅ **Single‑command installation** via Conda or Docker.
- ✅ **3‑4 hrs complete analysis** (100 samples, 16 cores).
- ✅ **100% local execution** – your data stays private.
- ✅ **Interactive HTML reports** with rotating educational facts about *E. cloacae*.
- ✅ **Gene‑centric tables** (scrollable, wrapped genome lists).
- ✅ **Automatic detection of critical resistance**: KPC, NDM, OXA‑48, VIM, IMP, mcr, tet(X).
- ✅ **Point mutation reporting** with grouping by MLST.

**Perfect for**: Clinical labs, outbreak investigations (CREC), public health surveillance, and research studies.

---

## 📊 Sample Output

See a complete interactive report generated by EnteroScope:

[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/enteroscope/#summary)

*The report includes AMR, Virulence, BACMET, Plasmid, and **Mutations** tabs, each with grouping by MLST, filter buttons, sortable headers, and scrollable genome lists.*

---

## ✨ Key Features

| Module | Purpose | Key Outputs | Speed |
|--------|---------|-------------|-------|
| **FASTA QC** | Assembly quality (N50, GC%, contigs) | HTML, TSV, JSON | <30 sec |
| **MLST** | Sequence typing (7 housekeeping genes) | ST, allele profile | <1 min |
| **AMR Profiling** | Resistance genes (AMRFinderPlus) | Carbapenemases, ESBLs, mcr, tet(X) | 2‑3 min |
| **ABRicate** | Multi‑database (9 DBs) – virulence/plasmids | VFDB, ResFinder, CARD, PlasmidFinder | 3‑4 min |
| **BACMET2** | Biocide & heavy metal resistance | qac, mer, ars, czc, etc. | 1‑2 min |
| **PlasmidFinder** | Replicon typing (IncL/M, IncFII, IncHI2) | Outbreak plasmid tracking | 1‑2 min |
| **Mutations (AMRfinderPlus)** | Point mutations (`gyrA`, `parC`, `rpoB`, 23S) | Gene‑centric table, grouped by ST | 1‑2 min |
| **Ultimate Reporter** | Interactive HTML (gene‑centric, scrollable lists, grouping) | All results in one place | <1 min |

### 🛡️ Enterobacter‑Specific Innovations
- **Intrinsic AmpC handling**: Automatic flagging of ACT/MIR genes (chromosomal) with clinical interpretation.
- **Carbapenemase detection**: Word‑boundary pattern matching for KPC, NDM, OXA‑48‑like, VIM, IMP, **IMI** (chromosomal).
- **Colistin & tigecycline resistance**: mcr‑1 to mcr‑10, tet(X) variants, efflux pumps (AcrAB‑TolC, OqxAB).
- **Plasmid epidemiology**: Highlights IncL/M (OXA‑48), IncFII (KPC), IncHI2 (NDM/mcr).
- **Educational facts**: 20 rotating facts about *E. cloacae* in the report header.
- **Point mutation surveillance**: Tracks resistance‑associated mutations in gyrA, parC, rpoB, 23S rRNA.

### 🚀 Performance Advantages

- **Low memory footprint**: runs on 4 GB RAM, scales to HPC clusters.
- **Dynamic resource allocation** (psutil).
- **Temporary‑directory execution** – no leftover files, HPC‑friendly.

---

## ⚡ Quick Start (CLI)

### Install in 60 seconds
```bash
# Conda (Recommended) Recipe submitted to Bioconda
conda create -n enteroscope -c conda-forge -c bioconda -c enteroscope -y
conda activate enteroscope

# Mamba (Faster)
mamba create -n enteroscope -c conda-forge -c bioconda -c enteroscope -y
mamba activate enteroscope
```

### Run your first analysis
```bash
# Single genome
enteroscope -i genome.fasta -o results/

# Batch processing (100 genomes)
enteroscope -i "*.fna" -o batch_results --threads 16
# Complete in ~4 hrs! 🎉
```

### Advanced examples
```bash
# AMR with custom identity & coverage, skip mutations
enteroscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations

# ABRicate with custom thresholds
enteroscope -i "*.fna" -o results --abricate-min-identity 85 --abricate-min-coverage 85

# Force AMR database update before analysis
enteroscope -i "*.fna" -o results --amr-force-update

# Keep temporary directories for debugging
enteroscope -i "*.fna" -o results --keep-temp
```

---

## 🔧 Installation (CLI)

### System Requirements
| Resource | Minimum | Recommended | Production |
|----------|---------|-------------|------------|
| CPU Cores | 2 | 8+ | 16+ |
| RAM | 4 GB | 8 GB | 16 GB |
| Storage | 2 GB | 10 GB | 50 GB+ |
| OS | Linux, macOS, WSL2 | Linux | Linux Cluster |

### Step‑by‑Step

```bash
# 1. Install Miniconda (if needed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# 2. Add channels
conda config --add channels conda-forge
conda config --add channels bioconda

# 3. Create environment
conda create -n enteroscope python=3.9 enteroscope -y
conda activate enteroscope

# 4. Update databases (recommended)
abricate --setupdb
enteroscope --update-amr-db   # Only needed once
```

---

## 🐳 Docker Usage

### Quick Start
```bash
# Pull the image
docker pull bbeckleyhub/enteroscope:latest

# Test
docker run --rm bbeckleyhub/enteroscope:latest --help

# Analyse your data
docker run --rm \
  -v $(pwd)/genomes:/data/input \
  -v $(pwd)/results:/data/output \
  bbeckleyhub/enteroscope:latest \
  -i "/data/input/*.fasta" -o /data/output -t 4

# Fix ownership (output files owned by root)
sudo chown -R $USER:$USER ./results
```

---

## 🖥️ Singularity for HPC Clusters

On HPC clusters that support [Singularity/Apptainer](https://sylabs.io/singularity/), you can run EnteroScope **without `sudo`** and output files will belong to your user automatically.

> **Important:** EnteroScope writes temporary files inside its own installation directory. Singularity mounts containers as read‑only by default, so you **must** add the `--writable-tmpfs` flag.

### Option A: Direct pull (if network allows)
```bash
singularity pull enteroscope.sif docker://bbeckleyhub/enteroscope:latest
singularity run --writable-tmpfs -B $(pwd):/data enteroscope.sif -i "/data/*.fasta" -o /data/output
```

### Option B: Convert from local Docker image (when pull fails on HPC)

**On a machine with Docker:**
```bash
docker pull bbeckleyhub/enteroscope:latest
docker save bbeckleyhub/enteroscope:latest -o enteroscope.tar
singularity build enteroscope.sif docker-archive://enteroscope.tar
```
Copy `enteroscope.sif` to the HPC. Then:

```bash
singularity run --writable-tmpfs -B $(pwd):/data enteroscope.sif -i "/data/*.fasta" -o /data/output
```

### Explanation of flags
| Flag | Purpose |
|------|---------|
| `--writable-tmpfs` | Creates a temporary writable overlay – **required** for EnteroScope. |
| `-B $(pwd):/data` | Binds current directory to `/data` inside the container. |
| `-i "/data/*.fasta"` | Input pattern (quotes prevent host expansion). |
| `-o /data/output` | Output directory (appears as `./output` on host). |

After the run, all results in `./output` are owned by your HPC user – **no `sudo chown` needed**.

---

## 🔗 Integrated External Tools & Dependencies

These tools are **automatically installed as Conda dependencies** (see `environment.yml`). EnteroScope’s MIT license does not cover them; each tool is used under its own license.

| Tool | Purpose | Source | License |
|------|---------|--------|---------|
| **MLST** | Sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | AMR gene & mutation detection | [ncbi/amr](https://github.com/ncbi/amr) | Public domain |
| **PlasmidFinder** | Replicon typing | [genomicepidemiology/plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder) | Apache‑2.0 |
| **BACMET2** | Biocide/metal resist. | [bacmet2](https://bacmet.biomedicine.gu.se/) | Academic free |

---

## 🚀 Usage Guide (CLI)

### Basic commands
```bash
# Single genome
enteroscope -i genome.fasta -o results/

# Batch processing with wildcards
enteroscope -i "*.fna" -o results_2025 --threads 8

# Skip modules
enteroscope -i sample.fna -o results --skip-abricate
```

### Full options
```
usage: enteroscope [-h] [-i INPUT] [-o OUTPUT] [-t THREADS] [--skip-qc] [--skip-mlst] [--skip-abricate] [--skip-amr] [--skip-summary]
                   [--amr-min-identity AMR_MIN_IDENTITY] [--amr-min-coverage AMR_MIN_COVERAGE] [--skip-amr-mutations]
                   [--amr-force-update] [--abricate-min-identity ABRICATE_MIN_IDENTITY]
                   [--abricate-min-coverage ABRICATE_MIN_COVERAGE] [--update-amr-db] [--force-update-amr-db] [--keep-temp]
                   [--clean-output] [--version]
```

---

## 📁 Output Structure

```
Enteroscope_results/
├── fasta_qc_results/                # QC metrics (HTML, TSV, JSON)
├── mlst_results/                    # MLST typing (ST, allele profile)
├── enteroscope_abricate_results/    # ABRicate outputs (9 databases)
│   ├── enteroscope_card_summary_report.html
│   ├── enteroscope_vfdb_summary_report.html
│   ├── enteroscope_plasmidfinder_summary_report.html
│   └── ...
├── enteroscope_amrfinder_results/   # AMRfinderPlus results
│   ├── enteroscope_amrfinder_summary_report.html
│   ├── mutation_summary.html        # Gene-centric mutation table (for internal use)
│   ├── mutation_summary.tsv
│   └── mutation_master_summary.json
├── ENTERO_ULTIMATE_REPORTS/         # Final HTML dashboard
│   ├── enteroscope_ultimate_report.html  # Main interactive report
│   ├── enteroscope_ultimate_report.json  # Full data for AI
│   ├── amr_genes.csv
│   ├── virulence_genes.csv
│   ├── bacmet_genes.csv
│   ├── plasmid_markers.csv
│   ├── mutations.csv                     # Exported point mutation table
│   ├── enterobacter_samples.csv
│   └── mlst_distribution.csv
└── enteroscope_run.log              # Detailed log file
```

---

## 🔍 Analytical Modules

### 1. FASTA QC
- Metrics: N50, N70, N90, L50, GC%, total length, contigs.
- Formats: HTML (interactive), TSV, JSON.

### 2. MLST (Enterobacter cloacae complex)
- Loci: *dnaA, fusA, gyrB, leuS, pyrG, rplB, rpoB*
- Database: PubMLST (updated via `mlst --make`)

### 3. AMR Profiling (AMRFinderPlus)
- >5,000 resistance genes; critical flags: carbapenemases (KPC, NDM, OXA‑48‑like, VIM, IMP), ESBLs, mcr, tet(X)
- Risk categories: Critical / High / Standard
- **Point mutations**: gyrA (quinolone), parC (quinolone), rpoB (rifampin), 23S rRNA (linezolid)

### 4. ABRicate Screening
- Databases: VFDB, ResFinder, CARD, PlasmidFinder, MegaRes, NCBI, ARG‑ANNOT, ECOH, EcoLi_VF
- Thresholds: user‑adjustable via `--abricate-min-identity` and `--abricate-min-coverage` (default: 80%)

### 5. BACMET2
- Biocides: qac, cepA, formA
- Heavy metals: mer, ars, cop, sil, chr, cad, znt, czc

### 6. PlasmidFinder
- Replicons: IncL/M (OXA‑48), IncFII (KPC), IncHI2 (NDM/mcr), Col‑like, etc.
- Category: “Incompatibility group X” or “Colicin plasmid”

### 7. Mutation Detection (AMRfinderPlus)
- AMRfinderPlus reports point mutations in genes like **gyrA**, **parC**, **rpoB**, and **23S rRNA**.
- EnteroScope creates a gene‑centric mutation table with grouping by MLST, shown in the Ultimate Reporter’s **Mutations** tab.
- Mutations can be exported as CSV for further analysis.

### 8. Ultimate Reporter
- Interactive HTML with scrollable genome lists, filter buttons, combination tables, rotating educational facts, AI‑ready JSON, and **grouping by MLST** on all gene‑centric tabs (AMR, Virulence, BACMET, Plasmids, Mutations).

---

## 🔬 Validation & Accuracy

### Reference strain validation (100% concordance)
| Strain | Expected | EnteroScope result |
|--------|----------|--------------------|
| E. cloacae ATCC 13047 | ST114 | ✅ ST114 |
| E. hormaechei subsp. oharae | ST78 | ✅ ST78 |
| KPC‑producing isolate | ST171, blaKPC‑2 | ✅ ST171, blaKPC‑2 |

### Clinical isolate pilot (n=24)
- Carbapenemase genes: 8 isolates (33%): KPC‑2 (5), NDM‑1 (2), OXA‑48 (1)
- Colistin resistance (mcr‑9): 2 isolates
- Dominant STs: ST171 (6), ST78 (5), ST114 (4)

---

## 🤖 AI Integration Guide

EnteroScope’s HTML reports are designed for AI analysis.

### Example questions:
- “Which samples carry carbapenemase genes (KPC, NDM, OXA‑48)?”
- “List all isolates with mcr genes (colistin resistance).”
- “What are the common plasmid replicons in ST171?”
- “Show me the co‑occurrence of blaKPC and IncFII plasmids.”
- “Summarise the AMR profile of sample XYZ.”
- “Which isolates have gyrA mutations?”
- “Show me the mutation table grouped by MLST.”

Upload the HTML or JSON file to ChatGPT, Claude, or Gemini and ask away.

---

## 🔮 Future Development

- Raw read support (FASTQ → assembly → typing)
- Real‑time database updates
- Machine learning module for AMR phenotype prediction
- Integration with local CREC outbreak databases

---

## ❓ Frequently Asked Questions

**Q: Does EnteroScope detect the intrinsic IMI carbapenemase?**  
A: Yes. The word‑boundary pattern matching in `categorize_gene()` specifically handles `blaIMI` without false positives (e.g., not matching “fimI”).

**Q: Can I use it for *E. coli* or *Klebsiella*?**  
A: No. Please use our dedicated tools: **EcoliTyper** and **Kleboscope**.

**Q: The AMR databases are huge. Why aren’t they in the repo?**  
A: The databases are downloaded on‑the‑fly via the `download_databases.sh` script. This keeps the repository small and the installation legal.

**Q: How do I disable point mutation reporting?**  
A: Use the `--skip-amr-mutations` flag.

**Q: Where are mutation summary files?**  
A: They are inside `enteroscope_amrfinder_results/` and are automatically copied to the summary module. The Ultimate Reporter’s **Mutations** tab displays them.

---

## 🐛 Troubleshooting

| Problem | Solution |
|---------|----------|
| `amrfinder` not found | Run `enteroscope --update-amr-db` |
| `abricate` databases missing | Run `abricate --setupdb` |
| Permission denied (Docker) | Use `sudo chown -R $USER:$USER ./output` |
| Singularity “No such file or directory” | Add `--writable-tmpfs` flag |
| MLST fails (unknown scheme) | Run `mlst --make` to update PubMLST data |
| Temporary directories not cleaning up | Use `--keep-temp` to debug; normally they are auto‑removed |

If problems persist, open an [issue on GitHub](https://github.com/bbeckley-hub/enteroscope/issues).

---

## 📚 Citation

If you use EnteroScope in your research, please cite:

> Beckley, B.(2026). EnteroScope: a species‑optimized computational pipeline for rapid and accessible *Enterobacter cloacae* genotyping and surveillance.
```bibtex
@article{beckley2026enteroscope,
  title={EnteroScope: A Tailored Computational Workflow Enabling Rapid, User-Friendly Genotyping and Epidemiological Surveillance of the Enterobacter cloacae Complex},
  author={Beckley Brown},
  journal={In preparation},
  volume={?},
  pages={??},
  year={??},
  doi={??} 
}
```

---

## 🙏 Acknowledgements

We thank the open‑source bioinformatics community, especially:
- **Torsten Seemann** (MLST, ABRicate)
- **NCBI AMRFinderPlus team**
- **CGE** (PlasmidFinder)
- **BACMET2** creators
- **University of Ghana Medical School** for support

---

## 👥 Authors & Contact

**Brown Beckley** (Primary developer)  
📧 brownbeckley94@gmail.com  
🐙 [bbeckley-hub](https://github.com/bbeckley-hub)  
🔗 [LinkedIn](https://www.linkedin.com/in/brown-beckley-190315319)

---

## 📄 License

- **EnteroScope pipeline code**: MIT License (see [LICENSE](LICENSE)).
- **Third‑party tools**: each under its own license (see section above).

---

## 📚 Third‑Party Tool Citations

<details>
<summary>Click to expand BibTeX citations</summary>

```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  url = {https://github.com/tseemann/mlst}
}

@article{jolley_pubmlst_2018,
  author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
  title = {Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications},
  journal = {Wellcome Open Research},
  volume = {3},
  pages = {124},
  year = {2018},
  doi = {10.12688/wellcomeopenres.14826.1}
}

@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  url = {https://github.com/tseemann/abricate}
}

@article{feldgarden_amrfinderplus_2019,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2019},
  doi = {10.1038/s41598-021-91456-0}
}

@article{carattoli_plasmidfinder_2014,
  author = {Carattoli, A. et al.},
  title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {7},
  pages = {3895-3903},
  year = {2014},
  doi = {10.1128/AAC.02412-14}
}
```
</details>

---

<div align="center">

## 🚀 Ready to fight CREC?

[![Get Started CLI](https://img.shields.io/badge/GET_STARTED_CLI-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/enteroscope#-quick-start-cli)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/enteroscope/issues)

**From fragmented data to actionable insights.**

*EnteroScope: Precision surveillance for the Enterobacter cloacae.*

⭐ **If this tool helps your research, please star the repository!** ⭐

</div>
