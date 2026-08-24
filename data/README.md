# data/

Raw data is not committed to this repository (large binary files: gene expression matrices, reference genomes, curated GEMs). This describes how to regenerate it.

## Source

Flores Ramos, S. et al. "Metatranscriptomics uncovers diurnal functional shifts in bacterial transgenes with profound metabolic effects." *Cell Host & Microbe* (2025). BioProject [PRJNA1258316](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1258316) (452 SRA experiments, 268 BioSamples, mouse cecal metagenome/metatranscriptome).

## Files expected in this directory

| File | Description | Source |
|---|---|---|
| `coco_gene_expr.pkl` | Per-sample gene expression matrix, indexed by real sample names (e.g. `cFT13a_S48`) | Derived from BioProject PRJNA1258316 metatranscriptomic counts, processed per the source study's `processing_metaT` pipeline |
| `gems_wol2/*.xml` | Genome-scale metabolic models: Muribaculaceae sp. CAG-873, Roseburia intestinalis, Faecalibacterium prausnitzii_C, Dubosiella newyorkensis, curated with real BSH reactions | Reconstructed with CarveMe from [WoL2 reference genomes](https://github.com/biocore/wol); see `src/gem_curation/` |
| `wol2_taxonomy/lineages.txt`, `genome_metaG/genome.tsv` | Genome-level metagenomic read counts and taxonomy lineage map, used to compute real per-taxon relative abundance | [ZarrinparLab/TRF-metaT](https://github.com/ZarrinparLab/TRF-metaT) — the source study's own code repository |
| `BSH_db_metadata.txt` | Metadata for the bile salt hydrolase (BSH) reference database used to identify and curate BSH reactions in the GEMs | Compiled during GEM curation; see `src/gem_curation/add_bsh_reactions_v2.py` |
| `coco_paper/` (external, not vendored) | CoCo's own code and data, used as the reference implementation for the cooperative-trade-off reaction-bound scaling method | [MetabioinfomicsLab/coco_paper](https://github.com/MetabioinfomicsLab/coco_paper) (Zampieri et al., 2023, *Cell Reports Methods*) — clone separately, gitignored here |

## Regenerating

1. Download raw reads from BioProject PRJNA1258316 via SRA (`src/data_acquisition/01_fetch_metagenome.sh`).
2. Process metatranscriptomic reads to a genome-level counts table (WoL2 reference, Woltka), matching the source study's own pipeline.
3. Run `src/gem_curation/add_bsh_reactions_v2.py` on the base CarveMe-reconstructed GEMs to add BSH reactions.
4. Confirm gene-identifier correspondence between the counts table and the GEMs before proceeding to `src/community_modeling/`.
5. Compute real per-taxon relative abundance from `TRF-metaT/data/genome_metaG/genome.tsv`, normalized per annotated gene feature (not raw read sum — genome/feature-catalog size varies by two orders of magnitude across taxa and biases raw sums); see `archive/exploratory/normalized_abundance.py` for the validated method.

Proteomic data (needed for a complete ecFBA implementation, see main README) was not part of the public release; a request has been sent to the source study's authors.
