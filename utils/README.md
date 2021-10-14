# glyBERT scripts

All scripts and programs used to perform this analysis

---
## Installation instructions

- Packages and dependencies
- Should set add this directory to python path to import stuff from sugarTrees.py
- Update the homeDir variable in sugarTrees.py the path in which this directory was cloned

---

## Brief overview of the contents of each directory and how they were used
1. Preprocessing
2. Postprocessing

### Pre-processing
General glycan encoding and SugarBase glycans
- `sugarTrees.py`: defines classes and helper functions for encoding glycan structures from IUPAC format
- `collapseTaxonomy.py`: Map listed genus species of origin from SugarBase to defined NCBI taxids and get out labels for taxonomic origin at each rank.
- `getStats.py`: pre-process SugarBase glycans, removing improperly formatted or duplicated glycans and calculating statistics of glycan structural arrangements and compositions.
- `getEncodings.py`: generate array-form glycan encodings with 7 descriptors per monosaccharide from SugarBase glycans

Glycans on CFG microarray
- `reformat_CFG_glycans.py`: Process IUPAC representations of CFG microarray glycans, reformat to match SugarBase glycans, and get encodings for unique glycans
- `coff_cfgGlys2glyIDs.py`: Re-format lectin-glycan interaction data files from Coff et al 2020 to use keys relating to encoding/embeddings for glyBERT (v5.x glycan ID numbers)

Matched structures for immunogenic glycans (generative exploration set-up)
- `matchedTree-immuno.py`: Explore toplogical relationships between glycan strutures to find sets of glycans with different immunogenicity labels reachable by a generative process that only subsititues monosaccharides

### Post-processing
- `gly_emb_UMAP.py`: Explore and visualize whole glycans embeddings in UMAP space
- `monosacc_emb_UMAP.py`: Explore and visualize monosacccharide-level embeddings in UMAP space
- `prediction_results.R`: Process and visualize results of predictive classification
- `immunoGen_plots.py`: Process results of generative exploration and reformat to read into R
- `immGenGraphs.R`: Visualize results of initial generative exploration (fig 4)
- `immunoGen_plots_run2.py`: Process results of re-run generative exploration saving out results for more generative paths than done in the initial run
- `immGenGraphs_run2.R`: Visualize results of complete generative exploration (fig S2)
