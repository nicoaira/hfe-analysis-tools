# HFE analysis tools

Browser-based sequence analysis tools for the Clinical Genetics practice (UIC Barcelona): assess whether a sequenced *HFE* amplicon carries the C282Y variant.

**Live site:** https://nicoaira.github.io/hfe-analysis-tools/

1. Sequence annotation – locate the amplicon (GRCh38, *HFE* locus) and mark exons/introns (MANE Select)
2. Six-frame translation
3. Alignment to the UniProt reference protein
4. 3D structure viewer (RCSB PDB) with subsequence highlighting
5. Disulfide bridge finder
6. Codon ↔ amino acid map
7. Sanger chromatogram (`.ab1`) viewer

Static site, no build step and no backend: everything runs in the browser. `genome.json` holds chr6:26,077,000–26,108,000 (GRCh38) with the MANE Select exons of the genes in that window.

`data/` holds offline copies of UniProt Q30201 and PDB 1DE4, used automatically if UniProt or RCSB cannot be reached.

Run locally: `python3 -m http.server` in this folder, then open http://localhost:8000.
