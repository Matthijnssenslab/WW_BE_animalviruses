# WW_BE_animalviruses

Scripts and data to reproduce figures from our wastewater animal virus surveillance study.

## Publication
Published article: https://doi.org/10.1016/j.envint.2025.109500

**Citation**: Mustafa Karatas, Mandy Bloemen, Jill Swinnen, Inge Roukaerts, Steven Van Gucht, Marc Van Ranst, Elke Wollants, Jelle Matthijnssens,
Untapped potential of wastewater for animal and potentially zoonotic virus surveillance: Pilot study to detect non-human animal viruses in urban settings,
Environment International,
Volume 199,
2025,
109500,
ISSN 0160-4120,
https://doi.org/10.1016/j.envint.2025.109500.
(https://www.sciencedirect.com/science/article/pii/S016041202500251X)

## Repository Structure
```
├── scripts/                    # R scripts for analysis and figure generation
│   ├── Figure_1A.R            # Completeness distribution analysis
│   ├── Figure_1B.R            # Coverage analysis
│   ├── Figure_1C.R            # Sample diversity analysis  
│   ├── Figure_2.R             # Taxonomic distribution plots
│   ├── Figure_3.R             # Phylogenetic analysis
│   └── plot_coverage_zoonosis.R # Coverage plotting utilities
├── input/                     # Input data files
│   ├── data.csv              # Main dataset with virus detection results
│   └── bam_fasta/            # Reference genomes and mapping files
│       ├── *.fasta           # Viral reference sequences
│       ├── *.gff3            # Genome annotations
│       └── bam_mapped.bed    # Read mapping coordinates
└── figures/                   # Generated publication figures (PDF)
    ├── Figure1.pdf
    ├── Figure2.pdf
    ├── Figure3.pdf
    ├── Figure4.pdf
    ├── Figure5.pdf
    └── Supplementary_information.pdf
```

## Usage
1. Ensure R is installed with required packages, see scripts for packages.
2. Run scripts individually to regenerate specific figures
3. Scripts assume `data.csv` is in the working directory

## Contact
For questions about the code or data, please open an issue on GitHub or you can reach me with my email mustafa.karatas@kuleuven.be

## Acknowledgments
We gratefully acknowledge the support of the Research Foundation – Flanders (FWO) for funding Mustafa Karataş with the PhD Fellowship fundamental research grant (11P7I24N). Additionally, we thank the National Reference Center for their support and providing HEV sequences.