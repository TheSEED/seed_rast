# SEED / RAST Code

## Overview

In 2004, [the SEED](http://theseed.org/) was created to provide consistent and accurate genome annotations across thousands of genomes and as a platform for discovering and developing de novo annotations. The SEED is a constantly updated integration of genomic data with a genome database, web front end, API and server scripts.

[RAST (Rapid Annotation using Subsystem Technology)](https://rast.nmpdr.org) is a fully-automated service for annotating complete or nearly complete bacterial and archaeal genomes. It provides high quality genome annotations for these genomes across the whole phylogenetic tree.

## About this module

This module contains the RAST code that was originally integrated with the SEED system. The code is used in both RAST and the BV-BRC genome annotation service.

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

## References

- Aziz RK, Bartels D, Best AA, DeJongh M, Disz T, Edwards RA, Formsma K, Gerdes S, Glass EM, Kubal M, Meyer F, Olsen GJ, Olson R, Osterman AL, Overbeek RA, McNeil LK, Paarmann D, Paczian T, Parrello B, Pusch GD, Reich C, Stevens R, Vassieva O, Vonstein V, Wilke A, Zagnitko O. The RAST Server: rapid annotations using subsystems technology. BMC Genomics. 2008 Feb 8;9:75. doi: 10.1186/1471-2164-9-75. PMID: 18261238; PMCID: PMC2265698.
- Overbeek R, Olson R, Pusch GD, Olsen GJ, Davis JJ, Disz T, Edwards RA, Gerdes S, Parrello B, Shukla M, Vonstein V, Wattam AR, Xia F, Stevens R. The SEED and the Rapid Annotation of microbial genomes using Subsystems Technology (RAST). Nucleic Acids Res. 2014 Jan;42(Database issue):D206-14. doi: 10.1093/nar/gkt1226. Epub 2013 Nov 29. PMID: 24293654; PMCID: PMC3965101.
- Brettin T, Davis JJ, Disz T, Edwards RA, Gerdes S, Olsen GJ, Olson R, Overbeek R, Parrello B, Pusch GD, Shukla M, Thomason JA 3rd, Stevens R, Vonstein V, Wattam AR, Xia F. RASTtk: a modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes. Sci Rep. 2015 Feb 10;5:8365. doi: 10.1038/srep08365. PMID: 25666585; PMCID: PMC4322359.
