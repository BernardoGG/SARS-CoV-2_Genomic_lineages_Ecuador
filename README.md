# SARS-CoV-2_Genomic_lineages_Ecuador
An analysis of the early transmission dynamics of SARS-CoV-2 in Ecuador inferred from genomic, epidemiological and human mobility data.

Genomic epidemiology of SARS-CoV-2 transmission lineages in Ecuador

Bernardo Gutierrez<sup>1,2</sup>, Sully Márquez3, Belén Prado-Vivar3, Mónica Becerra-Wong3, Juan José Guadalupe4, Darlan da Silva Candido1, Juan Carlos Fernandez-Cadena5, Gabriel Morey-Leon6, Rubén Armas-Gonzalez7, Derly Madeleiny Andrade-Molina5, Alfredo Bruno7,8, Domenica de Mora7, Maritza Olmedo7, Denisse Portugal7, Manuel Gonzalez7, Alberto Orlando7, Jan Felix Drexler9, Andres Moreira9, Anna-Lena Sander9, Nina Krause9, Leandro Patiño7, Andres Carrasco7, Orson Mestanza7, Jeannette Zurita10,11, Gabriela Sevillano11, Louis du Plessis1, John T. McCrone12, Josefina Coloma13, Gabriel Trueba3, Veronica Barragan3, Patricio Rojas-Silva3, Michelle Grunauer14, Moritz U.G. Kraemer1, Nuno R. Faria1,15, Marina Escalera-Zamudio1, Oliver G. Pybus1,16*, Paul Cárdenas3*

1Department of Zoology, University of Oxford, Oxford, UK
2Colegio de Ciencias Biológicas y Ambientales, Universidad San Francisco de Quito, Quito, Ecuador
3Instituto de Microbiología, Universidad San Francisco de Quito, Quito, Ecuador
4Laboratorio de Biotecnología Vegetal, Universidad San Francisco de Quito, Quito, Ecuador
5Omics Sciences Laboratory, Faculty of Medical Sciences, Universidad de Especialidades Espíritu Santo, Samborondón, Ecuador
6Faculty of Medical Sciences, Universidad de Guayaquil, Guayaquil, Ecuador
7Faculty of Sciences, Escuela Superior Politécnica del Litoral, Guayaquil, Ecuador
7Instituto Nacional de Investigación en Salud Pública, Guayaquil, Ecuador
8Universidad Agraria del Ecuador
9Charité – Universitätsmedizin Berlin, Berlin, Germany
10Facultad de Medicina, Pontificia Universidad Católica del Ecuador, Quito, Ecuador
11Unidad de Investigaciones en Biomedicina, Zurita & Zurita Laboratorios, Quito, Ecuador
12Institute of Evolutionary Biology, University of Edinburgh, Edinburgh, UK
13School of Public Health, University of California, Berkeley, USA
14Escuela de Medicina, Universidad San Francisco de Quito, Quito, Ecuador
15MRC Centre for Global Infectious Disease Analysis, J-IDEA, Imperial College London, London, UK 
16Department of Pathobiology and Population Sciences, Royal Veterinary College London, London, UK

## Repository usage and structure

The structure of this repository is shown below. It highlights the location of individual R scripts used to generate plots or analyses and input data files. Input and output data from the phylogenetic analyses are not listed in their entirety.

```
SARS-CoV-2_Genomic_lineages_Ecuador/
├── Data
│   ├── MSP_epi_data
│   ├── genetic_data
│   │   ├── GISAID_complete_highcoverage_20210107
│   │   └── metadata
│   └── mobility_data
├── Analyses
│   ├── spatial_dynamics
│   └── phylogenetics
│       ├── ML_SARS2_EC
│       └── BEAST_SARS2_EC
├── Plots
│   └── SARS2_phylo_plots.R
└── README.md
```
