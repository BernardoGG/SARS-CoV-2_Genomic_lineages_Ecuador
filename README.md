# SARS-CoV-2_Genomic_lineages_Ecuador
An analysis of the early transmission dynamics of SARS-CoV-2 in Ecuador inferred from genomic, epidemiological and human mobility data.

Bernardo Gutierrez<sup>1,2</sup>, Sully Márquez<sup>3</sup>, Belén Prado-Vivar<sup>3</sup>, Mónica Becerra-Wong<sup>3</sup>, Juan José Guadalupe<sup>4</sup>, Darlan da Silva Candido<sup>1</sup>, Juan Carlos Fernandez-Cadena<sup>5</sup>, Gabriel Morey-Leon<sup>6</sup>, Rubén Armas-Gonzalez<sup>7</sup>, Derly Madeleiny Andrade-Molina<sup>5</sup>, Alfredo Bruno<sup>8,9</sup>, Domenica de Mora<sup>8</sup>, Maritza Olmedo<sup>8</sup>, Denisse Portugal<sup>8</sup>, Manuel Gonzalez<sup>8</sup>, Alberto Orlando<sup>8</sup>, Jan Felix Drexler<sup>10</sup>, Andres Moreira<sup>10</sup>, Anna-Lena Sander<sup>10</sup>, Nina Krause<sup>10</sup>, Leandro Patiño<sup>8</sup>, Andres Carrasco<sup>8</sup>, Orson Mestanza<sup>8</sup>, Jeannette Zurita<sup>11,12</sup>, Gabriela Sevillano<sup>12</sup>, Louis du Plessis<sup>1</sup>, John T. McCrone<sup>13</sup>, Josefina Coloma<sup>14</sup>, Gabriel Trueba<sup>3</sup>, Veronica Barragan<sup>3</sup>, Patricio Rojas-Silva<sup>3</sup>, Michelle Grunauer<sup>15</sup>, Moritz U.G. Kraemer<sup>1</sup>, Nuno R. Faria<sup>1,16</sup>, Marina Escalera-Zamudio<sup>1</sup>, Oliver G. Pybus<sup>1,17*</sup>, Paul Cárdenas<sup>3*</sup>

<sup>1</sup>Department of Zoology, University of Oxford, Oxford, UK
<sup>2</sup>Colegio de Ciencias Biológicas y Ambientales, Universidad San Francisco de Quito, Quito, Ecuador
<sup>3</sup>Instituto de Microbiología, Universidad San Francisco de Quito, Quito, Ecuador
<sup>4</sup>Laboratorio de Biotecnología Vegetal, Universidad San Francisco de Quito, Quito, Ecuador
<sup>5</sup>Omics Sciences Laboratory, Faculty of Medical Sciences, Universidad de Especialidades Espíritu Santo, Samborondón, Ecuador
<sup>6</sup>Faculty of Medical Sciences, Universidad de Guayaquil, Guayaquil, Ecuador
<sup>7</sup>Faculty of Sciences, Escuela Superior Politécnica del Litoral, Guayaquil, Ecuador
<sup>8</sup>Instituto Nacional de Investigación en Salud Pública, Guayaquil, Ecuador
<sup>9</sup>Universidad Agraria del Ecuador
<sup>10</sup>Charité – Universitätsmedizin Berlin, Berlin, Germany
<sup>11</sup>Facultad de Medicina, Pontificia Universidad Católica del Ecuador, Quito, Ecuador
<sup>12</sup>Unidad de Investigaciones en Biomedicina, Zurita & Zurita Laboratorios, Quito, Ecuador
<sup>13</sup>Institute of Evolutionary Biology, University of Edinburgh, Edinburgh, UK
<sup>14</sup>School of Public Health, University of California, Berkeley, USA
<sup>15</sup>Escuela de Medicina, Universidad San Francisco de Quito, Quito, Ecuador
<sup>16</sup>MRC Centre for Global Infectious Disease Analysis, J-IDEA, Imperial College London, London, UK 
<sup>17</sup>Department of Pathobiology and Population Sciences, Royal Veterinary College London, London, UK

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
