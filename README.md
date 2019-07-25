## :fish: :seedling: :tada: ISB Aquaponics

Analysis code and scripts for the Aquaponics project with Jessica Day et. al.
available as:

**Negative plant-microbiome feedback limits productivity in aquaponics**<br>
Jessica A Day, Anne E Otwell, Christian Diener, Kourtney E Tams, Brad M Bebout,
Angela M Detweiler, Michael D Lee, Madeline T Scott, Wilson Ta, Monica Ha,
Shienna A Carreon, Kenny Tong, Abdirizak A Ali, Sean M Gibbons, Nitin S Baliga<br>
https://doi.org/10.1101/709162

Ask questions or report issues at https://github.com/gibbons-lab/aquaponics/issues.

## Install required software

We currently only support setups on MacOS or Linux.
You will need R ([installation instructions](https://cloud.r-project.org/)),
Rstudio ([installation instructions](https://www.rstudio.com/products/rstudio/download/))
and minimap2. To install minimap2 we recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html)
which allows you to install minimap2 with:

```bash
conda install -c bioconda minimap2
```

All analysis is performed by [mbtools](https://gibbons-lab.github.io/mbtools).
To install this package open your R console and use

```R
install.packages(c("BiocManager", "remotes"))
setRepositories(ind=1:2)
remotes::install_github("gibbons-lab/mbtools")
```

## Reproduce the study

### 0. Download or clone the repository.

You can download the repository using the green `Clone or download` on the
top right. Alternatively you can use git to clone the repository:

```bash
git clone https://github.com/gibbons-lab/aquaponics
```

### 1. Download the data and SILVA DB

*Coming soon* <br>
This is not required to reproduce the next steps since a final abundance
tables have already been placed in `data`.

### 2. Align full-length reads to the SILVA and database and count with the [EM algorithm](https://gibbons-lab.github.io/mbtools/articles/06_counting.html)

Again not required for the next steps. Everything is performed after downloading the
data and SILVA DB and running the `align.R` script. Note that this will require
large amounts of memory (~100GB) and several cores to be efficient. In case you have
less memory set `limited_memory = TRUE` in the `align.R` script. The number
of used threads will be inferred from `mc.cores` option in R. You can set it before
running the script. For instance by using

```bash
Rscript -e "options(mc.cores = 10); source('align.R')"
```

### 3. Reproduce study analyses

Each of the following steps can be run individually and out of order and
will reproduce the figures found in the manuscript.
To reproduce the steps open the top level of this repository and open any of
the `*.rmd` files mentioned below. The output of each step is linked as
well.

- Perform quality checks for tanks [notebook](tank_metrics.rmd) | [output](https://gibbons-lab.github.io/aquaponics/tank_metrics.nb.html)
- Analyze plant metrics [notebook](plant_metrics.rmd) | [output](https://gibbons-lab.github.io/aquaponics/plant_metrics.nb.html)
- Analyze general microbiome features [notebook](general.rmd) | [output](https://gibbons-lab.github.io/aquaponics/general.nb.html)
- Get associations between microbes and inocula/plant growth [notebook](associations.rmd) | [output](https://gibbons-lab.github.io/aquaponics/associations.nb.html)

## Nanopore benchmarks

Please see the [nanopore_benchmark](nanopore_benchmark) directory for details.


