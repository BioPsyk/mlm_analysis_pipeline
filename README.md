# mlm_analysis_pipeline

## This is a nextflow pipeline that performs mixed model association analysis using SAIGE

### Obtaining the pipeline

The pipeline can be cloned from git using the following command:

        git clone https://github.com/BioPsyk/mlm_analysis_pipeline.git

### Usage

The pipeline can be launched from your working directory using the following command:

        nextflow run main.nf \
        --bfile <path to pruned plink files for making the GRM in step 1 of SAIGE> \
        --genotypes <path to json file containing per chromosome VCF.gz files> \
        --phenoCovFile <path to file containing the phenotype and covariates> \
        --outcome <binary/continuous> \
        --phenotype <phenotype col name in the supplied phenoCovFile>
        -w <path to where you want the nextflow work directory to be present>

Running `nextflow run main.nf -h` displays the help message

The default values for all parameters and example file formats can be seen from the</br> 
config file: `nextlow.config`

I recommend making changes here for project name on slurm etc before running the pipeline.

### Helpful tip

Nextflow recommends not launching more than one workflow from the same directory. </br>
So if you are launching the pipeline on multiple phenotypes or datasets, create new </br>
directories and use `--chdir` option in the `sbatch` command to specify a launch directory

### Expected outputs

* Per chromosome GWAS association outputs
* Concatenated GWAS association output
* Manhattan and qq plots

### Bugs and error reporting

Any questions or concerns regarding the pipeline, raise an issue on github
(or) send me an e-mail/slack message `vivek.appadurai@regionh.dk`
