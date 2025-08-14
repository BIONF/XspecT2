# How to use the CLI

XspecT comes with a built-in command line interface (CLI), which enables quick classifications without the need to use the web interface. The command line interface can also be used to download and train models.

After installing XspecT, a list of available commands can be viewed by running:

```bash
xspect --help
```

In general, XspecT commands will prompt you for parameters if they are not provided. However, you can also provide them directly in the command line, for example when using scripts or tools such as Slurm. Simply run the command with the `--help` option to see all available parameters.

## Model Management

At its core, XspecT uses models to classify and filter samples. These models are based on kmer indices trained on publicly available genomes as well as, possibly, a support vector machine (SVM) classifier.

To manage models, the `xspect models` command can be used. This command allows you to download, train, and view available models.

### Viewing Available Models

To view a list of available models, run:

```bash
xspect models list
```
This will show a list of all available models, separated by their type (species, genus, MLST).

### Downloading Models

To download a basic set of pre-trained models (Acinetobacter and Salonella), run:

```bash
xspect models download
```

### Model Training

Models can be trained based on data from NCBI, which is automatically downloaded and processed by XspecT.

To train a model with NCBI data, run the following command:

```bash
xspect models train ncbi
```

If you would like to train models with manually curated data from a directory, you can use:

```bash
xspect models train directory
```

Your directory should have the following structure:
```
your-directory/
├── cobs
│   ├── species1
│   │   ├── genome1.fna
│   │   ├── genome2.fna
│   │   └── ...
│   ├── species2
│   │   ├── genome1.fna
│   │   ├── genome2.fna
│   │   └── ...
│   └── ...
├── svm
│   ├── species1
│   │   ├── genome1.fna
│   │   ├── genome2.fna
│   │   └── ...
│   ├── species2
│   │   ├── genome1.fna
│   │   ├── genome2.fna
│   │   └── ...
│   └── ...
```

To train models for MLST classifications, run:

```bash
xspect models train mlst
```

## Classification

To classify samples, the command `xspect classify` can be used. This command will classify the sample based on the models available in your XspecT installation.

### Genus Classification

To classify a sample based on its genus, run the following command:

```bash
xspect classify genus
```

XspecT will prompt you for the genus and path to your sample directory.

### Species Classification

To classify a sample based on its species, run the following command:

```bash
xspect classify species
```

XspecT will prompt you for the genus and path to your sample directory.

### Sparse Sampling
XspecT uses a kmer-based approach to classify samples. This means that the entire sample is analyzed, which can be time-consuming for large samples. To speed up the analysis, you can use the `--sparse-sampling-step` option to only consider every nth kmer:

**Example**:
```bash
xspect classify species --sparse-sampling-step 10 Acinetobacter path
```

This will only consider every 10th kmer in the sample.

### Inclusion of display names
By default, the classification results show only the taxonomy ID of each species along with its corresponding score for better readability. To display the full names associated with each taxonomy ID, you can use the `--display-names` (or `-n`) option:

```bash
xspect classify species --display-names Acinetobacter path
```
The output will then be formatted as: `Taxonomy_ID - Display_Name: Score` for each species.

### MLST Classification

Samples can also be classified based on Multi-locus sequence type schemas. To MLST-classify a sample, run:

```bash
xspect classify mlst
```

## Filtering
XspecT can also be used to filter samples based on their classification results. This is useful when analyzing metagenomic samples, for example when looking at genomic bycatch.

To filter samples, the command `xspect filter` can be used. This command will filter the samples based on the specified criteria.

### Filtering by Genus

To filter samples by genus, run the following command:

```bash
xspect filter genus
```
XspecT will prompt you for the genus and path to your sample directory, as well as for a threshold to use for filtering.

### Filtering by Species
To filter samples by species, run the following command:

```bash
xspect filter species
```

You will be prompted for the genus and path to your sample directory, as well for the species to filter by and for a threshold to use for filtering. Next to normal threshold-based filtering, you can also enter a threshold of `-1` to only include contigs if the selected species is the maximum scoring species.