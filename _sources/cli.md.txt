# How to use the CLI

XspecT comes with a built-in command line interface (CLI), which enables quick classifications without the need to use the web interface. The command line interface can also be used to download and train models.

After installing XspecT, a list of available commands can be viewed by running:

```bash
xspect --help
```

## Model downloads

A basic set of pre-trained models (Acinetobacter and Salonella) can be downloaded using the following command:

```bash
xspect download-models
```

For the moment, it is not possible to specify exactly which models should be downloaded.

## Classification

To classify samples, the command

```bash
xspect classify-species GENUS PATH
```

can be used, when `GENUS` refers to the NCBI genus name of your sample and `PATH` refers to the path to your sample *directory*. This command will classify the species of your sample within the given genus.

The following options are available:

```bash
-m, --meta / --no-meta  Metagenome classification.
-s, --step INTEGER      Sparse sampling step size (e. g. only every 500th
                        kmer for step=500).
--help                  Show this message and exit.
```

To speed up the analysis, only every nth kmer can be considered ("sparse sampling"). For example, to only consider every 10th kmer, run:

```bash
xspect classify-species -s 10 Acinetobacter path
```

### Metagenome Mode

To analyze a sample in metagenome mode, the `-m`/`--meta` (`--no-meta`) option can be used:

```bash
xspect classify-species -m Acinetobacter path
```

Compared to normal XspecT species classification, this mode first identifies reads belonging to the given genus and continues classification only with the resulting reads, It is thus more suitable for metagenomic samples as the resulting runtime is decreased.

### MLST Classification

Samples can also be classified based on Multi-locus sequence type schemas. To MLST-classify a sample, run:

```bash
xspect classify-mlst -p path
```

## Model Training

Models can be trained based on data from NCBI, which is automatically downloaded and processed by XspecT.

To train a model, run the following command:

```bash
xspect train-species your-ncbi-genus
```

`you-ncbi-genus` can be a genus name from NCBI or an NCBI taxonomy ID.

To train models for MLST classifications, run:

```bash
xspect train-mlst
```