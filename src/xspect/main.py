"""Project CLI"""

from pathlib import Path
from uuid import uuid4
from importlib import import_module
import click
from xspect.model_management import get_models

# inline imports lead to "invalid name" issues
# pylint: disable=invalid-name


@click.group()
@click.version_option()
def cli():
    """XspecT CLI."""


@cli.command()
def web():
    """Open the XspecT web application."""
    app = import_module("xspect.web").app
    run = import_module("uvicorn").run

    run(app, host="0.0.0.0", port=8000)


# # # # # # # # # # # # # # #
# Model management commands #
# # # # # # # # # # # # # # #
@cli.group()
def models():
    """Model management commands."""


@models.command(
    help="Download models from the internet.",
)
def download():
    """Download models."""
    click.echo("Downloading models, this may take a while...")
    download_test_models = import_module("xspect.download_models").download_test_models

    download_test_models(
        "https://assets.adrianromberg.com/science/xspect-models-07-08-2025.zip"
    )


@models.command(
    name="list",
    help="List all models in the model directory.",
)
def list_models():
    """List models."""
    available_models = get_models()
    if not available_models:
        click.echo("No models found.")
        return
    click.echo("Models found:")
    click.echo("--------------")
    for model_type, names in available_models.items():
        if not names:
            continue
        click.echo(f"  {model_type}:")
        for name in names:
            click.echo(f"    - {name}")


@models.group()
def train():
    """Train models."""


@train.command(
    name="ncbi",
    help="Train a species and a genus model based on NCBI data.",
)
@click.option("-g", "--genus", "model_genus", prompt=True)
@click.option("--svm_steps", type=int, default=1)
@click.option(
    "--author",
    help="Author of the model.",
    default=None,
)
@click.option(
    "--author-email",
    help="Email of the author.",
    default=None,
)
def train_ncbi(model_genus, svm_steps, author, author_email):
    """Train a species and a genus model based on NCBI data."""
    click.echo(f"Training {model_genus} species and genus metagenome model.")
    try:
        train_from_ncbi = import_module("xspect.train").train_from_ncbi

        train_from_ncbi(model_genus, svm_steps, author, author_email)
    except ValueError as e:
        click.echo(f"Error: {e}")
        return
    click.echo(f"Training of {model_genus} model finished.")


@train.command(
    name="directory",
    help="Train a species (and possibly a genus) model based on local data.",
)
@click.option("-g", "--genus", "model_genus", prompt=True)
@click.option(
    "-i",
    "--input-path",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
)
@click.option(
    "--meta",
    is_flag=True,
    help="Train a metagenome model for the genus.",
    default=True,
)
@click.option(
    "--svm-steps",
    type=int,
    help="SVM Sparse sampling step size (e. g. only every 500th kmer for step=500).",
    default=1,
)
@click.option(
    "--author",
    help="Author of the model.",
    default=None,
)
@click.option(
    "--author-email",
    help="Email of the author.",
    default=None,
)
def train_directory(model_genus, input_path, svm_steps, meta, author, author_email):
    """Train a model based on data from a directory for a given genus."""
    click.echo(f"Training {model_genus} model with {svm_steps} SVM steps.")
    train_from_directory = import_module("xspect.train").train_from_directory

    train_from_directory(
        model_genus,
        Path(input_path),
        svm_step=svm_steps,
        meta=meta,
        author=author,
        author_email=author_email,
    )


@train.command(
    name="mlst",
    help="Train a MLST model based on PubMLST data.",
)
@click.option(
    "-c",
    "--choose_schemes",
    is_flag=True,
    help="Choose your own schemes."
    "Default setting is Oxford and Pasteur scheme of A.baumannii.",
)
def train_mlst(choose_schemes):
    """Download alleles and train bloom filters."""
    click.echo("Updating alleles")
    mlst_helper = import_module("xspect.mlst_feature.mlst_helper")
    pick_scheme = mlst_helper.pick_scheme

    pub_mlst_handler = import_module("xspect.mlst_feature.pub_mlst_handler")
    PubMLSTHandler = pub_mlst_handler.PubMLSTHandler

    probabilistic_filter_mlst_model = import_module(
        "xspect.models.probabilistic_filter_mlst_model"
    )
    ProbabilisticFilterMlstSchemeModel = (
        probabilistic_filter_mlst_model.ProbabilisticFilterMlstSchemeModel
    )

    definitions = import_module("xspect.definitions")
    get_xspect_model_path = definitions.get_xspect_model_path

    handler = PubMLSTHandler()
    handler.download_alleles(choose_schemes)
    click.echo("Download finished")
    scheme_path = pick_scheme(handler.get_scheme_paths())
    species_name = str(scheme_path).split("/")[-2]
    scheme_name = str(scheme_path).rsplit("/", maxsplit=1)[-1]
    scheme_url = handler.scheme_mapping[str(scheme_path)]
    model = ProbabilisticFilterMlstSchemeModel(
        31, f"{species_name}:{scheme_name}", get_xspect_model_path(), scheme_url
    )
    click.echo("Creating mlst model")
    model.fit(scheme_path)
    model.save()
    click.echo(f"Saved at {model.cobs_path}")


# # # # # # # # # # # # # # #
# Classification commands   #
# # # # # # # # # # # # # # #
@cli.group(
    name="classify",
    help="Classify sequences using XspecT models.",
)
def classify_seqs():
    """Classification commands."""


@classify_seqs.command(
    name="genus",
    help="Classify samples using a genus model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to classify.",
    type=click.Choice(get_models().get("Genus", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"result_{uuid4()}.json",
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
def classify_genus(model_genus, input_path, output_path, sparse_sampling_step):
    """Classify samples using a genus model."""
    click.echo("Classifying...")
    classify = import_module("xspect.classify")

    classify.classify_genus(
        model_genus, Path(input_path), Path(output_path), sparse_sampling_step
    )


@classify_seqs.command(
    name="species",
    help="Classify samples using a species model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to classify.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"result_{uuid4()}.json",
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
@click.option(
    "-n",
    "--display-names",
    help="Includes the display names next to taxonomy-IDs.",
    is_flag=True,
)
def classify_species(
    model_genus, input_path, output_path, sparse_sampling_step, display_names
):
    """Classify samples using a species model."""
    click.echo("Classifying...")
    classify = import_module("xspect.classify")

    classify.classify_species(
        model_genus,
        Path(input_path),
        Path(output_path),
        sparse_sampling_step,
        display_names,
    )


@classify_seqs.command(
    name="mlst",
    help="Classify samples using a MLST model.",
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA-file for mlst identification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    default=Path(".") / f"MLST_result_{uuid4()}.json",
)
@click.option(
    "-l", "--limit", is_flag=True, help="Limit the output to 5 results for each locus."
)
def classify_mlst(input_path, output_path, limit):
    """MLST classify a sample."""
    click.echo("Classifying...")
    classify = import_module("xspect.classify")

    classify.classify_mlst(Path(input_path), Path(output_path), limit)


# # # # # # # # # # # # # # #
# Filtering commands        #
# # # # # # # # # # # # # # #
@cli.group(
    name="filter",
    help="Filter sequences using XspecT models.",
)
def filter_seqs():
    """Filter commands."""


@filter_seqs.command(
    name="genus",
    help="Filter sequences using a genus model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to use for filtering.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    prompt=True,
    default=Path(".") / f"genus_filtered_{uuid4()}.fasta",
)
@click.option(
    "--classification-output-path",
    help="Optional path to the classification output file.",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "-t",
    "--threshold",
    type=click.FloatRange(0, 1),
    help="Threshold for filtering (default: 0.7).",
    default=0.7,
    prompt=True,
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for '--sparse-sampling-step 500').",
    default=1,
)
def filter_genus(
    model_genus,
    input_path,
    output_path,
    classification_output_path,
    threshold,
    sparse_sampling_step,
):
    """Filter samples using a genus model."""
    click.echo("Filtering...")
    filter_sequences = import_module("xspect.filter_sequences")

    filter_sequences.filter_genus(
        model_genus,
        Path(input_path),
        Path(output_path),
        threshold,
        Path(classification_output_path) if classification_output_path else None,
        sparse_sampling_step=sparse_sampling_step,
    )


@filter_seqs.command(
    name="species",
    help="Filter sequences using a species model.",
)
@click.option(
    "-g",
    "--genus",
    "model_genus",
    help="Genus of the model to use for filtering.",
    type=click.Choice(get_models().get("Species", [])),
    prompt=True,
)
@click.option(
    "-s",
    "--species",
    "model_species",
    help="Species of the model to filter for.",
)
@click.option(
    "-i",
    "--input-path",
    help="Path to FASTA or FASTQ file for classification.",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    prompt=True,
    default=Path("."),
)
@click.option(
    "-o",
    "--output-path",
    help="Path to the output file.",
    type=click.Path(dir_okay=False, file_okay=True),
    prompt=True,
    default=Path(".") / f"species_filtered_{uuid4()}.fasta",
)
@click.option(
    "--classification-output-path",
    help="Optional path to the classification output file.",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "-t",
    "--threshold",
    type=float,
    help="Threshold for filtering (default: 0.7). Use -1 to filter for the highest scoring "
    "species.",
    default=0.7,
    prompt=True,
)
@click.option(
    "--sparse-sampling-step",
    type=int,
    help="Sparse sampling step (e. g. only every 500th kmer for "
    "'--sparse-sampling-step 500').",
    default=1,
)
def filter_species(
    model_genus,
    model_species,
    input_path,
    output_path,
    threshold,
    classification_output_path,
    sparse_sampling_step,
):
    """Filter a sample using the species model."""

    if threshold != -1 and (threshold < 0 or threshold > 1):
        raise click.BadParameter(
            "Threshold must be between 0 and 1, or -1 for filtering by the highest "
            "scoring species."
        )

    get_model_metadata = import_module("xspect.model_management").get_model_metadata

    available_species = get_model_metadata(f"{model_genus}-species")["display_names"]
    available_species = {
        id: name.replace(f"{model_genus} ", "")
        for id, name in available_species.items()
    }
    if not model_species:
        sorted_available_species = sorted(available_species.values())
        model_species = click.prompt(
            f"Please enter the species name: {model_genus}",
            type=click.Choice(sorted_available_species, case_sensitive=False),
        )
    if model_species not in available_species.values():
        raise click.BadParameter(
            f"Species '{model_species}' not found in the {model_genus} species model."
        )

    # get the species ID from the name
    model_species = [
        id
        for id, name in available_species.items()
        if name.lower() == model_species.lower()
    ][0]

    click.echo("Filtering...")
    filter_sequences = import_module("xspect.filter_sequences")

    filter_sequences.filter_species(
        model_genus,
        model_species,
        Path(input_path),
        Path(output_path),
        threshold,
        Path(classification_output_path) if classification_output_path else None,
        sparse_sampling_step=sparse_sampling_step,
    )


if __name__ == "__main__":
    cli()
