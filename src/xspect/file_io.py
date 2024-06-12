"""
File IO module.
"""

import os
from pathlib import Path
import zipfile
from Bio import SeqIO


def check_folder_structure():
    """Checks the folder structure and creates new folders if needed."""
    # Create list of all folder paths.
    root_path = Path(os.getcwd())
    filter_path = root_path / "filter"
    meta_path = root_path / "genus_metadata"
    filter_folder_names = [
        "array_sizes",
        "Metagenomes",
        "species_names",
        "translation_dicts",
    ]
    folder_paths = [filter_path, meta_path]
    for filter_folder_name in filter_folder_names:
        filter_folder_path = filter_path / filter_folder_name
        folder_paths.append(filter_folder_path)

    # Check if folders exist. If not create them.
    for folder_path in folder_paths:
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)


def delete_zip_files(dir_path):
    """Delete all zip files in the given directory."""
    files = os.listdir(dir_path)
    for file in files:
        if zipfile.is_zipfile(file):
            file_path = dir_path / str(file)
            os.remove(file_path)


def extract_zip(zip_path, unzipped_path):
    """Extracts all files from a directory with zip files."""
    # Make new directory.
    os.mkdir(unzipped_path)

    file_names = os.listdir(zip_path)
    for file in file_names:
        file_path = zip_path / file
        if zipfile.is_zipfile(file_path):
            with zipfile.ZipFile(file_path) as item:
                directory = unzipped_path / file.replace(".zip", "")
                item.extractall(directory)


def concatenate_meta(path: Path, genus: str):
    """Concatenates all species files to one fasta file.

    :param path: Path to the directory with the concatenated fasta files.
    :type path: Path
    :param genus: Genus name.
    :type genus: str
    """
    files_path = path / "concatenate"
    fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
    meta_path = path / (genus + ".fasta")
    files = os.listdir(files_path)

    with open(meta_path, "w", encoding="utf-8") as meta_file:
        # Write the header.
        meta_header = f">{genus} metagenome\n"
        meta_file.write(meta_header)

        # Open each concatenated species file and write the sequence in the meta file.
        for file in files:
            file_ending = str(file).rsplit(".", maxsplit=1)[-1]
            if file_ending in fasta_endings:
                with open(
                    (files_path / str(file)), "r", encoding="utf-8"
                ) as species_file:
                    for line in species_file:
                        if line[0] != ">":
                            meta_file.write(line.replace("\n", ""))


def get_record_iterator(file_path: Path):
    """Returns a record iterator for a fasta or fastq file."""
    if not isinstance(file_path, Path):
        raise ValueError("Path must be a Path object")

    if not file_path.exists():
        raise ValueError("File does not exist")

    if not file_path.is_file():
        raise ValueError("Path must be a file")

    if file_path.suffix in [".fasta", ".fna", ".fa"]:
        return SeqIO.parse(file_path, "fasta")

    if file_path.suffix in [".fastq", ".fq"]:
        return SeqIO.parse(file_path, "fastq")

    raise ValueError("Invalid file format, must be a fasta or fastq file")
