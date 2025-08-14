"""Probabilistic filter model for sequence data"""

import json
from math import ceil
from pathlib import Path
from typing import Any
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from slugify import slugify
import cobs_index as cobs
from xspect.definitions import fasta_endings, fastq_endings
from xspect.file_io import get_record_iterator
from xspect.models.result import ModelResult


class ProbabilisticFilterModel:
    """Probabilistic filter model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str | None,
        author_email: str | None,
        model_type: str,
        base_path: Path,
        fpr: float = 0.01,
        num_hashes: int = 7,
        training_accessions: dict[str, list[str]] | None = None,
    ) -> None:
        """
        Initializes the probabilistic filter model.

        This method sets up the model with the specified parameters, including the k-mer size,
        display name, author information, model type, base path for storage, false positive rate,
        number of hashes, and training accessions.

        Args:
            k (int): The size of the k-mers to be used in the model.
            model_display_name (str): The display name of the model.
            author (str | None): The name of the author of the model.
            author_email (str | None): The email of the author of the model.
            model_type (str): The type of the model.
            base_path (Path): The base path where the model will be stored.
            fpr (float): The false positive rate for the model. Default is 0.01.
            num_hashes (int): The number of hashes to use in the model. Default is 7.
            training_accessions (dict[str, list[str]] | None): A dictionary mapping filter IDs to
                lists of accession numbers used for training the model. Default is None.
        """
        if k < 1:
            raise ValueError("Invalid k value, must be greater than 0")
        if not model_display_name:
            raise ValueError("Invalid filter display name, must be a non-empty string")
        if not model_type:
            raise ValueError("Invalid filter type, must be a non-empty string")
        if not isinstance(base_path, Path):
            raise ValueError("Invalid base path, must be a pathlib.Path object")

        self.k = k
        self.model_display_name = model_display_name
        self.author = author
        self.author_email = author_email
        self.model_type = model_type
        self.base_path = base_path
        self.display_names = {}
        self.fpr = fpr
        self.num_hashes = num_hashes
        self.index = None
        self.training_accessions = training_accessions

    def get_cobs_index_path(self) -> str:
        """
        Returns the path to the cobs inde

        This method constructs the path where the cobs index file will be stored,
        based on the model's slug and the base path.

        Returns:
            str: The path to the cobs index file.
        """
        return str(self.base_path / self.slug() / "index.cobs_classic")

    def to_dict(self) -> dict:
        """
        Returns a dictionary representation of the model

        This method includes all relevant attributes of the model, such as k-mer size,
        display name, author information, model type, and other parameters.

        Returns:
            dict: A dictionary containing the model's attributes.
        """
        return {
            "model_slug": self.slug(),
            "k": self.k,
            "model_display_name": self.model_display_name,
            "author": self.author,
            "author_email": self.author_email,
            "model_type": self.model_type,
            "model_class": self.__class__.__name__,
            "display_names": self.display_names,
            "fpr": self.fpr,
            "num_hashes": self.num_hashes,
            "training_accessions": self.training_accessions,
        }

    def slug(self) -> str:
        """
        Returns a slug representation of the model

        This method generates a slug based on the model's display name and type,
        which can be used for file naming or identification purposes.

        Returns:
            str: A slug representation of the model.
        """
        return slugify(self.model_display_name + "-" + str(self.model_type))

    def fit(
        self,
        dir_path: Path,
        display_names: dict | None = None,
        training_accessions: dict[str, list[str]] | None = None,
    ) -> None:
        """
        Adds filters to the model

        This method constructs the model's index from sequence files in the specified directory.
        It reads files with specified extensions (fasta and fastq), constructs a document list,
        and builds a cobs index for efficient searching.

        Args:
            dir_path (Path): The directory containing sequence files to be indexed.
            display_names (dict | None): A dictionary mapping file names to display names.
                If None, uses file names as display names.
            training_accessions (dict[str, list[str]] | None): A dictionary mapping filter IDs to
                lists of accession numbers used for training the model. If None, no training
                accessions are set.
        Raises:
            ValueError: If the directory path is invalid, does not exist, or is not a directory.
        """

        if display_names is None:
            display_names = {}

        if not isinstance(dir_path, Path):
            raise ValueError("Invalid directory path, must be a pathlib.Path object")

        if not dir_path.exists():
            raise ValueError("Directory path does not exist")

        if not dir_path.is_dir():
            raise ValueError("Directory path must be a directory")

        self.training_accessions = training_accessions

        doclist = cobs.DocumentList()
        for file in dir_path.iterdir():
            if file.is_file() and file.suffix[1:] in fasta_endings + fastq_endings:
                # cobs only uses the file name to the first "." as the document name
                if file.stem in display_names:
                    self.display_names[file.stem.split(".")[0]] = display_names[
                        file.stem
                    ]
                else:
                    self.display_names[file.stem.split(".")[0]] = file.stem
                doclist.add(str(file))

        if len(doclist) == 0:
            raise ValueError(
                "No valid files found in directory. Must be fasta or fastq"
            )

        index_params = cobs.ClassicIndexParameters()
        index_params.term_size = self.k
        index_params.num_hashes = self.num_hashes
        index_params.false_positive_rate = self.fpr
        index_params.clobber = True

        cobs.classic_construct_list(doclist, self.get_cobs_index_path(), index_params)

        self.index = cobs.Search(self.get_cobs_index_path(), True)

    def calculate_hits(
        self, sequence: Seq, filter_ids: list[str] | None = None, step: int = 1
    ) -> dict:
        """
        Calculates the hits for a sequence

        This method searches the model's index for the given sequence and returns a dictionary
        of filter IDs and their corresponding scores. If filter_ids is provided, it filters the
        results to only include those IDs.

        Args:
            sequence (Seq): The sequence to search for in the model's index.
            filter_ids (list[str] | None): A list of filter IDs to filter the results. If None,
                all results are returned.
            step (int): The step size for the k-mer search. Default is 1.

        Returns:
            dict: A dictionary where keys are filter IDs and values are scores for the sequence.

        Raises:
            ValueError: If the sequence is not a valid Bio.Seq or Bio.SeqRecord object,
                        if the sequence length is not greater than k, or if the input is invalid.
        """
        if not isinstance(sequence, (Seq)):
            raise ValueError(
                "Invalid sequence, must be a Bio.Seq or a Bio.SeqRecord object"
            )

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        r = self.index.search(str(sequence), step=step)
        result_dict = self._convert_cobs_result_to_dict(r)
        if filter_ids:
            return {doc: result_dict[doc] for doc in filter_ids}
        return result_dict

    def predict(
        self,
        sequence_input: (
            SeqRecord
            | list[SeqRecord]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        filter_ids: list[str] = None,
        step: int = 1,
        display_name: bool = False,
    ) -> ModelResult:
        """
        Returns a model result object for the sequence(s) based on the filters in the model

        This method processes the input sequence(s) and calculates hits against the model's index.
        It supports various input types, including single sequences, lists of sequences,
        SeqIO iterators, and file paths. The results are returned as a ModelResult object.

        Args:
            sequence_input (SeqRecord | list[SeqRecord] | SeqIO.FastaIO.FastaIterator |
                            SeqIO.QualityIO.FastqPhredIterator | Path):
                The input sequence(s) to be processed. Can be a single SeqRecord, a list of
                SeqRecords, a SeqIO iterator, or a Path to a fasta/fastq file.
            filter_ids (list[str]): A list of filter IDs to filter the results. If None,
                all results are returned.
            step (int): The step size for the k-mer search. Default is 1.
            display_name (bool): Includes a display name for each tax_ID.

        Returns:
            ModelResult: An object containing the hits for each sequence, the number of kmers,
                         and the sparse sampling step.

        Raises:
            ValueError: If the input sequence is not valid, or if it is not a Seq object,
                        a list of Seq objects, a SeqIO iterator, or a Path object to a fasta/fastq
                        file.
        """
        if isinstance(sequence_input, (SeqRecord)):
            return ProbabilisticFilterModel.predict(
                self, [sequence_input], filter_ids, step=step, display_name=display_name
            )

        if self._is_sequence_list(sequence_input) | self._is_sequence_iterator(
            sequence_input
        ):
            hits = {}
            num_kmers = {}
            for individual_sequence in sequence_input:
                individual_hits = self.calculate_hits(
                    individual_sequence.seq, filter_ids, step=step
                )
                num_kmers[individual_sequence.id] = self._count_kmers(
                    individual_sequence, step=step
                )
                if display_name:
                    individual_hits.update(
                        {
                            f"{key} -{self.display_names.get(key, 'Unknown').replace( 
                                self.model_display_name, '', 1)}": individual_hits.pop(
                                key
                            )
                            for key in list(individual_hits.keys())
                        }
                    )
                hits[individual_sequence.id] = individual_hits
            return ModelResult(self.slug(), hits, num_kmers, sparse_sampling_step=step)

        if isinstance(sequence_input, Path):
            return ProbabilisticFilterModel.predict(
                self,
                get_record_iterator(sequence_input),
                step=step,
                display_name=display_name,
            )

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, a SeqIO FastqPhredIterator, or a Path object to a"
            " fasta/fastq file"
        )

    def save(self) -> None:
        """
        Saves the model to disk

        This method serializes the model's attributes to a JSON file and creates a directory
        for the model based on its slug. The JSON file contains all relevant information about
        the model, including k-mer size, display name, author information, model type, and
        other parameters. The directory structure is created if it does not already exist.
        """
        json_path = self.base_path / f"{self.slug()}.json"
        filter_path = self.base_path / self.slug()
        filter_path.mkdir(exist_ok=True, parents=True)

        json_object = json.dumps(self.to_dict(), indent=4)

        with open(json_path, "w", encoding="utf-8") as file:
            file.write(json_object)

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterModel":
        """
        Loads the model from a file

        This static method reads a JSON file containing the model's attributes and constructs
        a ProbabilisticFilterModel object. It also checks for the existence of the cobs index file
        and initializes the index if it exists.

        Args:
            path (Path): The path to the JSON file containing the model's attributes.

        Returns:
            ProbabilisticFilterModel: An instance of the ProbabilisticFilterModel class
            initialized with the attributes from the JSON file.

        Raises:
            FileNotFoundError: If the JSON file or the cobs index file does not exist.
        """
        with open(path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticFilterModel(
                model_json["k"],
                model_json["model_display_name"],
                model_json["author"],
                model_json["author_email"],
                model_json["model_type"],
                path.parent,
                model_json["fpr"],
                model_json["num_hashes"],
                model_json["training_accessions"],
            )
            model.display_names = model_json["display_names"]

            p = model.get_cobs_index_path()
            if not Path(p).exists():
                raise FileNotFoundError(f"Index file not found at {p}")
            model.index = cobs.Search(p, True)

            return model

    def _convert_cobs_result_to_dict(self, cobs_result: cobs.SearchResult) -> dict:
        """
        Converts a cobs SearchResult to a dictionary

        This method takes a cobs SearchResult object and converts it into a dictionary
        where the keys are document names and the values are their corresponding scores.

        Args:
            cobs_result (cobs.SearchResult): The result object from a cobs search.

        Returns:
            dict: A dictionary mapping document names to their scores.
        """
        return {
            individual_result.doc_name: individual_result.score
            for individual_result in cobs_result
        }

    def _count_kmers(
        self,
        sequence_input: (
            Seq
            | SeqRecord
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
        ),
        step: int = 1,
    ) -> int:
        """
        Counts the number of kmers in the sequence(s)

        This method calculates the number of k-mers in a given sequence or list of sequences.
        It supports various input types, including single sequences, SeqRecords, lists of sequences,
        and SeqIO iterators. The step size for the k-mer search can be specified.

        Args:
            sequence_input (Seq | SeqRecord | list[Seq] | SeqIO.FastaIO.FastaIterator |
                            SeqIO.QualityIO.FastqPhredIterator):
                The input sequence(s) to count k-mers in. Can be a single Seq, a SeqRecord,
                a list of Seq objects, or a SeqIO iterator.
            step (int): The step size for the k-mer search. Default is 1.

        Returns:
            int: The total number of k-mers in the input sequence(s).

        Raises:
            ValueError: If the input sequence is not valid, or if it is not a Seq object,
                        a SeqRecord, a list of Seq objects, or a SeqIO iterator.
        """
        if isinstance(sequence_input, Seq):
            return self._count_kmers([sequence_input], step=step)

        if isinstance(sequence_input, SeqRecord):
            return self._count_kmers(sequence_input.seq, step=step)

        is_sequence_list = isinstance(sequence_input, list) and all(
            isinstance(seq, Seq) for seq in sequence_input
        )
        is_iterator = isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        )

        if is_sequence_list | is_iterator:
            kmer_sum = 0
            for individual_sequence in sequence_input:
                # we need to look specifically at .seq for SeqIO iterators
                seq = individual_sequence.seq if is_iterator else individual_sequence
                num_kmers = ceil((len(seq) - self.k + 1) / step)
                kmer_sum += num_kmers
            return kmer_sum

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a"
            " SeqIO FastaIterator, or a SeqIO FastqPhredIterator"
        )

    def _is_sequence_list(self, sequence_input: Any) -> bool:
        """
        Checks if the input is a list of SeqRecord objects

        This method verifies if the input is a list and that all elements in the list
        are instances of SeqRecord. This is useful for ensuring that the input is a valid
        collection of sequence records.

        Args:
            sequence_input (Any): The input to check.

        Returns:
            bool: True if the input is a list of SeqRecord objects, False otherwise.
        """
        return isinstance(sequence_input, list) and all(
            isinstance(seq, (SeqRecord)) for seq in sequence_input
        )

    def _is_sequence_iterator(self, sequence_input: Any) -> bool:
        """
        Checks if the input is a SeqIO iterator

        This method verifies if the input is an instance of a SeqIO iterator, such as
        FastaIterator or FastqPhredIterator. This is useful for ensuring that the input
        is a valid sequence iterator that can be processed by the model.

        Args:
            sequence_input (Any): The input to check.

        Returns:
            bool: True if the input is a SeqIO iterator, False otherwise.
        """
        return isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        )
