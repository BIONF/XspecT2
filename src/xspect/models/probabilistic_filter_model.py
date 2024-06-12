"""Probabilistic filter model for sequence data"""

import json
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from slugify import slugify
import cobs_index as cobs
from xspect.file_io import get_record_iterator


class ProbabilisticFilterModel:
    """Probabilistic filter model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str,
        author_email: str,
        model_type: str,
        base_path: Path,
        fpr: float = 0.01,
        num_hashes: int = 7,
    ) -> None:
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

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the model"""
        return {
            "k": self.k,
            "model_display_name": self.model_display_name,
            "author": self.author,
            "author_email": self.author_email,
            "model_type": self.model_type,
            "display_names": self.display_names,
            "fpr": self.fpr,
            "num_hashes": self.num_hashes,
        }

    def __dict__(self) -> dict:
        """Returns a dictionary representation of the model"""
        return self.to_dict()

    def slug(self) -> str:
        """Returns a slug representation of the model"""
        return slugify(self.model_display_name + "-" + str(self.model_type))

    def fit(self, dir_path: Path, display_names: dict = None) -> None:
        """Adds filters to the model"""

        if display_names is None:
            display_names = {}

        if not isinstance(dir_path, Path):
            raise ValueError("Invalid directory path, must be a pathlib.Path object")

        if not dir_path.exists():
            raise ValueError("Directory path does not exist")

        if not dir_path.is_dir():
            raise ValueError("Directory path must be a directory")

        doclist = cobs.DocumentList()
        for file in dir_path.iterdir():
            if file.is_file() and file.suffix in [
                ".fasta",
                ".fna",
                ".fa",
                ".fastq",
                ".fq",
            ]:
                # cobs only uses the file name to the first "." as the document name
                if file.name in display_names:
                    self.display_names[file.name.split(".")[0]] = display_names[
                        file.name
                    ]
                else:
                    self.display_names[file.name.split(".")[0]] = file.stem
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

        cobs.classic_construct_list(
            doclist, str(self._get_cobs_index_path()), index_params
        )

    def _get_cobs_index_path(self):
        return self.base_path / self.slug() / "index.cobs_classic"

    def calculate_hits(
        self, sequence: Seq | SeqRecord, filter_ids: list[str] = None
    ) -> dict:
        """Calculates the hits for a sequence"""
        if isinstance(sequence, SeqRecord):
            sequence = sequence.seq

        if not isinstance(sequence, (Seq)):
            raise ValueError(
                "Invalid sequence, must be a Bio.Seq or a Bio.SeqRecord object"
            )

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        cobs_path = str(self._get_cobs_index_path())
        s = cobs.Search(cobs_path)
        r = s.search(str(sequence))
        result_dict = self._convert_cobs_result_to_dict(r)
        if filter_ids:
            return {doc: result_dict[doc] for doc in filter_ids}
        return result_dict

    def predict(
        self,
        sequence_input: (
            Seq
            | SeqRecord
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        filter_ids: list[str] = None,
    ) -> tuple[dict, dict]:
        """Returns scores for the sequence(s) based on the filters in the model"""
        if isinstance(sequence_input, (Seq, SeqRecord)):
            return ProbabilisticFilterModel.predict(self, [sequence_input], filter_ids)

        if self._is_sequence_list(sequence_input) | self._is_sequence_iterator(
            sequence_input
        ):
            # calculate scores and hits for all sequences combined,
            # by calculating hits for each sequence first
            kmer_sum = 0
            hits = {}
            for individual_sequence in sequence_input:
                individual_hits = self.calculate_hits(individual_sequence, filter_ids)
                for doc in individual_hits:
                    if doc in hits:
                        hits[doc] += individual_hits[doc]
                    else:
                        hits[doc] = individual_hits[doc]
                num_kmers = len(individual_sequence) - self.k + 1
                kmer_sum += num_kmers
            scores = {doc: round(hits[doc] / kmer_sum, 2) for doc in hits}
            return scores, hits

        if isinstance(sequence_input, Path):
            return ProbabilisticFilterModel.predict(
                self, get_record_iterator(sequence_input)
            )

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a SeqIO "
            "FastaIterator, a SeqIO FastqPhredIterator, or a Path object to a fasta/fastq file"
        )

    def filter(
        self,
        sequences: (
            Seq
            | SeqRecord
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Path
        ),
        threshold: float = 0.7,
        filter_ids: list[str] = None,
    ):
        """Filters the sequences"""

        if isinstance(sequences, (Seq, SeqRecord)):
            sequences = [sequences]

        if self._is_sequence_list(sequences) | self._is_sequence_iterator(sequences):
            filtered_sequences = {}
            for sequence in sequences:
                scores, _ = self.predict(sequence, filter_ids)
                for filter_id, score in scores.items():
                    if score >= threshold:
                        if filter_id in filtered_sequences:
                            filtered_sequences[filter_id].append(sequence)
                        else:
                            filtered_sequences[filter_id] = [sequence]

            return filtered_sequences

        if isinstance(sequences, Path):
            return self.filter(get_record_iterator(sequences), threshold, filter_ids)

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a SeqIO "
            "FastaIterator, a SeqIO FastqPhredIterator, or a Path object to a fasta/fastq file"
        )

    def save(self) -> None:
        """Saves the model to disk"""
        json_path = self.base_path / f"{self.slug()}.json"
        filter_path = self.base_path / self.slug()
        filter_path.mkdir(exist_ok=True, parents=True)

        json_object = json.dumps(self.to_dict(), indent=4)

        with open(json_path, "w", encoding="utf-8") as file:
            file.write(json_object)

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterModel":
        """Loads the model from a file"""
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
            )
            model.display_names = model_json["display_names"]
            return model

    def _convert_cobs_result_to_dict(self, cobs_result: cobs.SearchResult) -> dict:
        return {
            individual_result.doc_name: individual_result.score
            for individual_result in cobs_result
        }

    def _count_kmers(
        self,
        sequence_input: (
            Seq
            | list[Seq]
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
        ),
    ) -> int:
        """Counts the number of kmers in the sequence(s)"""
        if isinstance(sequence_input, Seq):
            return self._count_kmers([sequence_input])

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
                num_kmers = len(seq) - self.k + 1
                kmer_sum += num_kmers
            return kmer_sum

        raise ValueError(
            "Invalid sequence input, must be a Seq object, a list of Seq objects, a SeqIO "
            "FastaIterator, or a SeqIO FastqPhredIterator"
        )

    def _is_sequence_list(self, sequence_input):
        return isinstance(sequence_input, list) and all(
            isinstance(seq, (Seq, SeqRecord)) for seq in sequence_input
        )

    def _is_sequence_iterator(self, sequence_input):
        return isinstance(
            sequence_input,
            (SeqIO.FastaIO.FastaIterator, SeqIO.QualityIO.FastqPhredIterator),
        )
