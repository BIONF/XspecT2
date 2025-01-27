"""Base probabilistic filter model for sequence data"""

# pylint: disable=no-name-in-module, too-many-instance-attributes

import json
from math import ceil
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rbloom import Bloom
from xxhash import xxh3_64_intdigest
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel
from xspect.file_io import get_record_iterator


class ProbabilisticSingleFilterModel(ProbabilisticFilterModel):
    """Base probabilistic filter model for sequence data"""

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str,
        author_email: str,
        model_type: str,
        base_path: Path,
        fpr: float = 0.01,
    ) -> None:
        super().__init__(
            k=k,
            model_display_name=model_display_name,
            author=author,
            author_email=author_email,
            model_type=model_type,
            base_path=base_path,
            fpr=fpr,
            num_hashes=1,
        )
        self.bf = None

    def fit(self, file_path: Path, display_name: str) -> None:
        """Fit the cobs classic index to the sequences and labels"""
        # estimate number of kmers
        total_length = 0
        for record in get_record_iterator(file_path):
            total_length += len(record.seq)
        num_kmers = total_length - self.k + 1

        self.bf = Bloom(num_kmers, self.fpr, hash_func=xxh3_64_intdigest)
        for record in get_record_iterator(file_path):
            for kmer in self._generate_kmers(record.seq):
                self.bf.add(kmer)
        self.display_names[file_path.stem] = display_name

        bloom_path = self.base_path / self.slug() / "filter.bloom"
        bloom_path.parent.mkdir(parents=True, exist_ok=True)
        self.bf.save(str(bloom_path))

    def calculate_hits(
        self, sequence: Seq | SeqRecord, filter_ids=None, step: int = 1
    ) -> dict:
        """Calculate the hits for the sequence"""
        if isinstance(sequence, SeqRecord):
            sequence = sequence.seq

        if not isinstance(sequence, Seq):
            raise ValueError("Invalid sequence, must be a Bio.Seq object")

        if not len(sequence) > self.k:
            raise ValueError("Invalid sequence, must be longer than k")

        num_hits = sum(
            1 for kmer in self._generate_kmers(sequence, step=step) if kmer in self.bf
        )
        return {next(iter(self.display_names)): num_hits}

    @staticmethod
    def load(path: Path) -> "ProbabilisticSingleFilterModel":
        """Load the model from disk"""
        with open(path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticSingleFilterModel(
                model_json["k"],
                model_json["model_display_name"],
                model_json["author"],
                model_json["author_email"],
                model_json["model_type"],
                path.parent,
                fpr=model_json["fpr"],
            )
            model.display_names = model_json["display_names"]
            bloom_path = model.base_path / model.slug() / "filter.bloom"
            model.bf = Bloom.load(
                str(bloom_path),
                hash_func=xxh3_64_intdigest,
            )
            return model

    def _generate_kmers(self, sequence: Seq, step: int = 1):
        """Generate kmers from the sequence"""
        num_kmers = ceil((len(sequence) - self.k + 1) / step)
        for i in range(num_kmers):
            start_pos = i * step
            kmer = sequence[start_pos : start_pos + self.k]
            minimizer = min(kmer, str(kmer.reverse_complement()))
            yield str(minimizer)
