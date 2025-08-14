"""Probabilistic filter SVM model for sequence data"""

# pylint: disable=no-name-in-module, too-many-instance-attributes, arguments-renamed

import csv
import json
from pathlib import Path
from sklearn.svm import SVC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import cobs_index as cobs
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel
from xspect.definitions import fasta_endings, fastq_endings
from xspect.models.result import ModelResult


class ProbabilisticFilterSVMModel(ProbabilisticFilterModel):
    """
    Probabilistic filter SVM model for sequence data

    In addition to the standard probabilistic filter model, this model uses an SVM to predict
    labels based on their scores and training data. It requires the `scikit-learn` library
    to be installed.
    """

    def __init__(
        self,
        k: int,
        model_display_name: str,
        author: str | None,
        author_email: str | None,
        model_type: str,
        base_path: Path,
        kernel: str,
        c: float,
        fpr: float = 0.01,
        num_hashes: int = 7,
        training_accessions: dict[str, list[str]] | None = None,
        svm_accessions: dict[str, list[str]] | None = None,
    ) -> None:
        """
        Initialize the SVM model with the given parameters.

        In addition to the standard parameters, this model uses an SVM.
        Therefore, it requires the `kernel` and `C` parameters to be set.
        Furthermore, the `svm_accessions` parameter is used to store which accessions
        are used for training the SVM.

        Args:
            k (int): The k-mer size for the probabilistic filter.
            model_display_name (str): The display name of the model.
            author (str | None): The author of the model.
            author_email (str | None): The author's email address.
            model_type (str): The type of the model.
            base_path (Path): The base path where the model will be stored.
            kernel (str): The kernel type for the SVM (e.g., 'linear', 'rbf').
            c (float): Regularization parameter for the SVM.
            fpr (float, optional): False positive rate for the probabilistic filter.
            Defaults to 0.01.
            num_hashes (int, optional): Number of hashes for the probabilistic filter.
            Defaults to 7.
            training_accessions (dict[str, list[str]] | None, optional): Accessions used for
            training the probabilistic filter. Defaults to None.
            svm_accessions (dict[str, list[str]] | None, optional): Accessions used for
            training the SVM. Defaults to None.
        """
        super().__init__(
            k=k,
            model_display_name=model_display_name,
            author=author,
            author_email=author_email,
            model_type=model_type,
            base_path=base_path,
            fpr=fpr,
            num_hashes=num_hashes,
            training_accessions=training_accessions,
        )
        self.kernel = kernel
        self.c = c
        self.svm_accessions = svm_accessions

    def to_dict(self) -> dict:
        """
        Convert the model to a dictionary representation

        Returns:
            dict: A dictionary containing the model's parameters and state.
        """
        return super().to_dict() | {
            "kernel": self.kernel,
            "C": self.c,
            "svm_accessions": self.svm_accessions,
        }

    def set_svm_params(self, kernel: str, c: float) -> None:
        """
        Set the parameters for the SVM

        Args:
            kernel (str): The kernel type for the SVM (e.g., 'linear', 'rbf').
            c (float): Regularization parameter for the SVM.
        """
        self.kernel = kernel
        self.c = c
        self.save()

    def fit(
        self,
        dir_path: Path,
        svm_path: Path,
        display_names: dict[str, str] | None = None,
        svm_step: int = 1,
        training_accessions: dict[str, list[str]] | None = None,
        svm_accessions: dict[str, list[str]] | None = None,
    ) -> None:
        """
        Fit the SVM to the sequences and labels.

        This method first trains the probabilistic filter model and then calculates scores for
        the SVM training. It expects the sequences to be in the specified directory and the SVM
        training sequences to be in the specified SVM path. The scores are saved in a CSV file
        for later use.

        Args:
            dir_path (Path): The directory containing the training sequences.
            svm_path (Path): The directory containing the SVM training sequences.
            display_names (dict[str, str] | None): A mapping of accession IDs to display names.
            svm_step (int): Step size for sparse sampling in SVM training.
            training_accessions (dict[str, list[str]] | None): Accessions used for training the
            probabilistic filter.
            svm_accessions (dict[str, list[str]] | None): Accessions used for training the SVM.
        """

        # Since the SVM works with score data, we need to train
        # the underlying data structure for score generation first
        super().fit(
            dir_path,
            display_names=display_names,
            training_accessions=training_accessions,
        )

        self.svm_accessions = svm_accessions

        # calculate scores for SVM training
        score_list = []

        for species_folder in svm_path.iterdir():
            if not species_folder.is_dir():
                continue
            for file in species_folder.iterdir():
                if file.suffix[1:] not in fasta_endings + fastq_endings:
                    continue
                print(f"Calculating {file.name} scores for SVM training...")
                res = super().predict(file, step=svm_step)
                scores = res.get_scores()["total"]
                accession = file.stem
                label_id = species_folder.name

                # format scores for csv
                scores = dict(sorted(scores.items()))
                scores = ",".join([str(score) for score in scores.values()])
                scores = f"{accession},{scores},{label_id}"
                score_list.append(scores)

        # csv header
        keys = list(self.display_names.keys())
        keys.sort()
        score_list.insert(0, f"file,{','.join(keys)},label_id")

        with open(
            self.base_path / self.slug() / "scores.csv", "w", encoding="utf-8"
        ) as file:
            file.write("\n".join(score_list))

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
        Predict the labels of the sequences.

        This method uses the SVM to predict labels based on the scores generated
        from the sequences. It expects the sequences to be in a format compatible
        with the probabilistic filter model, and it will return a `ModelResult`.

        Args:
            sequence_input (SeqRecord | list[SeqRecord] | SeqIO.FastaIO.FastaIterator |
            SeqIO.QualityIO.FastqPhredIterator | Path): The input sequences to predict.
            filter_ids (list[str], optional): A list of IDs to filter the predictions.
            step (int, optional): Step size for sparse sampling. Defaults to 1.
            display_name (bool): Includes a display name for each tax_ID.

        Returns:
            ModelResult: The result of the prediction containing hits, number of kmers, and the
            predicted label.
        """
        # get scores and format them for the SVM
        res = super().predict(sequence_input, filter_ids, step, display_name)
        svm_scores = dict(sorted(res.get_scores()["total"].items()))
        svm_scores = [list(svm_scores.values())]

        svm = self._get_svm(filter_ids)
        return ModelResult(
            self.slug(),
            res.hits,
            res.num_kmers,
            sparse_sampling_step=step,
            prediction=str(svm.predict(svm_scores)[0]),
        )

    def _get_svm(self, id_keys) -> SVC:
        """
        Get the SVM for the given id keys.

        This method loads the SVM model from the scores CSV file and trains it
        using the scores from the CSV. If `id_keys` is provided, it filters the
        training data to only include those keys.

        Args:
            id_keys (list[str] | None): A list of IDs to filter the training data.
                If None, all data is used.

        Returns:
            SVC: The trained SVM model.
        """
        svm = SVC(kernel=self.kernel, C=self.c)
        # parse csv
        with open(
            self.base_path / self.slug() / "scores.csv", "r", encoding="utf-8"
        ) as file:
            file.readline()
            x_train = []
            y_train = []
            for row in csv.reader(file):
                if id_keys is None or row[-1] in id_keys:
                    x_train.append(row[1:-1])
                    y_train.append(row[-1])

        # train svm
        svm.fit(x_train, y_train)
        return svm

    @staticmethod
    def load(path: Path) -> "ProbabilisticFilterSVMModel":
        """
        Load the model from disk

        Loads the model from the specified path. The path should point to a JSON file
        containing the model's parameters and state. It also checks for the existence of
        the COBS index file.

        Args:
            path (Path): The path to the model JSON file.

        Returns:
            ProbabilisticFilterSVMModel: The loaded model instance.
        """
        with open(path, "r", encoding="utf-8") as file:
            json_object = file.read()
            model_json = json.loads(json_object)
            model = ProbabilisticFilterSVMModel(
                model_json["k"],
                model_json["model_display_name"],
                model_json["author"],
                model_json["author_email"],
                model_json["model_type"],
                path.parent,
                model_json["kernel"],
                model_json["C"],
                fpr=model_json["fpr"],
                num_hashes=model_json["num_hashes"],
                training_accessions=model_json["training_accessions"],
                svm_accessions=model_json["svm_accessions"],
            )
            model.display_names = model_json["display_names"]

            p = model.get_cobs_index_path()
            if not Path(p).exists():
                raise FileNotFoundError(f"Index file not found at {p}")
            model.index = cobs.Search(p, True)

            return model
