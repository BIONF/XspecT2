"""FastAPI application for XspecT."""

from shutil import copyfileobj
from fastapi import FastAPI, UploadFile, BackgroundTasks
from xspect.definitions import get_xspect_upload_path
from xspect.download_filters import download_test_filters
from xspect.file_io import get_records_by_id
import xspect.model_management as mm
from xspect.train import train_ncbi

app = FastAPI()


@app.get("/download-filters")
def download_filters():
    """Download filters."""
    download_test_filters("https://xspect2.s3.eu-central-1.amazonaws.com/models.zip")


@app.get("/classify")
def classify(genus: str, file: str, meta: bool = False, step: int = 500):
    """Classify uploaded sample."""

    sequence_input = get_xspect_upload_path() / file
    species_filter_model = mm.get_species_model(genus)

    if meta:
        genus_filter_model = mm.get_genus_model(genus)
        filter_res = genus_filter_model.predict(sequence_input, step=step)
        filtered_sequences = filter_res.get_filtered_subsequences(genus, 0.7)
        sequence_input = get_records_by_id(sequence_input, filtered_sequences)
    res = species_filter_model.predict(sequence_input, step=step)

    return {"prediction": res.prediction, "scores": res.get_scores()["total"]}


@app.post("/train")
def train(genus: str, background_tasks: BackgroundTasks, svm_steps: int = 1):
    """Train NCBI model."""
    background_tasks.add_task(train_ncbi, genus, svm_steps)

    return {"message": "Training started."}


@app.get("/list-models")
def list_models():
    """List available models."""
    return mm.get_models()


@app.get("/model-metadata")
def get_model_metadata(model_slug: str):
    """Get metadata of a model."""
    return mm.get_model_metadata(model_slug)


@app.post("/model-metadata")
def post_model_metadata(model_slug: str, author: str, author_email: str):
    """Update metadata of a model."""
    try:
        mm.update_model_metadata(model_slug, author, author_email)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Metadata updated."}


@app.post("/model-display-name")
def post_model_display_name(model_slug: str, filter_id: str, display_name: str):
    """Update display name of a filter in a model."""
    try:
        mm.update_model_display_name(model_slug, filter_id, display_name)
    except ValueError as e:
        return {"error": str(e)}
    return {"message": "Display name updated."}


@app.post("/upload-file")
def upload_file(file: UploadFile):
    """Upload file to the server."""
    upload_path = get_xspect_upload_path() / file.filename

    if not upload_path.exists():
        try:
            with upload_path.open("wb") as buffer:
                copyfileobj(file.file, buffer)
        finally:
            file.file.close()

    return {"filename": file.filename}
