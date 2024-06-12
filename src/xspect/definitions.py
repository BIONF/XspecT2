"""This module contains definitions for the XspecT package."""

from pathlib import Path
from os import getcwd


def get_xspect_root_path():
    """Return the root path for XspecT data."""
    root_path = Path(getcwd()) / "xspect-data"
    root_path.mkdir(exist_ok=True, parents=True)
    return root_path


def get_xspect_model_path():
    """Return the path to the XspecT models."""
    model_path = get_xspect_root_path() / "models"
    model_path.mkdir(exist_ok=True, parents=True)
    return model_path


def get_xspect_tmp_path():
    """Return the path to the XspecT temporary files."""
    tmp_path = get_xspect_root_path() / "tmp"
    tmp_path.mkdir(exist_ok=True, parents=True)
    return tmp_path


def get_xspect_upload_path():
    """Return the path to the XspecT upload directory."""
    upload_path = get_xspect_root_path() / "uploads"
    upload_path.mkdir(exist_ok=True, parents=True)
    return upload_path
