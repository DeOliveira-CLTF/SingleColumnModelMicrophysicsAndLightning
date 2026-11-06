"""YAML/JSON configuration loading; manifests remain strict JSON."""
from pathlib import Path
import yaml


def load_config(path):
    result=yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    if not isinstance(result, dict):
        raise ValueError('Configuration must be a mapping')
    return result
