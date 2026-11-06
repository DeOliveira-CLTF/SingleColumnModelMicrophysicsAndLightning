"""Public project entry point; run with --help for supported stages."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from scm_lightning.entrypoints import entry
if __name__ == "__main__":
    entry("reproduce")
