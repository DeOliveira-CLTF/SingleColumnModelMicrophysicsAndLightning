"""Eight lightweight entry points; scientific implementations live in the package."""
import argparse
import importlib
import sys

ACTIONS = {
    "download": {"era5": "cli:download-era5", "soundings": "cli:soundings", "cache": "cli:cache-inputs", "boundaries": "download.boundaries:main"},
    "prepare_forcing": {"assemble": "forcing.assemble:main", "build": "cli:forcing", "select-support": "forcing.select_support:main", "collocate-eda": "cli:collocate-eda"},
    "run_scm": {"controls": "cli:controls", "prepare-ensemble": "cli:prepare-ensemble", "ensemble": "cli:run-ensemble"},
    "process_lightning": {"observations": "cli:lightning", "aggregate": "observations.finalize:main", "diagnostics": "lightning.adaptations:main", "native": "evaluation.process_ensemble:main"},
    "calibrate": {"fit": "evaluation.process_ensemble:main", "adaptations": "lightning.adaptations:evaluate_main"},
    "evaluate": {"contracts": "evaluation.check_delivery:main", "stability": "evaluation.check_stability:main", "controls": "evaluation.check_controls:main", "vertical": "evaluation.vertical_audit:main", "migration": "evaluation.migration:main", "tables": "evaluation.finalize_execution:main"},
    "plot": {"all": "plotting.figures:main"},
    "reproduce": {"status": "cli:status", "audit": "cli:audit", "archive": "evaluation.archive_products:main", "utc": "evaluation.normalize_time:main", "inventory": "evaluation.prepare_products:main"},
}

def entry(kind):
    parser=argparse.ArgumentParser(description="Argentina 13–14 Dec 2018: "+kind)
    parser.add_argument("stage", choices=ACTIONS[kind])
    args, remainder=parser.parse_known_args()
    module, action=ACTIONS[kind][args.stage].split(":")
    sys.argv=[sys.argv[0]]+([action] if module=="cli" else [])+(["--evaluate"] if kind=="calibrate" and args.stage=="fit" else [])+remainder
    getattr(importlib.import_module("scm_lightning."+module), "main" if module=="cli" else action)()
