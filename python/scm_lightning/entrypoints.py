"""Eight lightweight entry points; scientific implementations live in the package."""
import argparse
import importlib
import sys

ACTIONS = {
    "download": {"merge": "download.merge:main", "era5": "cli:download-era5", "soundings": "cli:soundings", "cache": "cli:cache-inputs", "boundaries": "download.boundaries:main"},
    "prepare_forcing": {"event": "forcing.events:main", "assemble": "forcing.assemble:main", "build": "cli:forcing", "select-support": "forcing.select_support:main", "collocate-eda": "cli:collocate-eda"},
    "run_scm": {"independent": "evaluation.revision:main", "sensitivity": "ensemble.updraft_sensitivity:main", "controls": "cli:controls", "prepare-ensemble": "cli:prepare-ensemble", "ensemble": "cli:run-ensemble"},
    "process_lightning": {"profiles": "observations.profiles:main", "sensitivity": "evaluation.updraft_sensitivity:main", "observations": "cli:lightning", "aggregate": "observations.finalize:main", "diagnostics": "lightning.adaptations:main", "native": "evaluation.process_ensemble:main"},
    "calibrate": {"independent": "evaluation.revision:main", "fit": "evaluation.process_ensemble:main", "adaptations": "lightning.adaptations:evaluate_main"},
    "evaluate": {"merge": "evaluation.merge:main", "meteorology": "evaluation.meteorology:main", "independent": "evaluation.revision:main", "article": "evaluation.revision:evaluate", "contracts": "evaluation.revision:check", "stability": "evaluation.check_stability:main", "controls": "evaluation.check_controls:main", "vertical": "evaluation.vertical_audit:main", "migration": "evaluation.migration:main", "tables": "evaluation.finalize_execution:main"},
    "plot": {"article": "plotting.revision:main", "sensitivity": "evaluation.updraft_sensitivity:plot_main", "all": "plotting.figures:main"},
    "reproduce": {"independent": "evaluation.revision:main", "status": "cli:status", "audit": "cli:audit", "archive": "evaluation.archive_products:main", "utc": "evaluation.normalize_time:main", "inventory": "evaluation.prepare_products:main"},
}

def entry(kind):
    parser=argparse.ArgumentParser(description="Argentina independent calibration/evaluation: "+kind)
    parser.add_argument("stage", choices=ACTIONS[kind])
    parser.add_argument('--config')
    args, remainder=parser.parse_known_args()
    module, action=ACTIONS[kind][args.stage].split(":")
    config=['--config',args.config] if args.config else []
    sys.argv=[sys.argv[0]]+(config+[action] if module=="cli" else config)+(["--evaluate"] if kind=="calibrate" and args.stage=="fit" else [])+remainder
    getattr(importlib.import_module("scm_lightning."+module), "main" if module=="cli" else action)()
