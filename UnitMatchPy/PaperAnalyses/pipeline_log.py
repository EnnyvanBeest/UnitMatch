# Per-location event log for the paper pipeline, so a location that drops out
# of a stage leaves a record of why (instead of only a console message).
#
# Every stage appends one JSON line per event to
#   cfg.LOG_DIR/<stage>_<machine>.jsonl
# (one file per machine, so machines running concurrently never write to the
# same file on the network share). Fields: time, machine, stage, group
# ("mouse/probe/location"), condition (e.g. "DeepUnitMatch", "UMPy",
# "xval_m3_1"), status ("done" | "failed" | "skipped"), message, traceback.
#
# Logging is best effort: a failure to write the log prints a warning and
# never interrupts the pipeline itself.
#
# Read back with read_events(); check_pipeline_completeness.py shows the
# latest logged event next to every location missing from a stage.

import datetime
import functools
import glob
import json
import os
import socket
import traceback

import pipeline_config as cfg

STATUSES = ("done", "failed", "skipped")


def log_event(stage, group, condition, status, message="", tb=None):
    """Append one event; never raises."""
    assert status in STATUSES, status
    record = {
        "time": datetime.datetime.now().isoformat(timespec="seconds"),
        "machine": socket.gethostname(),
        "stage": stage,
        "group": group,
        "condition": condition,
        "status": status,
        "message": str(message),
        "traceback": tb or "",
    }
    try:
        os.makedirs(cfg.LOG_DIR, exist_ok=True)
        path = os.path.join(cfg.LOG_DIR, f"{stage}_{socket.gethostname()}.jsonl")
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps(record) + "\n")
    except Exception as e:  # logging must never stop the pipeline
        print(f"  WARNING: could not write pipeline log ({e})")


def logged_run(stage, group_of, condition_of, output_of):
    """
    Decorator for one matching run: logs 'failed' (with traceback) when the
    run raises -- the exception is re-raised unchanged -- and otherwise
    'done' only if the run actually wrote its output file, 'failed' if it
    returned without it (several runs print an error and return early).

    group_of / condition_of / output_of map the call's (args, kwargs) to the
    logged group, condition, and the output file that marks success.
    """

    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            group = condition = "?"
            output = None
            try:
                group = group_of(args, kwargs)
                condition = condition_of(args, kwargs)
                output = output_of(args, kwargs)
            except Exception:
                pass
            try:
                result = fn(*args, **kwargs)
            except Exception as e:
                log_event(stage, group, condition, "failed", f"{type(e).__name__}: {e}", traceback.format_exc())
                raise
            if output is None or os.path.isfile(output):
                log_event(stage, group, condition, "done")
            else:
                log_event(stage, group, condition, "failed",
                          f"returned without writing {os.path.basename(output)} (see console output)")
            return result

        return wrapper

    return decorate


def read_events(stage=None):
    """All logged events (optionally of one stage) as a list of dicts, oldest first."""
    pattern = os.path.join(cfg.LOG_DIR, f"{stage or '*'}_*.jsonl")
    events = []
    for path in glob.glob(pattern):
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        events.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass
    return sorted(events, key=lambda e: e["time"])


def latest_events(stage=None):
    """{(stage, group, condition): latest event}."""
    return {(e["stage"], e["group"], e["condition"]): e for e in read_events(stage)}
