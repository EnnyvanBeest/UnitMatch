# Sweep the UnitMatch batch-processing output trees for lingering
# ".processing*.lock" files (see batch_lock.py) and optionally remove them.
#
# A lock normally reclaims itself the next time some batch script tries to
# acquire that exact same lock_path and finds it older than
# batch_lock.STALE_AFTER_SECONDS -- see try_lock(). That doesn't help if
# nothing is currently trying to touch the group (e.g. no script has run
# since a crash), so this is a manual, one-off version of the same check:
# same 24h staleness threshold by default, but scans every lock under the
# given root(s) instead of waiting for a matching acquire attempt.
#
# Usage:
#   python clear_stale_locks.py                        # list-only, no changes
#   python clear_stale_locks.py --delete                # remove locks flagged STALE (>=24h old)
#   python clear_stale_locks.py --delete --min-age-hours 0   # remove every lock found, regardless of age
#   python clear_stale_locks.py --root <path> [--root <path> ...]  # sweep custom root(s)

import argparse
import fnmatch
import os
import time

import batch_lock

DEFAULT_ROOTS = [
    r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2",
    r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2_merged",
    r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026_V3_OnMergedData",
]

LOCK_GLOB = ".processing*.lock"


def find_locks(root):
    if not os.path.isdir(root):
        print(f"  (root not found, skipping: {root})")
        return
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if fnmatch.fnmatch(name, LOCK_GLOB):
                yield os.path.join(dirpath, name)


def describe(lock_path):
    """Return (age_seconds, one-line summary of the lock's contents) or None if it vanished."""
    try:
        age_s = time.time() - os.path.getmtime(lock_path)
    except OSError:
        return None
    try:
        with open(lock_path, "r") as f:
            content = f.read().strip().replace("\n", ", ")
    except OSError:
        content = "(unreadable)"
    return age_s, content


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--root", action="append", dest="roots", default=None,
        help="Root directory to sweep (repeatable). Defaults to the known BASE_OUTPUT trees.",
    )
    parser.add_argument(
        "--delete", action="store_true",
        help="Actually remove qualifying locks. Default is list-only (dry run).",
    )
    default_hours = batch_lock.STALE_AFTER_SECONDS / 3600
    parser.add_argument(
        "--min-age-hours", type=float, default=default_hours,
        help=f"Only delete locks at least this many hours old (default: {default_hours:.0f}, "
             "matching batch_lock's own stale threshold). Pass 0 with --delete to remove "
             "every lock found, including ones a run might still legitimately hold.",
    )
    args = parser.parse_args()
    roots = args.roots or DEFAULT_ROOTS

    found = removed = 0
    for root in roots:
        print(f"Scanning {root} ...")
        for lock_path in find_locks(root):
            found += 1
            info = describe(lock_path)
            if info is None:
                print(f"  ? {lock_path} (vanished while scanning)")
                continue
            age_s, content = info
            age_h = age_s / 3600
            qualifies = age_h >= args.min_age_hours
            tag = "STALE" if qualifies else "fresh"
            print(f"  [{tag}] {lock_path}  ({age_h:.1f}h old)  {content}")
            if args.delete and qualifies:
                try:
                    os.remove(lock_path)
                    removed += 1
                    print("    -> removed")
                except OSError as e:
                    print(f"    -> WARNING: could not remove: {e}")

    print(f"\n{found} lock(s) found.")
    if args.delete:
        print(f"{removed} removed.")
    else:
        print("Dry run -- nothing removed. Re-run with --delete to remove locks flagged STALE.")


if __name__ == "__main__":
    main()
