# Lightweight cross-machine file lock for batch scripts that read/write a
# shared network folder (e.g. \\znas.cortexlab.net\...) from several
# computers at once.
#
# Locking is advisory and file-based: "acquiring" means creating a lock file
# with os.O_CREAT | os.O_EXCL, which maps to CreateFile(..., CREATE_NEW) on
# Windows and is *supposed* to be atomic both locally and on SMB network
# shares -- in practice this has been observed to fail on at least one share
# (two processes both got a successful open() for the same path within a
# short window of each other), so try_lock() also does a write-then-verify
# check after "winning" the open() -- see its docstring.
# There's no cross-machine way to tell whether the process that created a
# lock is still alive, so a lock is instead reclaimed once it's older than
# `stale_after` seconds -- long enough that it's very unlikely a real run is
# still going, short enough that a crashed machine doesn't block the others
# indefinitely.

import contextlib
import datetime
import os
import random
import socket
import time
import uuid

STALE_AFTER_SECONDS = 24 * 3600  # reclaim locks older than this
EARLY_RECLAIM_GRACE_SECONDS = 2 * 3600  # see try_lock()'s redo_from_date param
LOCK_VERIFY_DELAY_SECONDS = 1.5  # see try_lock()'s write-then-verify step


def sentinel_is_fresh(sentinel_path, redo_from_date=None):
    """
    Return True when `sentinel_path` exists and doesn't need redoing.

    This is the shared "should we skip this unit of work" check used in place
    of a plain os.path.isfile() so batch scripts can recover from an algorithm
    change without a blanket REDO=True (which forces every group through the
    lock/reprocess path, including ones a concurrent run already redid since
    the change -- see run_deepunitmatch_batch_onMerged.py's module docstring
    for the incident this was added for).

    redo_from_date : datetime.datetime or None
        If None, only existence is checked (classic "skip if present"
        behaviour). If set, a sentinel older than this datetime is treated as
        stale (needs redoing) even though it exists -- e.g. set it to the date
        a bug fix landed so only output computed before the fix gets redone.
        This is independent of, and composes with, try_lock()/STALE_AFTER_SECONDS
        above: staleness decides *whether* a group needs work; the lock still
        decides whether *this* process is allowed to do that work right now.
    """
    if not os.path.isfile(sentinel_path):
        return False
    if redo_from_date is None:
        return True
    mtime = datetime.datetime.fromtimestamp(os.path.getmtime(sentinel_path))
    return mtime >= redo_from_date


@contextlib.contextmanager
def try_lock(lock_path, stale_after=STALE_AFTER_SECONDS, redo_from_date=None,
             early_reclaim_grace=EARLY_RECLAIM_GRACE_SECONDS):
    """
    Attempt to acquire an exclusive lock at `lock_path`.

    Yields True if the lock was acquired -- the caller should do the work;
    the lock file is removed automatically on exit, whether the work
    succeeds or raises. Yields False if another run already holds the lock
    (caller should skip this unit of work).

    os.O_CREAT | os.O_EXCL is supposed to make lock creation atomic even on
    SMB network shares (see module docstring), but that's been observed to
    not hold in practice on at least one share: two processes both got a
    successful open() for the same lock_path within a short window of each
    other and both proceeded to do the work. To guard against that, after
    "winning" the open() we write a unique token, wait briefly for any
    near-simultaneous writer to land its own write, then re-read the file
    back -- only the process whose token is still there actually holds the
    lock; the loser yields False without touching the file (it now belongs
    to whoever really won).

    redo_from_date : datetime.datetime or None
        If set (typically the same value passed to sentinel_is_fresh()), a
        lock older than `early_reclaim_grace` whose mtime predates this date
        is reclaimed immediately instead of waiting the full `stale_after`
        window -- it was necessarily taken by a process running code from
        before a known fix, so its output needs redoing regardless of
        whether that process is still running. `early_reclaim_grace` still
        guards against yanking the lock out from under a run that started
        shortly before the fix landed and is legitimately still in
        progress -- there's no cross-machine way to check liveness directly
        (see module docstring), so this is a smaller, but still real, safety
        margin than the full stale_after window.
    """
    os.makedirs(os.path.dirname(lock_path), exist_ok=True)

    try:
        lock_mtime_epoch = os.path.getmtime(lock_path)
        age = time.time() - lock_mtime_epoch
        reclaim, reason = False, ""
        if age > stale_after:
            reclaim, reason = True, f"{age / 3600:.1f}h old"
        elif redo_from_date is not None and age > early_reclaim_grace:
            lock_mtime = datetime.datetime.fromtimestamp(lock_mtime_epoch)
            if lock_mtime < redo_from_date:
                reclaim, reason = True, f"{age / 3600:.1f}h old, predates redo_from_date={redo_from_date}"
        if reclaim:
            print(f"  Reclaiming lock ({reason}): {lock_path}")
            os.remove(lock_path)
    except FileNotFoundError:
        pass
    except OSError as e:
        print(f"  WARNING: could not inspect/remove lock {lock_path}: {e}")

    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        yield False
        return
    except OSError as e:
        print(f"  WARNING: could not create lock {lock_path}: {e}")
        yield False
        return

    # Unique per-attempt token to identify our own write below, distinct from
    # whatever a concurrent racer writes to the same path.
    token = f"{socket.gethostname()}:{os.getpid()}:{time.time_ns()}:{uuid.uuid4().hex}"

    try:
        with os.fdopen(fd, "w") as f:
            f.write(
                f"host={socket.gethostname()}\n"
                f"pid={os.getpid()}\n"
                f"started={time.ctime()}\n"
                f"token={token}\n"
            )
    except OSError as e:
        print(f"  WARNING: could not write lock {lock_path}: {e}")
        yield False
        return

    # Give any near-simultaneous writer time to land its own write, then
    # re-read: whichever process wrote last "wins" the file's contents on
    # disk, so only proceed if it's still us.
    time.sleep(LOCK_VERIFY_DELAY_SECONDS + random.uniform(0, LOCK_VERIFY_DELAY_SECONDS))
    try:
        with open(lock_path, "r") as f:
            on_disk = f.read()
    except OSError as e:
        print(f"  WARNING: could not verify lock {lock_path}: {e}")
        yield False
        return

    if token not in on_disk:
        print(f"  Lost lock race for {lock_path} (another process's write won); not proceeding.")
        yield False
        return

    try:
        yield True
    finally:
        # Only remove the lock if it's still ours -- e.g. a stale-reclaim by
        # another process while we were running would mean it's no longer
        # our lock to delete.
        try:
            with open(lock_path, "r") as f:
                on_disk = f.read()
            if token in on_disk:
                os.remove(lock_path)
            else:
                print(f"  WARNING: lock {lock_path} was reclaimed by another process before we finished; not removing it.")
        except FileNotFoundError:
            pass
        except OSError as e:
            print(f"  WARNING: could not remove lock {lock_path}: {e}")
