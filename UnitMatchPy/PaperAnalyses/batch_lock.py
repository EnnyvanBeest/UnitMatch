# Lightweight cross-machine file lock for batch scripts that read/write a
# shared network folder (e.g. \\znas.cortexlab.net\...) from several
# computers at once.
#
# Locking is advisory and file-based: "acquiring" means atomically creating
# a lock file with os.O_CREAT | os.O_EXCL, which is atomic both locally and
# on SMB network shares (maps to CreateFile(..., CREATE_NEW) on Windows).
# There's no cross-machine way to tell whether the process that created a
# lock is still alive, so a lock is instead reclaimed once it's older than
# `stale_after` seconds -- long enough that it's very unlikely a real run is
# still going, short enough that a crashed machine doesn't block the others
# indefinitely.

import contextlib
import os
import socket
import time

STALE_AFTER_SECONDS = 24 * 3600  # reclaim locks older than this


@contextlib.contextmanager
def try_lock(lock_path, stale_after=STALE_AFTER_SECONDS):
    """
    Attempt to acquire an exclusive lock at `lock_path`.

    Yields True if the lock was acquired -- the caller should do the work;
    the lock file is removed automatically on exit, whether the work
    succeeds or raises. Yields False if another run already holds the lock
    (caller should skip this unit of work).
    """
    os.makedirs(os.path.dirname(lock_path), exist_ok=True)

    try:
        age = time.time() - os.path.getmtime(lock_path)
        if age > stale_after:
            print(f"  Reclaiming stale lock ({age / 3600:.1f}h old): {lock_path}")
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

    try:
        with os.fdopen(fd, "w") as f:
            f.write(
                f"host={socket.gethostname()}\n"
                f"pid={os.getpid()}\n"
                f"started={time.ctime()}\n"
            )
        yield True
    finally:
        try:
            os.remove(lock_path)
        except OSError as e:
            print(f"  WARNING: could not remove lock {lock_path}: {e}")
