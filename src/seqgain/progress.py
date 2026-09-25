"""Throttled terminal heartbeat; one reporting thread, no worker processes."""
import sys
import shutil
import threading
import time


def duration(seconds):
    seconds = int(max(0, seconds))
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


class Progress:
    """Phase-local counts and approximate ETA when a trustworthy total is known.

    Disabled instances are also the library API's silent default. A context
    manager owns the single heartbeat thread, including exception cleanup.
    """

    def __init__(self, enabled=False, interval=2.0, stream=None):
        self.enabled = enabled
        self.interval = interval
        self.stream = stream if stream is not None else sys.stderr
        self.interactive = bool(getattr(self.stream, "isatty", lambda: False)())
        self.label = None
        self.total = None
        self.count = 0
        self.unit = "items"
        self.started = time.monotonic()
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._thread = None

    def __enter__(self):
        if self.enabled and self.interactive:
            self._thread = threading.Thread(target=self._heartbeat, name="seqgain-progress", daemon=True)
            self._thread.start()
        return self

    def __exit__(self, exc_type, exc, traceback):
        self._stop.set()
        if self._thread:
            self._thread.join()
        self.finish("stopped" if exc_type else "done")

    def finish(self, status="done"):
        """End the active line before warnings, the result path, or shutdown."""
        with self._lock:
            self._emit(status)
            self.label = None

    def phase(self, label, total=None, unit="items"):
        if not self.enabled:
            return
        with self._lock:
            self._emit("done")
            self.label, self.total, self.unit = label, total, unit
            self.count, self.started = 0, time.monotonic()
            self._emit()

    def advance(self, amount=1):
        if self.enabled:
            with self._lock:
                self.count += amount

    def update(self, count):
        if self.enabled:
            with self._lock:
                self.count = count

    def _heartbeat(self):
        while not self._stop.wait(self.interval):
            with self._lock:
                self._emit()

    def _emit(self, status=None):
        if not self.enabled or self.label is None:
            return
        if not self.interactive and status is None:
            return
        elapsed = time.monotonic() - self.started
        details = [f"elapsed {duration(elapsed)}"]
        if self.total is not None:
            percent = 100 * self.count / self.total if self.total else 100
            details.insert(0, f"{self.count:,}/{self.total:,} {self.unit} ({percent:.1f}%)")
        elif self.count:
            details.insert(0, f"{self.count:,} {self.unit}")
        if self.count and elapsed > 0:
            details.append(f"{self.count / elapsed:,.1f} {self.unit}/s")
        if status:
            details.append(status)
        elif self.total is not None and self.count > 0:
            details.append(f"ETA ~{duration(elapsed * max(0, self.total-self.count) / self.count)}")
        else:
            details.append("running; ETA unavailable")
        label = self.label.replace("\n", " ").replace("\r", " ")
        line = f"[{label}] " + " · ".join(details)
        if self.interactive:
            # Reserve one column to avoid terminal auto-wrap; shorten long
            # contig labels first so the percentage remains visible.
            width = max(1, shutil.get_terminal_size().columns - 1)
            if len(line) > width:
                limit = max(1, width // 3)
                label = label[:limit-1] + "…" if len(label) > limit else label
                line = f"[{label}] " + " · ".join(details)
                line = line[:width-1] + "…" if len(line) > width else line
            print("\r\033[2K" + line, end="\n" if status else "", file=self.stream, flush=True)
        else:
            print(line, file=self.stream, flush=True)
