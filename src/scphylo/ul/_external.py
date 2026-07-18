"""Discovery and checked execution of optional external tools."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import IO


class ExternalToolError(RuntimeError):
    """Base error raised by an optional external-tool integration."""


class ExternalToolNotFoundError(ExternalToolError):
    """An executable or runtime artifact required by a tool was not found."""


class ExternalToolExecutionError(ExternalToolError):
    """An external tool failed or did not produce a usable result."""


def _directory_values(value):
    """Return normalized directories from a path-list setting."""
    if value is None:
        return []
    if isinstance(value, (str, os.PathLike)):
        values = os.fspath(value).split(os.pathsep)
    else:
        values = [os.fspath(item) for item in value]
    return [Path(item).expanduser() for item in values if item]


def _search_directories(tools_dir=None):
    """Build the external-tool search path in precedence order."""
    from scphylo import settings

    directories = []
    directories.extend(_directory_values(tools_dir))
    directories.extend(_directory_values(os.environ.get("SCPHYLO_TOOLS_DIR")))
    directories.extend(_directory_values(getattr(settings, "tools_dir", "")))
    directories.extend(_directory_values(os.environ.get("PATH")))

    unique = []
    seen = set()
    for directory in directories:
        key = os.path.abspath(directory)
        if key not in seen:
            seen.add(key)
            unique.append(directory)
    return unique


def _is_explicit_path(name):
    """Return whether a tool name includes a directory component."""
    path = Path(name).expanduser()
    return path.is_absolute() or len(path.parts) > 1


def _missing_message(name, appname, kind, path_option):
    """Build an actionable missing-tool message."""
    override = f" Pass `{path_option}=...`" if path_option else " Pass an explicit path"
    return (
        f"Cannot run {appname}: required {kind} `{name}` was not found."
        f"{override}, set `SCPHYLO_TOOLS_DIR`, set "
        "`scphylo.settings.tools_dir`, or add it to `PATH`."
    )


def resolve_executable(
    name,
    appname,
    *,
    tools_dir=None,
    path_option=None,
):
    """Resolve an executable and validate that it can be run.

    Resolution uses an explicit path first, then ``tools_dir``,
    ``SCPHYLO_TOOLS_DIR``, ``scphylo.settings.tools_dir``, and finally ``PATH``.
    """
    name = os.fspath(name)
    if _is_explicit_path(name):
        candidate = Path(name).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    else:
        search_path = os.pathsep.join(map(os.fspath, _search_directories(tools_dir)))
        candidate = shutil.which(name, path=search_path)
        if candidate is not None and Path(candidate).is_file():
            return str(Path(candidate).resolve())

    raise ExternalToolNotFoundError(
        _missing_message(name, appname, "executable", path_option)
    )


def resolve_external_file(
    name,
    appname,
    *,
    tools_dir=None,
    path_option=None,
):
    """Resolve a readable external-tool artifact such as a Java JAR."""
    name = os.fspath(name)
    if _is_explicit_path(name):
        candidates = [Path(name).expanduser()]
    else:
        candidates = [directory / name for directory in _search_directories(tools_dir)]

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.R_OK):
            return str(candidate.resolve())

    raise ExternalToolNotFoundError(
        _missing_message(name, appname, "file", path_option)
    )


def _log_tail(log_path, limit=4000):
    """Read a bounded tail from a tool log, if one is available."""
    if log_path is None:
        return ""
    try:
        text = Path(log_path).read_text(errors="replace")
    except OSError:
        return ""
    return text[-limit:].strip()


def run_external(
    args: Sequence[str | os.PathLike],
    appname: str,
    *,
    cwd: str | os.PathLike | None = None,
    stdout: int | IO | None = None,
    stderr: int | IO | None = subprocess.PIPE,
    timeout: float | None = None,
    log_path: str | os.PathLike | None = None,
):
    """Run an external tool without a shell and check its exit status."""
    command = [os.fspath(arg) for arg in args]
    if not command:
        raise ValueError("External command must contain at least one argument.")

    display_command = shlex.join(command)
    try:
        return subprocess.run(
            command,
            cwd=cwd,
            stdout=stdout,
            stderr=stderr,
            check=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as error:
        raise ExternalToolNotFoundError(
            f"Cannot run {appname}: executable `{command[0]}` was not found."
        ) from error
    except PermissionError as error:
        raise ExternalToolExecutionError(
            f"Cannot run {appname}: executable `{command[0]}` is not executable."
        ) from error
    except OSError as error:
        raise ExternalToolExecutionError(
            f"Cannot start {appname} with `{command[0]}`: {error}"
        ) from error
    except subprocess.TimeoutExpired as error:
        raise ExternalToolExecutionError(
            f"{appname} timed out after {error.timeout}s: {display_command}"
        ) from error
    except subprocess.CalledProcessError as error:
        detail = _log_tail(log_path)
        if not detail and isinstance(error.stderr, str):
            detail = error.stderr[-4000:].strip()
        suffix = f"\nLast tool output:\n{detail}" if detail else ""
        raise ExternalToolExecutionError(
            f"{appname} failed with exit code {error.returncode}: "
            f"{display_command}{suffix}"
        ) from error
