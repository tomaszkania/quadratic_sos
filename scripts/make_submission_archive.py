"""Build the TOMS software-component submission archive."""

from __future__ import annotations

from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_INIT = ROOT / "src" / "quadratic_diagonal" / "__init__.py"

EXCLUDED_DIR_NAMES = {
    ".git",
    ".pytest_cache",
    "__pycache__",
    "build",
    "dist",
    "submission",
}
EXCLUDED_SUFFIXES = {
    ".aux",
    ".bbl",
    ".blg",
    ".fdb_latexmk",
    ".fls",
    ".log",
    ".out",
    ".pyc",
    ".pyo",
    ".synctex.gz",
    ".toc",
}


def read_version() -> str:
    """Read the package version.

    Returns
    -------
    str
        Version string exported by ``quadratic_diagonal.__version__``.

    Raises
    ------
    RuntimeError
        If the version string cannot be found.
    """
    text = PACKAGE_INIT.read_text(encoding="utf-8")
    match = re.search(r'__version__\s*=\s*"([^"]+)"', text)
    if match is None:
        raise RuntimeError("Could not read package version from __init__.py")
    return match.group(1)


def is_excluded(path: Path) -> bool:
    """Return whether ``path`` should be omitted from the archive.

    Parameters
    ----------
    path : pathlib.Path
        Absolute path under the repository root.

    Returns
    -------
    bool
        ``True`` if the path is a cache, build artefact, or generated
        submission archive that should not be included.
    """
    relative = path.relative_to(ROOT)
    parts = set(relative.parts)
    if parts & EXCLUDED_DIR_NAMES:
        return True
    if any(part.endswith(".egg-info") for part in relative.parts):
        return True
    name = path.name
    if any(name.endswith(suffix) for suffix in EXCLUDED_SUFFIXES):
        return True
    return False


def iter_archive_files() -> list[Path]:
    """Return repository files to include in the submission archive.

    Returns
    -------
    list[pathlib.Path]
        Sorted absolute paths to regular files under the repository root.
    """
    files: list[Path] = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and not is_excluded(path):
            files.append(path)
    return files


def build_archive() -> Path:
    """Create the software-component archive.

    Returns
    -------
    pathlib.Path
        Path to the generated zip file.
    """
    version = read_version()
    out_dir = ROOT / "submission"
    out_dir.mkdir(exist_ok=True)
    archive_path = out_dir / f"quadratic_diagonal-{version}-toms-software.zip"
    prefix = f"quadratic_diagonal-{version}"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in iter_archive_files():
            relative = path.relative_to(ROOT)
            archive.write(path, Path(prefix) / relative)
    return archive_path


def main() -> None:
    """Build the archive and print its path."""
    archive_path = build_archive()
    print(archive_path)


if __name__ == "__main__":
    main()
