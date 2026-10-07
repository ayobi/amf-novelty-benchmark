#!/usr/bin/env python3
"""Find UNITE FASTA files on local drives and say which one crosswalk.py should use.

Walks the given roots (default: home, /media, /mnt, /data) and inspects every FASTA whose name
or path mentions UNITE or follows UNITE's release naming (sh_general_release*, sh_refs_qiime*),
plain, gzipped, or inside .tgz/.tar.gz/.zip archives (read in place, nothing is extracted).
For each file it counts records, records carrying an SH code, records with a taxonomy string in
the header, and Glomeromycota records, and reads the UNITE version from the SH codes (.10FU = 10).

crosswalk.py needs SH codes and the taxonomy in the header, which the general FASTA release has
and the QIIME release does not (QIIME keeps taxonomy in a separate file). Ranking among usable
files: not a developer ('_dev') set, highest UNITE version, dynamic thresholds, most Glomeromycota,
then a plain file over an archive member.

Usage:
  python3 find_unite.py                      # default roots
  python3 find_unite.py ~/Downloads /media   # specific roots
"""
import argparse
import gzip
import io
import os
import re
import sys
import tarfile
import time
import zipfile

FASTA_EXT = (".fasta", ".fa", ".fas", ".fna", ".fasta.gz", ".fa.gz", ".fas.gz", ".fna.gz")
ARCHIVE_EXT = (".tgz", ".tar.gz", ".tar", ".zip")
NAME_HINT = re.compile(r"unite|sh_general_release|sh_refs|sh_qiime|utax_reference", re.I)
SKIP_DIRS = {"proc", "sys", "dev", ".git", "node_modules", "__pycache__", ".Trash", ".cache"}
SH = re.compile(r"SH\d+\.(\d+)FU")


def inspect(fh):
    """Stream a text FASTA; return counts and the first header."""
    n = sh = tax = glom = 0
    versions, first = set(), ""
    for line in fh:
        if not line.startswith(">"):
            continue
        n += 1
        if not first:
            first = line[1:].strip()
        m = SH.search(line)
        if m:
            sh += 1
            versions.add(int(m.group(1)))
        if "k__" in line or "p__" in line:
            tax += 1
        if "Glomeromycota" in line:
            glom += 1
    return {"records": n, "sh": sh, "tax": tax, "glom": glom, "versions": sorted(versions), "first": first}


def candidates(roots):
    seen = set()
    for root in roots:
        root = os.path.expanduser(root)
        if not os.path.isdir(root):
            continue
        for d, dirs, files in os.walk(root, onerror=lambda e: None):
            dirs[:] = [x for x in dirs if x not in SKIP_DIRS and not x.startswith(".Trash")]
            hinted_dir = bool(NAME_HINT.search(d))
            for f in files:
                low = f.lower()
                if not (low.endswith(FASTA_EXT) or low.endswith(ARCHIVE_EXT)):
                    continue
                if not (hinted_dir or NAME_HINT.search(f)):
                    continue
                p = os.path.realpath(os.path.join(d, f))
                if p not in seen:
                    seen.add(p)
                    yield p


def scan(path):
    """Yield (display path, extract hint, stats) for a FASTA or for each FASTA member of an archive."""
    low = path.lower()
    try:
        if low.endswith(FASTA_EXT):
            opener = gzip.open if low.endswith(".gz") else open
            with opener(path, "rt", errors="replace") as fh:
                yield path, "", inspect(fh)
        elif low.endswith(".zip"):
            with zipfile.ZipFile(path) as z:
                for m in z.namelist():
                    if m.lower().endswith(FASTA_EXT) and not m.endswith("/"):
                        with z.open(m) as raw:
                            fh = io.TextIOWrapper(gzip.GzipFile(fileobj=raw) if m.lower().endswith(".gz") else raw,
                                                  errors="replace")
                            yield f"{path}::{m}", (f"unzip -j '{path}' '{m}' -d ../ref/unite", os.path.basename(m)), inspect(fh)
        else:
            with tarfile.open(path, "r:*") as t:
                for m in t:
                    if m.isfile() and m.name.lower().endswith(FASTA_EXT):
                        raw = t.extractfile(m)
                        fh = io.TextIOWrapper(gzip.GzipFile(fileobj=raw) if m.name.lower().endswith(".gz") else raw,
                                              errors="replace")
                        yield f"{path}::{m.name}", (f"tar -xf '{path}' -C ../ref/unite '{m.name}'", m.name), inspect(fh)
    except (OSError, EOFError, tarfile.TarError, zipfile.BadZipFile) as e:
        print(f"  (could not read {path}: {e})", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="*", default=["~", "/media", "/mnt", "/data"])
    a = ap.parse_args()
    t0 = time.time()
    rows = []
    for p in candidates(a.roots):
        for shown, hint, st in scan(p):
            if st["records"]:
                rows.append((shown, hint, st))
    print(f"searched {', '.join(a.roots)} in {time.time() - t0:.0f} s: {len(rows)} FASTA file(s) look like UNITE\n")
    if not rows:
        print("nothing found; try other roots, or: locate -i unite | grep -iE 'fasta|tgz|zip'")
        return
    usable = []
    for shown, hint, st in rows:
        ok = st["sh"] > 0 and st["tax"] > 0 and st["glom"] > 0
        name = os.path.basename(shown).lower()
        why = [] if ok else [w for w, bad in (("no SH codes", not st["sh"]), ("no taxonomy in headers", not st["tax"]),
                                              ("no Glomeromycota", not st["glom"])) if bad]
        print(f"{'USABLE ' if ok else 'skip   '} {shown}")
        print(f"         {st['records']:,} records, {st['sh']:,} with SH codes, {st['glom']:,} Glomeromycota; "
              f"UNITE version {','.join(map(str, st['versions'])) or '?'}"
              + (f"; {', '.join(why)}" if why else ""))
        print(f"         first header: {st['first'][:110]}")
        if ok:
            score = (0 if "_dev" in name or "developer" in shown.lower() else 1,
                     max(st["versions"]) if st["versions"] else 0, 1 if "dynamic" in name else 0, st["glom"],
                     0 if hint else 1)  # at a tie, a plain file beats an archive member
            usable.append((score, shown, hint, st))
    if not usable:
        print("\nno usable file: crosswalk.py needs the UNITE general FASTA release (SH codes and taxonomy in headers)")
        return
    usable.sort(key=lambda x: x[0], reverse=True)
    _, shown, hint, st = usable[0]
    print(f"\nsuggested: {shown}")
    if hint:
        cmd, dest = hint
        print(f"  extract:  mkdir -p ../ref/unite && {cmd}")
        print(f"  then:     --unite ../ref/unite/{dest}")
    else:
        print(f"  use:      --unite '{shown}'   (or copy it into ../ref/unite/ first)")


if __name__ == "__main__":
    sys.exit(main())
