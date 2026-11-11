"""Command-line interface: ``barycon list | info ID | validate FILE... | bib``."""
from __future__ import annotations

import argparse
import sys

from .collection import load
from .dataset import read_entry


def main(argv=None):
    p = argparse.ArgumentParser(prog="barycon", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list", help="list entries")
    s.add_argument("--quantity")
    s.add_argument("--delta")
    s.add_argument("--pending", action="store_true", help="include pending entries")
    s = sub.add_parser("info", help="show one entry")
    s.add_argument("id")
    s = sub.add_parser("validate", help="validate data files against the schema")
    s.add_argument("files", nargs="+")
    sub.add_parser("bib", help="print BibTeX for all entries")
    a = p.parse_args(argv)

    if a.cmd == "list":
        col = load(include_pending=a.pending).filter(quantity=a.quantity, delta=a.delta)
        col.summary().pprint(max_lines=-1, max_width=-1)
    elif a.cmd == "info":
        print(load(include_pending=True)[a.id].info())
    elif a.cmd == "validate":
        bad = 0
        for f in a.files:
            try:
                e = read_entry(f, validate=True)
                print(f"OK    {f}  ({e.id})")
            except Exception as err:
                bad += 1
                print(f"FAIL  {f}: {err}")
        sys.exit(1 if bad else 0)
    elif a.cmd == "bib":
        print(load().bibtex())


if __name__ == "__main__":  # pragma: no cover
    main()
