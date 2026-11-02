#!/usr/bin/env python3
"""audit_run.py - audited execution of a hashed LEAD-chi script (registration Annex A.7).

  audit_run.py [--aggregate PATH ...] [--out PATH ...] -- <interpreter> <script> [args ...]

Appends one JSON line per run to $LX/logs/audit.jsonl (hash-chained: each line carries the SHA-256 of the previous
line): UTC start/end, host, user, the command, the SHA-256 of the script and of every existing file named in the
arguments (inputs, before the run), the SHA-256 of every declared output (after the run), the git commit of $LX
(if it is a repository) and the exit status.

Technical enforcement of "no per-code values in visible outputs": every --aggregate output must be JSON and must
(1) contain no coded subject ID (read from blind/processing_order.txt) anywhere, as key or value, and
(2) contain no numeric array longer than 12 elements (a per-session vector cannot hide there).
An aggregate output that fails is renamed <name>.REJECTED, the run is logged as "rejected" and audit_run.py exits
with status 3. --out outputs (for example hidden code-keyed tables) are only hashed.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

MAX_ARRAY = 12


def coded_ids(lx):
    p = os.path.join(lx, "blind", "processing_order.txt")
    if not os.path.exists(p):
        return set()
    return {ln.split("\t")[0] for ln in open(p).read().split("\n") if ln.strip()}


def check_aggregate(path, ids):
    try:
        obj = json.load(open(path))
    except (ValueError, OSError) as e:
        return [f"not valid JSON: {e}"]
    bad = []
    txt = json.dumps(obj)
    for c in ids:
        if re.search(rf"(?<![0-9a-f]){re.escape(c)}(?![0-9a-f])", txt):
            bad.append("contains a coded subject ID")
            break

    def walk(o, where="$"):
        if isinstance(o, dict):
            for k, v in o.items():
                walk(v, f"{where}.{k}")
        elif isinstance(o, list):
            if sum(isinstance(v, (int, float)) and not isinstance(v, bool) for v in o) > MAX_ARRAY:
                bad.append(f"numeric array longer than {MAX_ARRAY} at {where}")
            for i, v in enumerate(o):
                walk(v, f"{where}[{i}]")
    walk(obj)
    return bad


def git_head(lx):
    try:
        return subprocess.run(["git", "-C", lx, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
    except OSError:
        return None


def main():
    if "--" not in sys.argv:
        sys.exit("usage: audit_run.py [--aggregate PATH ...] [--out PATH ...] -- <interpreter> <script> [args ...]")
    i = sys.argv.index("--")
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--aggregate", action="append", default=[])
    a.add_argument("--out", action="append", default=[])
    args = a.parse_args(sys.argv[1:i])
    cmd = sys.argv[i + 1:]
    if len(cmd) < 2:
        sys.exit("give the interpreter and the script")
    lx = L.lx_root()
    for p in args.aggregate + args.out:
        L.assert_inside(p)
    rec = {"utc_start": L.utcnow(), "host": L.socket.gethostname(), "user": L.getpass.getuser(), "cmd": cmd,
           "script_sha256": L.sha256_file(cmd[1]) if os.path.isfile(cmd[1]) else None, "git_head": git_head(lx),
           "inputs": {x: L.sha256_file(x) for x in cmd[2:] if os.path.isfile(x)}}
    res = subprocess.run(cmd)
    rec["exit_status"] = res.returncode
    rec["outputs"] = {p: (L.sha256_file(p) if os.path.isfile(p) else None) for p in args.aggregate + args.out}
    ids = coded_ids(lx)
    problems = {}
    for p in args.aggregate:
        if os.path.isfile(p):
            b = check_aggregate(p, ids)
            if b:
                problems[p] = b
                os.replace(p, p + ".REJECTED")
    rec["status"] = "rejected" if problems else ("ok" if res.returncode == 0 else "failed")
    rec["rejections"] = problems
    rec["utc_end"] = L.utcnow()
    log = os.path.join(lx, "logs", "audit.jsonl")
    L.ensure_parent(log)
    prev = None
    if os.path.exists(log):
        lines = open(log, "rb").read().splitlines()
        prev = L.sha256_bytes(lines[-1]) if lines else None
    rec["prev_line_sha256"] = prev
    with open(log, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")
    if problems:
        print(f"AUDIT REJECTED: {json.dumps(problems)}", file=sys.stderr)
        sys.exit(3)
    sys.exit(res.returncode)


if __name__ == "__main__":
    main()
