#!/usr/bin/env python3
"""Static checks on a WDL file's shell and embedded Python.

Three independent checks, because the first real execution of these workflows
will be on Terra and a syntax error there costs a full debug cycle:

1. Every task command block is shell-checked with `bash -n`, after WDL
   placeholders are replaced by a harmless literal, so unbalanced quotes,
   `if`/`fi` and heredoc terminators are caught here.
2. Every embedded Python heredoc is compiled with `ast.parse`.
3. The whole file is scanned for shell `${...}` nested inside a WDL
   `~{...}` expression, which Cromwell resolves as a WDL value and rejects.
   `bash -n` cannot see this, because placeholders are substituted first.
4. Every embedded Python heredoc is scanned for leaked `~{...}` placeholders.
   A quoted heredoc delimiter means WDL does not interpolate inside it, so a
   `~{...}` left in Python would reach the interpreter verbatim.
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PLACEHOLDER = re.compile(r"~\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}")
HEREDOC = re.compile(r"^\s*cat\s*>\s*(\S+)\s*<<\s*'([A-Za-z_][A-Za-z0-9_]*)'\s*$")


def command_blocks(text: str) -> list[tuple[str, str]]:
    """Return (task_name, command_body) for each <<< ... >>> command block."""
    blocks = []
    task = "<unknown>"
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^\s*task\s+(\w+)\s*\{", lines[i])
        if m:
            task = m.group(1)
        if re.match(r"^\s*command\s*<<<\s*$", lines[i]):
            body, i = [], i + 1
            while i < len(lines) and not re.match(r"^\s*>>>\s*$", lines[i]):
                body.append(lines[i])
                i += 1
            blocks.append((task, "\n".join(body)))
        i += 1
    return blocks


def python_heredocs(body: str) -> list[tuple[str, str]]:
    """Return (delimiter, python_source) for each quoted-heredoc Python block."""
    found = []
    lines = body.splitlines()
    i = 0
    while i < len(lines):
        m = HEREDOC.match(lines[i])
        if m and m.group(1).endswith(".py"):
            delim = m.group(2)
            src, i = [], i + 1
            while i < len(lines) and lines[i].strip() != delim:
                src.append(lines[i])
                i += 1
            found.append((delim, "\n".join(src)))
        i += 1
    return found



def nested_placeholders(text: str) -> list[tuple[int, str]]:
    """Find shell-style ${...} nested inside a WDL ~{...} expression.

    Inside a WDL string literal, ${...} is a WDL placeholder, not shell. So a
    line such as

        ~{if flag then "EXTRA=\"${EXTRA} --snp_vars\"" else "true"}

    makes Cromwell resolve EXTRA as a WDL value and reject the workflow with
    "Cannot lookup value 'EXTRA', it is never declared". This also applies to
    comments inside a command block, which are shell text and are still
    scanned for WDL placeholders. bash -n cannot catch this, because the
    placeholders are substituted away before the shell ever sees them.
    """
    hits = []
    lines = text.split("\n")
    i = 0
    while i < len(text):
        if text.startswith("~{", i):
            depth, j = 0, i + 1
            while j < len(text):
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            expr = text[i:j + 1]
            if "${" in expr:
                ln = text.count("\n", 0, i) + 1
                hits.append((ln, lines[ln - 1].strip()[:90]))
            i = j + 1
        else:
            i += 1
    return hits


def main(argv: list[str]) -> int:
    problems = 0
    for path in (Path(p) for p in argv[1:]):
        text = path.read_text()
        print(f"=== {path.name} ===")

        nested = nested_placeholders(text)
        if nested:
            problems += len(nested)
            for ln, src in nested:
                print(f"  FAIL line {ln}: shell ${{...}} nested in a WDL ~{{...}}: {src}")
        else:
            print("  ok   no shell ${...} nested inside a WDL ~{...}")
        for task, body in command_blocks(text):
            # bash -n on the placeholder-substituted body.
            shell = PLACEHOLDER.sub("PLACEHOLDER", body)
            with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
                fh.write(shell)
                tmp = fh.name
            proc = subprocess.run(["bash", "-n", tmp], capture_output=True, text=True)
            Path(tmp).unlink(missing_ok=True)
            if proc.returncode != 0:
                problems += 1
                print(f"  FAIL {task}: bash -n\n{proc.stderr.strip()}")
            else:
                print(f"  ok   {task}: bash -n")

            for delim, src in python_heredocs(body):
                n = len(src.splitlines())
                try:
                    ast.parse(src)
                except SyntaxError as exc:
                    problems += 1
                    print(f"  FAIL {task}/{delim}: {exc}")
                    continue
                leaked = PLACEHOLDER.findall(src)
                if leaked:
                    problems += 1
                    print(f"  FAIL {task}/{delim}: {len(leaked)} leaked placeholder(s): "
                          f"{leaked[:3]}")
                else:
                    print(f"  ok   {task}/{delim}: {n} lines, compiles, no interpolation")
    print(f"\n{problems} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
