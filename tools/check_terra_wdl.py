#!/usr/bin/env python3
"""Reject file-writing functions in workflow expressions, including imported workflows."""

import argparse
from pathlib import Path

import WDL


def workflow_write_errors(document):
    errors = []

    def visit(node):
        if isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("write_"):
            errors.append(f"{node.pos.uri}:{node.pos.line}:{node.pos.column}: {node.function_name} is not allowed at workflow scope")
        for child in node.children:
            visit(child)

    if document.workflow:
        visit(document.workflow)
    for imported in document.imports:
        errors.extend(workflow_write_errors(imported.doc))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wdls", nargs="+", type=Path)
    args = parser.parse_args()
    errors = []
    for path in args.wdls:
        document = WDL.load(str(path))
        errors.extend(workflow_write_errors(document))
    if errors:
        print("\n".join(errors))
        return 1
    print("Terra scope check passed: no workflow-scope file-writing functions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
