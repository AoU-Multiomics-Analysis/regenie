#!/usr/bin/env python3
"""Check Terra workflow expressions and cloud task output collection."""

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


def manifest_output_errors(document):
    """Reject File outputs obtained by reading paths from a generated manifest."""
    errors = []

    def contains_file(kind):
        if isinstance(kind, WDL.Type.File):
            return True
        if isinstance(kind, WDL.Type.StructInstance):
            return any(contains_file(member) for member in kind.members.values())
        return any(contains_file(parameter) for parameter in kind.parameters)

    for task in document.tasks:
        declarations = {decl.name: decl for decl in [*(task.inputs or []), *task.postinputs, *task.outputs]}

        def reads_manifest(node, seen):
            if isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("read_"):
                return True
            if isinstance(node, WDL.Expr.Ident) and node.name in declarations and node.name not in seen:
                expr = declarations[node.name].expr
                if expr is not None and reads_manifest(expr, seen | {node.name}):
                    return True
            return any(reads_manifest(child, seen) for child in node.children)

        for output in task.outputs:
            if contains_file(output.type) and reads_manifest(output.expr, {output.name}):
                errors.append(f"{output.pos.uri}:{output.pos.line}:{output.pos.column}: File outputs cannot come from manifest contents; use explicit File or glob outputs")
    for imported in document.imports:
        errors.extend(manifest_output_errors(imported.doc))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wdls", nargs="+", type=Path)
    args = parser.parse_args()
    errors = []
    for path in args.wdls:
        document = WDL.load(str(path))
        errors.extend(workflow_write_errors(document))
        errors.extend(manifest_output_errors(document))
    if errors:
        print("\n".join(errors))
        return 1
    print("Terra checks passed: no workflow-scope file writing or manifest-derived File outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
