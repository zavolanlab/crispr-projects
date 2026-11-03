"""
Code fingerprint of a workflow step (Snakefile code()): a hash of the code that computes the step's outputs, so
that a step reruns when that code changes, and only then.

The fingerprint covers
- the lincdisease modules the step uses, with everything they import (from . import x), and the input readers
  (inputs, resources);
- the parts of __main__.py that belong to the step: cmd_<step> and the top-level functions and constants it uses
  (directly or through each other), and the step's command-line options (its add_parser block in main(), with the
  helpers that block calls), plus the options shared by all steps.
It is computed on the syntax tree without docstrings, so edits of comments, docstrings or formatting, and edits of
other steps' code in __main__.py, do not change it.
"""
import ast
import hashlib
from pathlib import Path


def _strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                node.body = body[1:] or [ast.Pass()]
    return tree


def _parse(path):
    return _strip_docstrings(ast.parse(Path(path).read_text()))


def imports(modules_dir, module):
    """lincdisease modules imported by `module` (from . import x / from .x import y, also inside functions)."""
    out = set()
    for node in ast.walk(ast.parse((Path(modules_dir) / f"{module}.py").read_text())):
        if isinstance(node, ast.ImportFrom) and node.level == 1:
            out.update([node.module] if node.module else [a.name for a in node.names])
    return out


def module_closure(modules_dir, modules):
    """`modules` with everything they import, and inputs and resources."""
    todo, seen = ["inputs", "resources", *modules], set()
    while todo:
        m = todo.pop()
        if m not in seen:
            seen.add(m)
            todo.extend(imports(modules_dir, m) - seen)
    return seen


def _step_named(stmt):
    """Step name of a statement that creates a subcommand (sub.add_parser("x", ...), _cohort_step(sub, "x", ...))."""
    for call in (n for n in ast.walk(stmt) if isinstance(n, ast.Call)):
        f = call.func
        name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
        if name in ("add_parser", "_cohort_step"):
            names = [a.value for a in call.args if isinstance(a, ast.Constant) and isinstance(a.value, str)]
            if names:
                return names[0]
    return None


def _parser_blocks(main):
    """main() statements per step (from the statement creating its subcommand to its set_defaults) and the other
    statements (options shared by all steps, parsing and dispatch)."""
    blocks, shared, current = {}, [], None
    for stmt in main.body:
        current = _step_named(stmt) or current
        if current is None:
            shared.append(stmt)
            continue
        blocks.setdefault(current, []).append(stmt)
        if any(isinstance(n, ast.Attribute) and n.attr == "set_defaults" for n in ast.walk(stmt)):
            current = None
    return blocks, shared


def main_parts(main_py, step):
    """The syntax trees of __main__.py that belong to `step` (module docstring)."""
    tree = _parse(main_py)
    defs = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defs[node.name] = node
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for t in targets:
                for n in ast.walk(t):
                    if isinstance(n, ast.Name):
                        defs[n.id] = node
    blocks, shared = _parser_blocks(defs["main"])
    if f"cmd_{step}" not in defs:
        raise KeyError(f"__main__.py has no cmd_{step}")
    roots = [defs[f"cmd_{step}"], *blocks.get(step, []), *shared]
    parts, seen, todo = [], set(), list(roots)
    while todo:
        node = todo.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        parts.append(node)
        todo.extend(defs[n.id] for n in ast.walk(node)
                    if isinstance(n, ast.Name) and n.id in defs and n.id != "main" and id(defs[n.id]) not in seen)
    return parts


def fingerprint(modules_dir, step, modules):
    """Hash of the step's code (module docstring); `modules`: the lincdisease modules the step uses."""
    modules_dir = Path(modules_dir)
    h = hashlib.sha1()
    for m in sorted(module_closure(modules_dir, modules)):
        h.update(m.encode() + ast.dump(_parse(modules_dir / f"{m}.py")).encode())
    for d in sorted(ast.dump(n) for n in main_parts(modules_dir / "__main__.py", step)):
        h.update(d.encode())
    return h.hexdigest()[:16]
