"""The src() modules of each Snakefile rule cover the lincdisease modules its command uses, so that an edit of one
of them reruns the step."""
import ast
import re
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
MODULES = {p.stem for p in (PKG / "lincdisease").glob("*.py")} - {"__init__", "__main__"}


def _imports(module):
    tree = ast.parse((PKG / "lincdisease" / f"{module}.py").read_text())
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 1:
            out.update([node.module] if node.module else [a.name for a in node.names])
    return out


def _closure(modules):
    todo, seen = list(modules), set()
    while todo:
        m = todo.pop()
        if m not in seen:
            seen.add(m)
            todo.extend(_imports(m) - seen)
    return seen


def _used_by_commands():
    """cmd_<step> -> lincdisease modules it uses, directly or through the helper functions of __main__."""
    tree = ast.parse((PKG / "lincdisease" / "__main__.py").read_text())
    funcs = {f.name: f for f in tree.body if isinstance(f, ast.FunctionDef)}

    def uses(name, seen):
        out = set()
        for n in ast.walk(funcs[name]):
            if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name) and n.value.id in MODULES:
                out.add(n.value.id)
            elif isinstance(n, ast.Name) and n.id in funcs and n.id not in seen:
                seen.add(n.id)
                out |= uses(n.id, seen)
        return out
    return {name: uses(name, {name}) for name in funcs if name.startswith("cmd_")}


def test_src_covers_modules_used():
    used = _used_by_commands()
    snakefile = (PKG / "workflow" / "Snakefile").read_text()
    checked = 0
    for block in re.split(r"\nrule ", snakefile)[1:]:
        m = re.search(r"src\(([^)]*)\)", block)
        if not m:
            continue
        rule = block.split(":")[0]
        declared = _closure(["inputs", "resources"] + [x.strip().strip('"') for x in m.group(1).split(",")])
        assert _closure(used[f"cmd_{rule}"]) <= declared, f"rule {rule}: add {_closure(used[f'cmd_{rule}']) - declared}"
        checked += 1
    assert checked >= 15
