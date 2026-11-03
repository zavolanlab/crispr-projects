"""The code() fingerprint of each Snakefile rule covers the lincdisease modules its command uses, so that an edit of
one of them reruns the step; it ignores docstrings, comments and other steps' code (workflow/codehash.py)."""
import ast
import re
import shutil
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
MODULES = {p.stem for p in (PKG / "lincdisease").glob("*.py")} - {"__init__", "__main__"}
sys.path.insert(0, str(PKG / "workflow"))
import codehash  # noqa: E402


def _closure(modules):
    return codehash.module_closure(PKG / "lincdisease", modules)


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


def _rule_code_args():
    """rule -> the arguments of its code() call (step name first)."""
    out = {}
    for block in re.split(r"\nrule ", (PKG / "workflow" / "Snakefile").read_text())[1:]:
        m = re.search(r"code=code\(([^)]*)\)", block)
        if m:
            out[block.split(":")[0]] = [x.strip().strip('"') for x in m.group(1).split(",")]
    return out


def test_code_covers_modules_used():
    used = _used_by_commands()
    rules = _rule_code_args()
    for rule, (step, *modules) in rules.items():
        assert step == rule
        own = [step] if step in MODULES else []
        declared = _closure(own + modules)
        assert _closure(used[f"cmd_{rule}"]) <= declared, f"rule {rule}: add {_closure(used[f'cmd_{rule}']) - declared}"
    assert len(rules) >= 15


def _copy_package(tmp_path):
    d = tmp_path / "lincdisease"
    shutil.copytree(PKG / "lincdisease", d, ignore=shutil.ignore_patterns("__pycache__"))
    return d


def test_fingerprint_ignores_docstrings_comments_and_other_steps(tmp_path):
    d = _copy_package(tmp_path)
    before = {s: codehash.fingerprint(d, s, [m]) for s, m in (("hcc", "hcc"), ("weights", "weights"))}
    main = (d / "__main__.py").read_text()
    main = main.replace('def cmd_weights(args):\n', 'def cmd_weights(args):\n    """Changed docstring."""\n    # comment\n')
    main = main.replace('s.add_argument("--min-papers", type=int',
                        's.add_argument("--min-papers", type=int, metavar="N"')
    (d / "__main__.py").write_text(main)
    hcc_py = (d / "hcc.py").read_text()
    (d / "hcc.py").write_text(hcc_py.replace('"""', '"""Edited. ', 1))
    assert codehash.fingerprint(d, "hcc", ["hcc"]) == before["hcc"]          # docstring of its module, other step
    assert codehash.fingerprint(d, "weights", ["weights"]) != before["weights"]   # its option changed
    (d / "hcc.py").write_text(hcc_py + "\nX_NEW = 1\n")
    assert codehash.fingerprint(d, "hcc", ["hcc"]) != before["hcc"]


def test_fingerprint_follows_main_helpers(tmp_path):
    d = _copy_package(tmp_path)
    before = codehash.fingerprint(d, "scaling", ["expression"])
    main = (d / "__main__.py").read_text()
    assert "def gene_sums(" in main
    (d / "__main__.py").write_text(main.replace("def gene_sums(data_dir, projects, samples):\n",
                                                "def gene_sums(data_dir, projects, samples):\n    _ = 0\n"))
    assert codehash.fingerprint(d, "scaling", ["expression"]) != before
    assert codehash.fingerprint(d, "weights", ["weights"]) == codehash.fingerprint(PKG / "lincdisease", "weights",
                                                                                   ["weights"])
