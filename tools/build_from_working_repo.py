"""Build this submission repository from the working repository.

The working repository (default D:/thesis-research) is only READ, never written.
This script copies the code behind every thesis display item into a structure
that follows the thesis, and every other code/Markdown file into drafts/ with its
original path.

Copied thesis code is changed in three ways only, each checked afterwards:
  1. the hardcoded project root (D:/thesis-research) becomes PROJECT_ROOT_STR
     from thesis_research/config.py (Python) or PROJECT_ROOT (R), both read
     from the THESIS_PROJECT_ROOT environment variable;
  2. references to renamed scripts (file names, imports, Run: lines) follow
     the new names;
  3. the few explicit edits in SPECIAL, which fix imports between scripts that
     now sit in different folders.
Every changed line is written to tools/build_report.md for review.

Hand-written files (README.md, TABLE_S1.md, requirements.txt, pyproject.toml,
.gitignore, thesis_research/config.py, tools/) are never touched.

Run:  python tools/build_from_working_repo.py            (first build)
      python tools/build_from_working_repo.py --update   (re-sync after changes)
"""
import argparse
import ast
import difflib
import os
import pathlib
import re
import shutil
import subprocess
import sys

DST = pathlib.Path(__file__).resolve().parents[1]
ROOT_LITERAL = re.compile(r"D:[/\\]thesis-research")
CONFIG_NAMES = ("PROJECT_ROOT_STR", "L321_TX_FILE")

A1, A2, A3 = ("analysis/01_quality_control", "analysis/02_tumor_identification",
              "analysis/03_data_quality")
PD, DX, FR = f"{A3}/probe_detection", f"{A3}/contamination_decontx", f"{A3}/segmentation_fastreseg"

# old path (working repo) -> new path (this repo)
MAP = {
    # shared library (paths unchanged)
    **{p: p for p in """
        thesis_research/__init__.py
        thesis_research/pipeline/__init__.py
        thesis_research/pipeline/run_pipeline.py
        thesis_research/pipeline/position_plots.py
        thesis_research/pipeline/cell_qc_plots.py
        thesis_research/pipeline/filters.py
        thesis_research/pipeline/fov_qc_plots.py
        thesis_research/pipeline/sample_slice.py
        thesis_research/pipeline/utils.py
        thesis_research/pipeline/cell_type_annotation/__init__.py
        thesis_research/pipeline/cell_type_annotation/tumor_cells/__init__.py
        thesis_research/pipeline/cell_type_annotation/tumor_cells/identify_tumor_cells.py
        thesis_research/pipeline/cell_type_annotation/tumor_cells/refine_annotation_classifiers.py
        thesis_research/pipeline/cell_type_annotation/tumor_cells/classifiers.py
        thesis_research/utils/__init__.py
        thesis_research/utils/columns.py
        thesis_research/utils/constants.py
        thesis_research/utils/entity_type.py""".split()},
    # Chapter 2: quality control
    "thesis_plots/make_qc_metrics_overview_fig.py": f"{A1}/fig02_qc_metrics_overview.py",
    "thesis_plots/make_qc_count_threshold_fig.py": f"{A1}/fig03_count_threshold.py",
    # Chapter 2: tumor-cell identification
    "outputs/cell_annotation/annotate.R": f"{A2}/01_singler_annotate.R",
    "outputs/cell_annotation/convert_annotations_to_df.R": f"{A2}/02_singler_scores_to_table.R",
    "thesis_plots/make_fig_singler_sets.py": f"{A2}/fig04_07_singler_sets.py",
    "thesis_plots/make_table_singler_sweep.py": f"{A2}/table02_singler_threshold_sweep.py",
    "thesis_plots/figure_3_model_comparison.py": f"{A2}/fig08_classifier_comparison.py",
    "thesis_plots/xgb_default_sensitivity.py": f"{A2}/xgboost_default_sensitivity.py",
    "thesis_plots/stage1_threshold_sensitivity.py": f"{A2}/stage1_threshold_sensitivity.py",
    "thesis_plots/stage1_threshold_heldout.py": f"{A2}/stage1_threshold_heldout.py",
    "thesis_plots/stage1_anchor_floor_heldout.py": f"{A2}/stage1_anchor_floor_heldout.py",
    "thesis_plots/figure_4_spatial_refinement.py": f"{A2}/fig09_14_spatial_refinement.py",
    "thesis_plots/control_specificity_validation.py": f"{A2}/table03_heldout_specificity.py",
    "thesis_plots/final_xgboost_refinement.py": f"{A2}/fig15_final_tumor_calls.py",
    # Chapter 3: data-quality assessment
    "thesis_plots/update_tumor_prediction_cache.py": f"{A3}/write_final_tumor_calls.py",
    "thesis_plots/make_detection_reliability_6slice.py": f"{PD}/table05_detection_strength.py",
    "thesis_plots/make_dq_fig1_detection.py": f"{PD}/fig16_probe_detection.py",
    "thesis_plots/rerun_chapter3_final_calls.py": f"{PD}/run_on_final_tumor_calls.py",
    "thesis_plots/make_dq_fig_reporter.py": f"{A3}/fig17_tdtomato_prevalence.py",
    "thesis_plots/make_dq_fig_lyve1.py": f"{A3}/fig18_lyve1_bam_markers.py",
    "score_genes/write_decontx_clusters.py": f"{DX}/write_cluster_labels.py",
    "score_genes/run_decontx_correct.py": f"{DX}/run_decontx.py",
    "score_genes/run_decontx.R": f"{DX}/decontx_model.R",
    "thesis_plots/make_dq_fig_decontx.py": f"{DX}/fig19_decontx_before_after.py",
    "thesis_plots/decontx_partition_sensitivity.py": f"{DX}/partition_sensitivity.py",
    "agents/segmentation/03_filter_tx.py": f"{FR}/01_extract_reference_fovs.py",
    "agents/segmentation/fastreseg/01_prep_inputs.py": f"{FR}/02_prepare_reference_inputs.py",
    "agents/segmentation/fastreseg/02_run_fastreseg.R": f"{FR}/03_run_fastreseg_reference_fovs.R",
    "agents/segmentation/fastreseg/04_export_scores.R": f"{FR}/04_export_reference_scores.R",
    "agents/segmentation/fastreseg/slice1_01_prep_inputs.py": f"{FR}/05_prepare_slice1_inputs.py",
    "agents/segmentation/fastreseg/slice1_02_run_fastreseg.R": f"{FR}/06_run_fastreseg_slice1.R",
    "agents/segmentation/fastreseg/slice1_04_evaluate_qcpassed.py": f"{FR}/07_evaluate_slice1.py",
    "agents/segmentation/fastreseg/slice1_05_sparsity.py": f"{FR}/08_sparsity_slice1.py",
    "agents/segmentation/fastreseg/slice1_06_figure20.py": f"{FR}/fig20_reference_heatmap.py",
    "agents/segmentation/fastreseg/06_figure_A2.py": f"{FR}/fig20_heatmap_base.py",
    "agents/segmentation/fastreseg/slice1_07_figS1_groups.py": f"{FR}/figS1_cell_groups.py",
    "thesis_plots/make_dq_fig_tumor_spatial.py": f"{A3}/fig21_tumor_spatial.py",
}
# committed input data copied with the code that reads it
DATA = {f"score_genes/decontx_clusters/slice_{i}_clusters.csv": f"{DX}/decontx_clusters/slice_{i}_clusters.csv"
        for i in range(1, 7)}
DATA["score_genes/decontx_clusters/sensitivity/slice_1_leiden5_clusters.csv"] = \
    f"{DX}/decontx_clusters/sensitivity/slice_1_leiden5_clusters.csv"

# explicit edits, applied to the original text: new path -> [(old, new), ...]
_TO_A2 = 'sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "02_tumor_identification"))'
SPECIAL = {
    "thesis_research/utils/constants.py": [
        ("PROJECT_ROOT = Path(__file__).resolve().parents[2]",
         "from thesis_research.config import PROJECT_ROOT  # set by THESIS_PROJECT_ROOT"),
    ],
    f"{A3}/write_final_tumor_calls.py": [
        ("sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))", _TO_A2),
    ],
    f"{A3}/fig21_tumor_spatial.py": [
        ("from final_xgboost_refinement import",
         "import pathlib  # noqa: E402\nimport sys  # noqa: E402\n" + _TO_A2
         + "\nfrom final_xgboost_refinement import"),
    ],
    f"{DX}/partition_sensitivity.py": [
        ('sys.path.insert(0, os.path.join(ROOT, "score_genes"))',
         "sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))"),
        ('CLUSTERS = os.path.join(ROOT, "score_genes", "decontx_clusters")',
         'CLUSTERS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "decontx_clusters")'),
    ],
    f"{PD}/run_on_final_tumor_calls.py": [
        ("OUT_DIR = HERE / OUT_SUB",
         'OUT_DIR = pathlib.Path(PROJECT_ROOT_STR) / "thesis_plots" / OUT_SUB'),
        ("make_detection_reliability_6slice.py (Table 4)",
         "make_detection_reliability_6slice.py (Table 5)"),
    ],
    f"{FR}/06_run_fastreseg_slice1.R": [
        ("# Run: Rscript agents/segmentation/fastreseg/02_run_fastreseg.R",
         "# Run: Rscript agents/segmentation/fastreseg/slice1_02_run_fastreseg.R"),
    ],
    f"{FR}/01_extract_reference_fovs.py": [
        ('TX = (r"D:\\20251214_CosMx_ReuvenStein\\20251214_CosMx_ReuvenStein.tar\\Analysis"\n'
         '      r"\\L321__1__31_12_2025_12_32_59_204\\flatFiles\\L321\\L321_tx_file.csv")',
         "TX = L321_TX_FILE"),
    ],
    f"{FR}/05_prepare_slice1_inputs.py": [
        ('TX = pathlib.Path("D:/20251214_CosMx_ReuvenStein/20251214_CosMx_ReuvenStein.tar/Analysis/"\n'
         '                  "L321__1__31_12_2025_12_32_59_204/flatFiles/L321/L321_tx_file.csv")',
         "TX = pathlib.Path(L321_TX_FILE)"),
    ],
}
for _f in (f"{FR}/03_run_fastreseg_reference_fovs.R", f"{FR}/04_export_reference_scores.R",
           f"{FR}/06_run_fastreseg_slice1.R"):
    SPECIAL.setdefault(_f, []).append(('.libPaths(c("D:/R-libs/fastreseg", .libPaths()))',
                                       '.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))'))

# never copied anywhere: local AI-assistant instruction files, not research material
SKIP = {"CLAUDE.md", "agents/AGENT_CONTEXT.md", "agents/RESEARCH_PLAN.md"}
DRAFT_EXT = (".py", ".R", ".md")
HANDWRITTEN = {"README.md", "TABLE_S1.md", "requirements.txt", "pyproject.toml", ".gitignore",
               "thesis_research/config.py", "drafts/README.md"}

# ---------------------------------------------------------------- transforms
NEW_NAME = {pathlib.Path(o).name: pathlib.Path(n).name for o, n in MAP.items()
            if pathlib.Path(o).name != pathlib.Path(n).name}
NEW_STEM = {o[:-3]: n[:-3] for o, n in NEW_NAME.items() if o.endswith(".py")}
OLD_PATHS = {o: n for o, n in MAP.items() if o != n}
TOKEN = re.compile(r"(?<![\w.-])((?:[\w-]+/)*[\w-]+\.(?:py|R))\b")


def rename_refs(text):
    def sub(m):
        tok = m.group(1)
        if tok in OLD_PATHS:
            return OLD_PATHS[tok]
        if "/" not in tok and tok in NEW_NAME:
            return NEW_NAME[tok]
        return tok
    text = TOKEN.sub(sub, text)
    for old, new in NEW_STEM.items():
        text = re.sub(rf"(?m)^(\s*(?:from|import)\s+){re.escape(old)}\b", rf"\g<1>{new}", text)
    return text


PY_LIT = re.compile(r"""(?P<pre>\b[rRfF]{1,2})?(?P<q>["'])(?P<body>[^"'\n]*?D:[/\\]thesis-research[^"'\n]*?)(?P=q)""")
R_LIT = re.compile(r'"(?P<body>[^"\n]*?D:/thesis-research[^"\n]*?)"')


def py_root(text):
    def sub(m):
        pre, q = m.group("pre") or "", m.group("q")
        body = m.group("body")
        if "\\" in body and "r" not in pre.lower():
            raise ValueError(f"non-raw backslash path literal: {m.group(0)}")
        body = body.replace("\\", "/")
        if "f" in pre.lower():
            return pre + q + ROOT_LITERAL.sub("{PROJECT_ROOT_STR}", body) + q
        parts = ROOT_LITERAL.split(body)
        out = []
        for i, part in enumerate(parts):
            if part:
                out.append(pre + q + part + q)
            if i < len(parts) - 1:
                out.append("PROJECT_ROOT_STR")
        return " + ".join(out)
    return PY_LIT.sub(sub, text)


def r_root(text):
    def sub(m):
        rest = m.group("body").split("D:/thesis-research", 1)[1]
        return f'paste0(PROJECT_ROOT, "{rest}")' if rest else "PROJECT_ROOT"
    return R_LIT.sub(sub, text)


def add_config_import(text):
    names = [n for n in CONFIG_NAMES if re.search(rf"\b{n}\b", text)
             and not re.search(rf"import[^\n]*\b{n}\b", text)]
    if not names:
        return text
    tree = ast.parse(text)
    body = tree.body
    i = 0
    if body and isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], "value", None), ast.Constant):
        i = 1
    while i < len(body) and isinstance(body[i], ast.ImportFrom) and body[i].module == "__future__":
        i += 1
    line = body[i].lineno - 1 if i < len(body) else len(text.splitlines())
    lines = text.splitlines(keepends=True)
    lines.insert(line, f"from thesis_research.config import {', '.join(names)}  # noqa: E402\n")
    return "".join(lines)


R_HEADER = ('PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")\n'
            'if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")\n')


def transform(old, new, text):
    for a, b in SPECIAL.get(new, []):
        if text.count(a) != 1:
            raise ValueError(f"{new}: special edit target found {text.count(a)} times: {a[:60]!r}")
        text = text.replace(a, b)
    text = rename_refs(text)
    if new.endswith(".py"):
        text = py_root(text)
        text = add_config_import(text)
    else:
        if "D:/thesis-research" in text and re.search(r'"[^"\n]*D:/thesis-research', text):
            text = r_root(text)
            lines = text.splitlines(keepends=True)
            k = 0
            while k < len(lines) and lines[k].lstrip().startswith("#"):
                k += 1
            lines.insert(k, R_HEADER)
            text = "".join(lines)
    # output folders live under automation/outputs/ in the project data folder
    # (a link to the working repo's agents/outputs/)
    text = text.replace("agents/outputs", "automation/outputs")
    # remaining references to working-repo code point at its copy in drafts/
    text = text.replace("agents/segmentation/", "drafts/automation/segmentation/")
    # what is left of the old root is prose in docstrings/comments (check() verifies)
    return text.replace("D:/thesis-research", "$THESIS_PROJECT_ROOT")


# ------------------------------------------------------------------- checks
def code_strings(text):
    """String constants in Python code, excluding docstrings."""
    tree = ast.parse(text)
    docs = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant):
                docs.add(id(node.body[0].value))
    return [n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs]


def check(new, text):
    problems = []
    if new.endswith(".py"):
        for s in code_strings(text):
            if (ROOT_LITERAL.search(s) or "$THESIS_PROJECT_ROOT" in s or "R-libs" in s
                    or "20251214_CosMx" in s):
                problems.append(f"hardcoded path left in code: {s[:80]!r}")
        for old in NEW_STEM:
            if re.search(rf"(?m)^\s*(from|import)\s+{re.escape(old)}\b", text):
                problems.append(f"old module import left: {old}")
    else:
        code = "\n".join(l.split("#", 1)[0] for l in text.splitlines())
        if ROOT_LITERAL.search(code) or "$THESIS_PROJECT_ROOT" in code or "R-libs" in code:
            problems.append("hardcoded path left in R code")
    return problems


# --------------------------------------------------------------------- main
def git_files(src, args):
    out = subprocess.run(["git", "-C", str(src)] + args, capture_output=True, text=True, check=True)
    return [l for l in out.stdout.splitlines() if l]


def write(path, data, update, written):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() == data:
            written.append(path)
            return
        if not update:
            raise FileExistsError(f"{path} exists with different content; use --update")
    path.write_bytes(data)
    written.append(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="D:/thesis-research")
    ap.add_argument("--update", action="store_true")
    a = ap.parse_args()
    src = pathlib.Path(a.src).resolve()
    if src == DST or DST in src.parents or src in DST.parents:
        raise SystemExit("source and destination must be separate folders")

    tracked = git_files(src, ["ls-files"])
    untracked = git_files(src, ["ls-files", "--others", "--exclude-standard"])
    for old in list(MAP) + list(DATA):
        if old not in tracked:
            raise SystemExit(f"mapped file not tracked in the working repo: {old}")

    written, report, problems = [], ["# Build report\n", f"Source: `{src}`\n"], []
    for old, new in MAP.items():
        text = (src / old).read_text(encoding="utf-8")
        out = transform(old, new, text)
        problems += [f"{new}: {p}" for p in check(new, out)]
        if new.endswith(".py"):
            try:
                ast.parse(out)
            except SyntaxError as e:
                problems.append(f"{new}: syntax error {e}")
        write(DST / new, out.encode("utf-8"), a.update, written)
        diff = list(difflib.unified_diff(text.splitlines(), out.splitlines(), old, new, n=0, lineterm=""))
        report.append(f"\n## `{old}` -> `{new}`\n")
        report.append("\n```diff\n" + "\n".join(diff[2:]) + "\n```\n" if diff else "\nunchanged\n")
    for old, new in DATA.items():
        write(DST / new, (src / old).read_bytes(), a.update, written)

    thesis = set(MAP) | set(DATA)
    drafts = sorted(f for f in set(tracked) | set(untracked)
                    if f.endswith(DRAFT_EXT) and f not in thesis and f not in SKIP)
    for old in drafts:
        rel = "automation/" + old[len("agents/"):] if old.startswith("agents/") else old
        new = "drafts/WORKING_REPO_README.md" if old == "README.md" else f"drafts/{rel}"
        write(DST / new, (src / old).read_bytes(), a.update, written)

    report.insert(2, f"\n{len(MAP)} thesis files, {len(DATA)} data files, {len(drafts)} draft files.\n")
    (DST / "tools" / "build_report.md").write_text("".join(report), encoding="utf-8")
    managed = {p.resolve() for p in written}
    stale = [p for p in DST.rglob("*") if p.is_file() and p.resolve() not in managed
             and ".git" not in p.parts and p.relative_to(DST).as_posix() not in HANDWRITTEN
             and p.relative_to(DST).parts[0] != "tools"]
    print(f"thesis files: {len(MAP)}, data files: {len(DATA)}, drafts: {len(drafts)}")
    if stale:
        print("files not produced by this build (left in place):", *stale, sep="\n  ")
    if problems:
        print("PROBLEMS:", *problems, sep="\n  ")
        sys.exit(1)
    print("checks passed; changed lines listed in tools/build_report.md")


if __name__ == "__main__":
    main()
