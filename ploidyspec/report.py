"""
Self-contained HTML report per species: embeds already-generated plots
(base64 data URIs, no external file references) and renders the key TSV
tables. Degrades gracefully when a stage hasn't been run -- the report must
never error on a partial pipeline, since most species in the panel don't
have every stage run (te-markers in particular is opt-in).
"""

import base64
import csv
import glob
import html
import os

from .common import homeologs_dir, matrix_dir, subgenomes_dir, windowed_dir

RANKED_CANDIDATES_LIMIT = 20


def _read_tsv(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _embed_image(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{data}"


def collect_species_data(outdir):
    """Gathers whatever's present per stage directory into a plain dict of
    tables (list-of-dict, or None if the file doesn't exist) and images
    (data URI, or None). Never raises on missing files/stages -- the report
    must work on however much of the pipeline has actually been run."""
    mdir = matrix_dir(outdir)
    wdir = windowed_dir(outdir)
    hdir = homeologs_dir(outdir)
    sdir = subgenomes_dir(outdir)

    data = {
        "species": os.path.basename(os.path.normpath(outdir)),
        "matrix": {
            "ploidy_summary": _read_tsv(os.path.join(mdir, "ploidy_summary.tsv")),
            "homologous_chromosomes": _read_tsv(
                os.path.join(mdir, "homologous_chromosomes.tsv")
            ),
            "heatmap": _embed_image(
                os.path.join(mdir, "whole_chrom_distance_heatmap.png")
            ),
            "heatmap_contrast": _embed_image(
                os.path.join(mdir, "whole_chrom_distance_heatmap_contrast.png")
            ),
            "k_resolution": _embed_image(
                os.path.join(mdir, "k_resolution_diagnostic.png")
            ),
        },
        "windowed": {
            "overview": _embed_image(
                os.path.join(wdir, "windowed_genome_overview.png")
            ),
            "overview_heatmap": _embed_image(
                os.path.join(wdir, "windowed_genome_overview_heatmap.png")
            ),
        },
        "homeologs": {
            "homeolog_pairs": _read_tsv(os.path.join(hdir, "homeolog_pairs.tsv")),
            "ploidy_ancestry_summary": _read_tsv(
                os.path.join(hdir, "ploidy_ancestry_summary.tsv")
            ),
            "ranked_candidates": _read_tsv(
                os.path.join(hdir, "homeolog_candidates_ranked.tsv")
            ),
            "homeolog_pairs_plot": _embed_image(
                os.path.join(hdir, "homeolog_pairs.png")
            ),
        },
        "subgenomes": None,
    }

    auto_allo = _read_tsv(os.path.join(sdir, "auto_allo_index.tsv"))
    if auto_allo is not None:
        windowed_pngs = sorted(
            glob.glob(os.path.join(sdir, "te_markers_windowed_*.png"))
        )
        data["subgenomes"] = {
            "auto_allo_index": auto_allo,
            "windows_summary": _read_tsv(
                os.path.join(sdir, "subgenome_windows_summary.tsv")
            ),
            "windowed_plots": [_embed_image(p) for p in windowed_pngs],
        }

    return data


def _table_html(rows, max_rows=None):
    if not rows:
        return '<p class="empty">No data.</p>'
    columns = list(rows[0].keys())
    limited = rows if max_rows is None else rows[:max_rows]
    header = "".join(f"<th>{html.escape(c)}</th>" for c in columns)
    body_rows = []
    for r in limited:
        cells = "".join(f"<td>{html.escape(str(r.get(c, '')))}</td>" for c in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    out = f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"
    if max_rows is not None and len(rows) > max_rows:
        out += f'<p class="note">showing top {max_rows} of {len(rows)} rows</p>'
    return out


def _img_html(data_uri, alt):
    if not data_uri:
        return ""
    return f'<img src="{data_uri}" alt="{html.escape(alt)}">'


def render_report_html(data):
    """Pure: takes a dict shaped like collect_species_data's output and
    returns a complete HTML document as a string. No file I/O -- testable
    with hand-constructed dicts."""
    species = html.escape(data.get("species", ""))
    sections = []

    m = data.get("matrix") or {}
    if m.get("ploidy_summary") is not None or m.get("heatmap"):
        sections.append(
            f"""<section>
  <h2>Whole-chromosome matrix</h2>
  {_img_html(m.get('heatmap'), 'distance heatmap')}
  {_img_html(m.get('heatmap_contrast'), 'contrast-stretched distance heatmap')}
  {_img_html(m.get('k_resolution'), 'k-resolution diagnostic')}
  <h3>Ploidy summary</h3>
  {_table_html(m.get('ploidy_summary'))}
  <h3>Homologous chromosome pairs</h3>
  {_table_html(m.get('homologous_chromosomes'))}
</section>"""
        )
    else:
        sections.append(
            '<section><h2>Whole-chromosome matrix</h2>'
            '<p class="empty">Not run yet.</p></section>'
        )

    w = data.get("windowed") or {}
    if w.get("overview"):
        sections.append(
            f"""<section>
  <h2>Windowed divergence</h2>
  {_img_html(w.get('overview'), 'windowed genome overview')}
  {_img_html(w.get('overview_heatmap'), 'windowed genome overview heatmap')}
</section>"""
        )
    else:
        sections.append(
            '<section><h2>Windowed divergence</h2>'
            '<p class="empty">Not run yet.</p></section>'
        )

    h = data.get("homeologs") or {}
    if h.get("homeolog_pairs") is not None or h.get("ranked_candidates") is not None:
        sections.append(
            f"""<section>
  <h2>Ancient homeolog pairing</h2>
  {_img_html(h.get('homeolog_pairs_plot'), 'homeolog pairs')}
  <h3>Accepted pairs</h3>
  {_table_html(h.get('homeolog_pairs'))}
  <h3>Ploidy / ancestry summary</h3>
  {_table_html(h.get('ploidy_ancestry_summary'))}
  <h3>All candidates (ranked by distance)</h3>
  {_table_html(h.get('ranked_candidates'), max_rows=RANKED_CANDIDATES_LIMIT)}
</section>"""
        )
    else:
        sections.append(
            '<section><h2>Ancient homeolog pairing</h2>'
            '<p class="empty">Not run yet.</p></section>'
        )

    s = data.get("subgenomes")
    if s:
        plot_imgs = "".join(
            _img_html(p, "te-marker windowed track")
            for p in s.get("windowed_plots", [])
            if p
        )
        sections.append(
            f"""<section>
  <h2>Subgenomes / auto-allo index</h2>
  <h3>Auto/allo index</h3>
  {_table_html(s.get('auto_allo_index'))}
  <h3>Window assignment summary</h3>
  {_table_html(s.get('windows_summary'))}
  {plot_imgs}
</section>"""
        )
    else:
        sections.append(
            '<section><h2>Subgenomes / auto-allo index</h2>'
            '<p class="empty">te-markers not run for this species.</p></section>'
        )

    body = "\n".join(sections)
    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>ploidyspec report: {species}</title>
<style>
{_CSS}
</style>
</head>
<body>
<h1>ploidyspec report: {species}</h1>
{body}
</body>
</html>
"""


_CSS = """
body { font-family: -apple-system, sans-serif; max-width: 1000px; margin: 2rem auto; padding: 0 1rem; color: #222; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.5rem; }
section { margin-bottom: 2.5rem; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem; margin-bottom: 1rem; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
th { background: #f0f0f0; }
img { max-width: 100%; display: block; margin: 0.5rem 0; }
.empty { color: #888; font-style: italic; }
.note { color: #888; font-size: 0.85rem; }
"""


def generate_report(outdir):
    data = collect_species_data(outdir)
    html_doc = render_report_html(data)
    path = os.path.join(outdir, "report.html")
    with open(path, "w") as f:
        f.write(html_doc)
    return path
