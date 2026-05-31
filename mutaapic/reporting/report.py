import os
import shutil
import html
import json
import pandas as pd


def _df_to_html(df):
    try:
        return df.to_html(index=False, classes='dataframe')
    except Exception:
        return f"<pre>{str(df)}</pre>"


def _interpret_structure_change(comparison_results):
    tm_score = comparison_results.get("tm_score")
    rmsd = comparison_results.get("rmsd")

    try:
        tm_score = float(tm_score) if tm_score is not None else None
    except Exception:
        tm_score = None

    try:
        rmsd = float(rmsd) if rmsd is not None else None
    except Exception:
        rmsd = None

    if tm_score is None and rmsd is None:
        return (
            "Preliminary interpretation: unclear. "
            "The report does not have enough structural evidence to estimate whether the mutation is likely tolerated."
        )

    if tm_score is not None and rmsd is not None:
        if tm_score >= 0.90 and rmsd <= 2.0:
            return (
                "Preliminary interpretation: probably tolerated / lower risk. "
                "The structures are highly similar, so the mutation may be less likely to disrupt the protein, "
                "but this is not a functional guarantee."
            )
        if tm_score >= 0.75 and rmsd <= 4.0:
            return (
                "Preliminary interpretation: uncertain, with moderate structural similarity. "
                "The mutation might be tolerated, but there is still a meaningful chance of structural or functional impact."
            )
        return (
            "Preliminary interpretation: possible structural impact. "
            "The alignment suggests noticeable change, so the mutation may affect stability or function, "
            "but additional evidence would be needed to be confident."
        )

    if tm_score is not None:
        if tm_score >= 0.90:
            return (
                "Preliminary interpretation: likely tolerated, based mainly on TM-score. "
                "This is a soft signal only and does not prove the protein function is unchanged."
            )
        if tm_score >= 0.75:
            return (
                "Preliminary interpretation: uncertain. The TM-score suggests some structural preservation, "
                "but the mutation could still have functional consequences."
            )
        return (
            "Preliminary interpretation: possible impact. The TM-score is not especially high, so the mutation may alter structure or function."
        )

    if rmsd is not None:
        if rmsd <= 2.0:
            return (
                "Preliminary interpretation: likely tolerated, based mainly on a low RMSD. "
                "This is only a heuristic and should not be treated as a definite functional prediction."
            )
        if rmsd <= 4.0:
            return (
                "Preliminary interpretation: uncertain. The RMSD is moderate, so the mutation may or may not change protein function."
            )
        return (
            "Preliminary interpretation: possible impact. The RMSD is relatively high, which can be consistent with a structural effect."
        )

    return (
        "Preliminary interpretation: unclear. The available structural metrics are incomplete, so no cautious assessment can be made."
    )


def generate_report(out_dir: str,
                    orig_pdb: str,
                    mut_pdb: str,
                    comparison_results: dict,
                    custom_results_df=None,
                    custom_structure_ids=None,
                    pdb_results_df=None,
                    af_results_df=None,
                    custom_downloads=None,
                    pdb_downloads=None,
                    af_downloads=None):
    """Generate a simple HTML report summarizing the analysis results.

    The report will be written to `out_dir/report/report.html` and will include
    the TM-align comparison, Foldseek search tables (if present) and links to
    downloaded structures.
    """
    report_dir = os.path.join(out_dir, "report")
    os.makedirs(report_dir, exist_ok=True)

    # Copy PDBs into report folder for convenient linking
    def _copy_to_report(path):
        if not path:
            return None
        try:
            base = os.path.basename(path)
            dst = os.path.join(report_dir, base)
            shutil.copy(path, dst)
            return base
        except Exception:
            return None

    orig_link = _copy_to_report(orig_pdb)
    mut_link = _copy_to_report(mut_pdb)

    parts = []
    
    parts.append("<html><head><meta charset='utf-8'><title>MutAAP-IC Report</title>")
    parts.append("<style>body{font-family:Arial,Helvetica,sans-serif;padding:18px}h1,h2{color:#2b2b2b}table.dataframe{border-collapse:collapse;width:100%}table.dataframe th, table.dataframe td{border:1px solid #ddd;padding:8px;text-align:left}tr:nth-child(even){background:#f9f9f9}</style>")
    parts.append("</head><body>")
    parts.append("<h1>MutAAP-IC Report</h1>")

    parts.append("<h2>Inputs</h2>")
    parts.append("<ul>")
    if orig_link:
        parts.append(f"<li>Original structure: <a href=\"{orig_link}\">{orig_link}</a></li>")
    else:
        parts.append(f"<li>Original structure: {orig_pdb}</li>")
    if mut_link:
        parts.append(f"<li>Mutant structure: <a href=\"{mut_link}\">{mut_link}</a></li>")
    else:
        parts.append(f"<li>Mutant structure: {mut_pdb}</li>")
    parts.append("</ul>")

 
    parts.append("<h2>TM-align Comparison</h2>")
    if comparison_results:
        parts.append("<h3>Preliminary interpretation of the mutations</h3>")
        parts.append(f"<p>{html.escape(_interpret_structure_change(comparison_results))}</p>")

        parts.append("<table>")
        # Only include tm_score, rmsd and alignment; omit raw_output
        tm = comparison_results.get("tm_score")
        if tm is not None:
            try:
                tm_val = f"{float(tm):.4f}"
            except Exception:
                tm_val = str(tm)
            parts.append(f"<tr><th style='text-align:left;padding-right:12px'>tm_score</th><td>{html.escape(tm_val)}</td></tr>")

        rmsd = comparison_results.get("rmsd")
        if rmsd is not None:
            try:
                rmsd_val = f"{float(rmsd):.3f}"
            except Exception:
                rmsd_val = str(rmsd)
            parts.append(f"<tr><th style='text-align:left;padding-right:12px'>rmsd</th><td>{html.escape(rmsd_val)}</td></tr>")

        alignment = comparison_results.get("alignment")
        if alignment:
            # alignment expected as dict with seq1, similarity, seq2
            if isinstance(alignment, dict):
                seq1 = alignment.get("seq1") or ""
                sim = alignment.get("similarity") or ""
                seq2 = alignment.get("seq2") or ""
                align_text = "\n".join([ln for ln in (seq1, sim, seq2) if ln])
            else:
                align_text = str(alignment)
            if len(align_text) > 1000:
                align_text = align_text[:1000] + "\n... (truncated)"
            parts.append(f"<tr><th style='text-align:left;padding-right:12px'>alignment</th><td><pre>{html.escape(align_text)}</pre></td></tr>")

        parts.append("</table>")
    else:
        parts.append("<p>No comparison results available.</p>")

    # Combined 3Dmol viewer showing original (blue) and mutant (red) structures
    parts.append("<h2>Structure View</h2>")
    parts.append("<div style='width:100%;margin-bottom:24px'>"
                "<div id='structure_viewer' style='width:100%;height:400px;max-height:400px;overflow:hidden;border:1px solid #ddd;position:relative'></div>"
                "<div style='margin-top:8px'>"
                "<span style='display:inline-block;width:16px;height:16px;background:#007bff;margin-right:6px;vertical-align:middle'></span>Original (blue)&nbsp;&nbsp;"
                "<span style='display:inline-block;width:16px;height:16px;background:#ff4136;margin-right:6px;margin-left:12px;vertical-align:middle'></span>Mutant (red)"
                "</div>"
                "</div>")
    parts.append("<style>body{font-family:Arial,Helvetica,sans-serif;padding:18px}h1,h2{color:#2b2b2b}table.dataframe{border-collapse:collapse;width:100%}table.dataframe th, table.dataframe td{border:1px solid #ddd;padding:8px;text-align:left}tr:nth-child(even){background:#f9f9f9}"
                "#structure_viewer canvas{position:absolute;top:0;left:0}""</style>")
    script = (
        "<script>var origFile=" + json.dumps(orig_link) + ";var mutFile=" + json.dumps(mut_link) + ";"
        "function loadModels(){"
        "if(!origFile && !mutFile){document.getElementById('structure_viewer').innerHTML='<p style=\"padding:8px\">No structures available</p>';return;}"
        "var doLoad=function(){var v=$3Dmol.createViewer('structure_viewer',{defaultcolors:$3Dmol.rasmolElementColors});"
        "Promise.all([origFile?fetch(origFile).then(r=>r.text()):Promise.resolve(null), mutFile?fetch(mutFile).then(r=>r.text()):Promise.resolve(null)])"
        ".then(([origTxt,mutTxt])=>{var idx=0; if(origTxt){v.addModel(origTxt,'pdb'); v.setStyle({model:idx},{cartoon:{color:'blue'}}); idx++;} if(mutTxt){v.addModel(mutTxt,'pdb'); v.setStyle({model:idx},{cartoon:{color:'red',opacity:0.8}}); idx++;} v.zoomTo(); v.render();}).catch(e=>{document.getElementById('structure_viewer').innerHTML='<p style=\"padding:8px\">Failed to load structures</p>';}); };"
        "if(typeof $3Dmol==='undefined'){var s=document.createElement('script');s.src='https://3dmol.csb.pitt.edu/build/3Dmol-min.js';s.onload=doLoad;document.head.appendChild(s);}else{doLoad();}"
        "}"
        "loadModels();"
        "</script>"
    )
    parts.append(script)

    def _add_search_section(title, df, downloads):
        parts.append(f"<h2>{title}</h2>")
        if df is None or (hasattr(df, 'empty') and df.empty):
            parts.append("<p>No results.</p>")
            return
        parts.append(_df_to_html(df))
        # Do not add links to downloaded/top-k structures here; tables are sufficient

    _add_search_section("Custom DB Foldseek Results", custom_results_df, custom_downloads)
    # If custom structure ids were provided, list them (no download links)
    if custom_structure_ids:
        parts.append("<h3>Custom structure IDs</h3>")
        parts.append("<ul>")
        for sid in custom_structure_ids:
            parts.append(f"<li>{html.escape(str(sid))}</li>")
        parts.append("</ul>")
    _add_search_section("PDB Foldseek Results", pdb_results_df, pdb_downloads)
    _add_search_section("AlphaFold Foldseek Results", af_results_df, af_downloads)

    parts.append("</body></html>")

    html_content = "\n".join(parts)

    out_path = os.path.join(report_dir, "report.html")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)

    return out_path
