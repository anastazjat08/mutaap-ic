import requests
import time
import json
from typing import Optional
import pandas as pd

IPRSCAN_BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
POLL_INTERVAL = 15  # seconds
MAX_WAIT = 600      # 10 minutes timeout

GO_ASPECT_MAP = {
    "molecular_function": "mf",
    "biological_process": "bp",
    "cellular_component": "cc",
}


def submit_iprscan(sequence: str, email: str, label: str = "") -> Optional[str]:
    """Submit a sequence to InterProScan and return the job ID."""

    # Clean and validate sequence
    clean_seq = sequence.replace(" ", "").replace("\n", "").upper()
    if not all(c in "ACDEFGHIKLMNPQRSTVWXY*" for c in clean_seq):
        raise ValueError(f"Invalid characters in sequence for '{label}'.")

    # API requires FASTA format
    fasta_seq = f">{label or 'sequence'}\n{clean_seq}"

    payload = {
        "email": email,
        "title": label or "sequence",
        "sequence": fasta_seq,
        "goterms":  "true",
        "pathways": "false",
    }

    print(f"[InterProScan] Submitting {label or 'sequence'}...")
    r = requests.post(f"{IPRSCAN_BASE}/run", data=payload)

    if r.status_code != 200:
        print(f"[InterProScan] Server response: {r.text}")  # print actual error
        r.raise_for_status()

    job_id = r.text.strip()
    print(f"[InterProScan] Job ID: {job_id}")
    return job_id


def poll_iprscan(job_id: str) -> bool:
    """Poll until job is finished. Returns True on success, False on failure/timeout."""
    elapsed = 0
    while elapsed < MAX_WAIT:
        status = requests.get(f"{IPRSCAN_BASE}/status/{job_id}").text.strip()
        print(f"[InterProScan] Status: {status} ({elapsed}s elapsed)")
        if status == "FINISHED":
            return True
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            print(f"[InterProScan] Job failed with status: {status}")
            return False
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
    print("[InterProScan] Timed out waiting for results.")
    return False


def fetch_go_terms(job_id: str) -> dict[str, list[dict]]:
    r = requests.get(f"{IPRSCAN_BASE}/result/{job_id}/json")
    r.raise_for_status()
    data = r.json()

    go_terms = {"mf": [], "bp": [], "cc": []}
    seen = set()

    for match in data.get("results", [{}])[0].get("matches", []):
        # GO terms are under signature.entry, not entry directly
        entry = match.get("signature", {}).get("entry") or {}
        for xref in entry.get("goXRefs") or []:
            go_id   = xref.get("id", "")
            go_name = xref.get("name", "")
            aspect  = xref.get("category", "").lower()
            ont_key = GO_ASPECT_MAP.get(aspect)

            if ont_key and go_id and go_id not in seen:
                seen.add(go_id)
                go_terms[ont_key].append({"GO_term": go_id, "description": go_name})

    print(f"[InterProScan] GO terms found: { {k: len(v) for k, v in go_terms.items()} }")
    return go_terms


def predict_function(sequence: str, email: str, label: str = "") -> Optional[dict[str, list[dict]]]:
    """
    Full pipeline: submit → poll → fetch GO terms.
    Returns dict of GO terms per ontology, or None on failure.
    """
    try:
        job_id = submit_iprscan(sequence, email, label)
        success = poll_iprscan(job_id)
        if not success:
            return None
        return fetch_go_terms(job_id)
    except Exception as e:
        print(f"[InterProScan] ERROR for {label}: {e}")
        return None


def compare_function(orig_terms: dict, mut_terms: dict) -> dict[str, pd.DataFrame]:
    results = {}

    for ont in ["mf", "bp", "cc"]:
        orig_set = {t["GO_term"]: t["description"] for t in (orig_terms.get(ont) or [])}
        mut_set  = {t["GO_term"]: t["description"] for t in (mut_terms.get(ont) or [])}
        all_terms = set(orig_set) | set(mut_set)

        rows = []
        for term in all_terms:
            in_orig = term in orig_set
            in_mut  = term in mut_set
            if in_orig and not in_mut:
                status = "lost"
            elif in_mut and not in_orig:
                status = "gained"
            else:
                status = "conserved"

            rows.append({
                "GO_term":     term,
                "description": orig_set.get(term) or mut_set.get(term, ""),
                "in_orig":     in_orig,
                "in_mut":      in_mut,
                "status":      status,
            })

        if not rows:
            results[ont] = pd.DataFrame(columns=["GO_term", "description", "in_orig", "in_mut", "status"])
            continue

        df = pd.DataFrame(rows)
        # Use Categorical ordering instead of key mapping
        df["status"] = pd.Categorical(df["status"], categories=["lost", "gained", "conserved"], ordered=True)
        df = df.sort_values("status").reset_index(drop=True)

        results[ont] = df

    return results