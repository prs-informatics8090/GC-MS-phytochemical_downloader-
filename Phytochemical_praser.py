#!/usr/bin/env python3
"""
GC-MS REPORT → Phytochemical_Downloader  
=================================
1. Extracts compound names & CAS numbers from a GC-MS PDF report.
2. Resolves each compound to a PubChem CID.
3. Fetches the canonical SMILES from PubChem.
4. Downloads the 3D SDF structure file.
5. Writes a master summary CSV combining everything.

Usage
-----
    python phytochemical_praser.py  path/to/report.pdf
    python phytochemical_praser.py  path/to/report.pdf  -o results_folder  -d 0.3

Requirements
------------
    pip install pdfplumber requests
"""

import re
import time
import sys
import logging
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Missing: pip install requests")
try:
    import pdfplumber
except ImportError:
    sys.exit("Missing: pip install pdfplumber")

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────
PUBCHEM_BASE  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 PubChemScraper/2.0 "
        "(research automation; pdfplumber+requests)"
    )
}


# ══════════════════════════════════════════════════════════════════════════════
#  PubChem helpers
# ══════════════════════════════════════════════════════════════════════════════

def _get(url, retries=3, backoff=2.0):
    """Robust GET with retry on 429 / 5xx."""
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, headers=HEADERS, timeout=20)
            if r.status_code == 429:
                time.sleep(backoff * attempt)
                continue
            return r
        except requests.RequestException as exc:
            log.debug("GET error attempt %d: %s", attempt, exc)
            time.sleep(backoff)
    return None


def get_cid_by_name(name):
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    r = _get(url)
    if r and r.status_code == 200:
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    return None


def get_cid_by_cas(cas):
    url = f"{PUBCHEM_BASE}/compound/name/{cas}/cids/JSON"
    r = _get(url)
    if r and r.status_code == 200:
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    return None


def get_smiles(cid):
    """Return the canonical SMILES for a CID, or empty string."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    r = _get(url)
    if r and r.status_code == 200:
        props = r.json().get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("CanonicalSMILES", "")
    return ""


def download_sdf_3d(cid, out_path):
    """Download 3D SDF; fall back to 2D if unavailable."""
    for rtype in ("3d", "2d"):
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type={rtype}"
        r = _get(url)
        if r and r.status_code == 200:
            out_path.write_bytes(r.content)
            return True
        if r and r.status_code == 404 and rtype == "3d":
            log.warning("  No 3D conformer for CID %s – trying 2D", cid)
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  PDF parsing
# ══════════════════════════════════════════════════════════════════════════════
CAS_RE = re.compile(r"\b\d{2,7}-\d{2}-\d\b")


def extract_compounds(pdf_path):
    """
    Parse the GC-MS PDF. Returns list of:
        {"rt": float|None, "name": str, "cas": str, "formula": str}
    """
    compounds = []
    seen = set()

    with pdfplumber.open(pdf_path) as pdf:
        lines = []
        for page in pdf.pages:
            text = page.extract_text()
            if text:
                lines.extend(text.splitlines())

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        cas_match = CAS_RE.search(line)
        if cas_match:
            cas = cas_match.group()
            before = line[: cas_match.start()].strip()
            tokens = before.split()
            rt_val, name_toks = None, []
            for tok in tokens:
                if rt_val is None:
                    try:
                        rt_val = float(tok); continue
                    except ValueError:
                        pass
                name_toks.append(tok)
            compound_name = " ".join(name_toks).strip()

            if not compound_name:
                j, parts = i - 1, []
                while j >= 0 and not CAS_RE.search(lines[j]):
                    prev = lines[j].strip()
                    if re.match(r"^(RT\s+Compound|Unknown Analysis|Batch|Sample)", prev):
                        break
                    if prev:
                        parts.insert(0, prev)
                    j -= 1
                compound_name = " ".join(parts).strip()
                toks = compound_name.split()
                try:
                    float(toks[0]); compound_name = " ".join(toks[1:])
                except (ValueError, IndexError):
                    pass

            compound_name = re.sub(r"\s+", " ", compound_name).strip()
            after = line[cas_match.end():].strip()
            fm = re.match(r"(C\d+H\d+\S*)", after)
            formula = fm.group(1) if fm else ""

            if (compound_name and len(compound_name) > 3
                    and compound_name not in seen
                    and not compound_name.lower().startswith("rt ")):
                seen.add(compound_name)
                compounds.append({
                    "rt": rt_val, "name": compound_name,
                    "cas": cas, "formula": formula,
                })
        i += 1
    return compounds


# ══════════════════════════════════════════════════════════════════════════════
#  Main pipeline
# ══════════════════════════════════════════════════════════════════════════════

def scrape(pdf_path, output_dir="results", delay=0.3):
    """
    Full pipeline:
      1. PDF  → compound names + CAS
      2. PubChem → CID + SMILES + 3D SDF
      3. Master summary CSV
    """
    out = Path(output_dir)
    sdf_dir = out / "sdf_files"
    sdf_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Extract compounds ──────────────────────────────────────────
    log.info("Reading PDF: %s", pdf_path)
    compounds = extract_compounds(pdf_path)
    if not compounds:
        log.error("No compounds found. Check PDF format."); return
    log.info("Extracted %d unique compounds from PDF.", len(compounds))

    results = []

    # ── Step 2: PubChem CID + SMILES + SDF ────────────────────────────────
    for idx, comp in enumerate(compounds, 1):
        name, cas = comp["name"], comp["cas"]
        log.info("[%d/%d] '%s'", idx, len(compounds), name)

        cid = get_cid_by_name(name);  time.sleep(delay)
        if cid is None and cas:
            log.info("  Name lookup failed – trying CAS %s", cas)
            cid = get_cid_by_cas(cas); time.sleep(delay)

        if cid is None:
            log.warning("  ✗ Not found on PubChem")
            results.append({**comp, "cid": "", "smiles": "",
                             "sdf_status": "not_found", "sdf_file": ""})
            continue

        log.info("  CID: %s", cid)

        # SMILES
        smiles = get_smiles(cid);  time.sleep(delay)
        if smiles:
            log.info("  SMILES: %s", smiles[:60] + ("…" if len(smiles) > 60 else ""))
        else:
            log.warning("  No SMILES available")

        # SDF
        safe = re.sub(r"[^\w\-]", "_", name)[:80]
        sdf_path = sdf_dir / f"{safe}_CID{cid}_3D.sdf"
        if sdf_path.exists():
            log.info("  ✓ SDF already exists")
            sdf_status = "already_exists"
        else:
            ok = download_sdf_3d(cid, sdf_path); time.sleep(delay)
            sdf_status = "downloaded" if ok else "no_conformer"
            log.info("  %s SDF: %s", "✓" if ok else "✗", sdf_path.name)

        results.append({
            **comp,
            "cid":        cid,
            "smiles":     smiles,
            "sdf_status": sdf_status,
            "sdf_file":   str(sdf_path) if sdf_path.exists() else "",
        })

    # ── Step 3: Master summary CSV ────────────────────────────────────────
    summary_csv = out / "summary.csv"
    with open(summary_csv, "w", encoding="utf-8") as f:
        header = "RT,Compound Name,CAS,Formula,CID,SMILES,SDF Status\n"
        f.write(header)
        for r in results:
            name_esc  = r["name"].replace(",", ";")
            smiles_esc = r["smiles"].replace(",", ";") if r["smiles"] else ""
            f.write(
                f"{r.get('rt','')},{name_esc},{r['cas']},{r['formula']},"
                f"{r.get('cid','')},{smiles_esc},{r['sdf_status']}\n"
            )
    log.info("Master summary saved: %s", summary_csv)

    # ── Final summary ─────────────────────────────────────────────────────
    total      = len(results)
    found      = sum(1 for r in results if r["cid"])
    with_smi   = sum(1 for r in results if r["smiles"])
    sdf_ok     = sum(1 for r in results if r["sdf_status"] in ("downloaded", "already_exists"))

    print("\n" + "=" * 64)
    print("  PIPELINE SUMMARY")
    print("=" * 64)
    print(f"  Compounds in PDF             : {total}")
    print(f"  Found on PubChem (CID)       : {found}")
    print(f"  SMILES retrieved             : {with_smi}")
    print(f"  3D SDF files saved           : {sdf_ok}")
    print(f"  Output directory             : {out.resolve()}")
    print(f"  Master summary               : {summary_csv.name}")
    print("=" * 64 + "\n")


# ══════════════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "GC-MS PDF → PubChem (CID + SMILES + 3D SDF)"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pdf",  help="Path to the GC-MS PDF report")
    parser.add_argument("-o", "--output",  default="results",
                        help="Output folder for SDF files and CSVs")
    parser.add_argument("-d", "--delay",   type=float, default=0.3,
                        help="Delay (s) between PubChem API calls")
    args = parser.parse_args()

    if not Path(args.pdf).is_file():
        sys.exit(f"Error: file not found: '{args.pdf}'")

    scrape(args.pdf, output_dir=args.output, delay=args.delay)

## Developed by Paras Dhiman
