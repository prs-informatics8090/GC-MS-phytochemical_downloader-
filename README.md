# GC-MS-phytochemical_downloader-
This pipeline parses GC-MS PDF reports to extract compound names, retention times, CAS numbers, and molecular formulas. It then resolves each compound to its corresponding PubChem Compound ID (CID) using name-based searches with CAS number fallback. The pipeline retrieves the canonical SMILES representation and downloads 3D SDF structures
# GC-MS → Phytochemical_downloader

An pipeline that takes a GC-MS (Gas Chromatography–Mass Spectrometry) report in PDF form, extracts every identified compound, and automatically enriches each one with structural data from PubChem.

Given a standard GC-MS PDF report, the script:

1. Parses the PDF and extracts each compound's retention time (RT), name, CAS number, and molecular formula.
2. Resolves each compound to PubChem — first by name, then by CAS number as a fallback — to get its PubChem CID.
3. Fetches the canonical SMILES string for each resolved compound.
4. Downloads the 3D structure file (SDF) for each compound (falls back to 2D if no 3D conformer exists).
5. Writes a master `summary.csv` combining the GC-MS data with the PubChem identifiers — ready for downstream filtering, docking prioritization, or reporting.

Why

GC-MS reports typically list dozens to hundreds of compounds by name/CAS with no structural information attached. Manually looking each one up on PubChem is slow and error-prone at scale (this pipeline was built to process a GC-MS report of 260+ compounds). This script automates that lookup chain.

Installation

```bash
git clone https://github.com/<your-username>/gcms-pubchem-pipeline.git
cd gcms-pubchem-pipeline
pip install -r requirements.txt
```
Requires Python 3.8+.

Usage

```bash
python phytochemical_praser.py path/to/report.pdf
```
Optional arguments:

```bash
python phytochemical_praser.py path/to/report.pdf -o results_folder -d 0.3
```

| Flag | Description | Default |
|------|-------------|---------|
| `-o`, `--output` | Output folder for SDF files and CSVs | `results` |
| `-d`, `--delay`  | Delay (seconds) between PubChem API calls | `0.3` |

 Output structure

```
results/
├── sdf_files/
│   ├── Compound_A_CID1234_3D.sdf
│   └── ...
└── summary.csv   # master table: GC-MS data + PubChem identifiers
```

`summary.csv` columns: `RT, Compound Name, CAS, Formula, CID, SMILES, SDF Status`

Notes and limitations

- PDF parsing is format-sensitive. `extract_compounds()` is tuned for GC-MS reports that list `RT  Compound Name  CAS-Number  Formula` per line (typical NIST-library-style output). Reports in a different layout may need the regex/parsing logic adjusted.
- PubChem name lookup can miss synonyms. Compounds not found by name are retried by CAS number; those still unresolved are recorded with an empty CID/SMILES in the summary rather than dropped.
- Not affiliated with or endorsed by PubChem/NCBI — please respect their terms of use and rate limits.

REQUIREMENTS

- `pdfplumber` — PDF text extraction
- `requests` — HTTP calls to PubChem

License

MIT — see [LICENSE](LICENSE).

