#!/usr/bin/env python3
"""
Pull Bangladesh-origin bacterial data from NCBI, using BioSample as the hub.

Stages (run in order, or use --stage all):
  biosample : search + download BioSample metadata, verify country from the
              geo_loc_name attribute, write clean table
  sra       : fetch SRA run tables for the verified BioSamples
  assembly  : fetch genome assembly summaries for the verified BioSamples

Setup:
  pip install biopython pandas
  export NCBI_EMAIL="you@example.com"
  export NCBI_API_KEY="your_key"      # optional but recommended (10 req/s)

Usage:
  python ncbi_bangladesh.py --stage biosample --dry-run   # counts only
  python ncbi_bangladesh.py --stage all
"""
import argparse
import datetime
import io
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from Bio import Entrez, __version__ as BIOPYTHON_VERSION


# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
ORGANISMS = {
    "Escherichia coli": 562,
    "Klebsiella pneumoniae": 573,
    "Salmonella enterica": 28901,        # includes all serovars (and Typhi)
    "Salmonella Typhi": 90370,           # subset of 28901; deduplicated later
    "Acinetobacter baumannii": 470,
    "Enterococcus faecium": 1352,
    "Staphylococcus aureus": 1280,
    "Pseudomonas aeruginosa": 287,
    "Vibrio cholerae": 666,
    "Mammaliicoccus sciuri": 1296,
}

# BioSample attributes to pull into their own columns (harmonized names).
# Everything else is kept in the all_attributes JSON column.
KEEP_ATTRS = [
    "collection_date", "host", "isolation_source", "strain", "isolate",
    "host_disease", "sample_type", "serovar",
]

BATCH = 200
ROOT = Path("ncbi_bangladesh")
RAW = ROOT / "raw"
CLEAN = ROOT / "clean"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def setup():
    email = os.environ.get("NCBI_EMAIL","sajirhtml@gmail.com")
    if not email:
        sys.exit("Please set NCBI_EMAIL (NCBI asks for a contact address).")
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY","6dd72a8488426759cd7a22fdad3a39047508")
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 5
    for d in (RAW / "biosample", RAW / "sra", RAW / "assembly", CLEAN):
        d.mkdir(parents=True, exist_ok=True)


def search_history(db, term):
    """Run esearch with history; return (count, WebEnv, QueryKey)."""
    with Entrez.esearch(db=db, term=term, usehistory="y", retmax=0) as h:
        r = Entrez.read(h)
    return int(r["Count"]), r["WebEnv"], r["QueryKey"]


def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def split_location(raw):
    """'Bangladesh: Dhaka' / 'Dhaka, Bangladesh' / 'BD' -> (country, region)."""
    raw = (raw or "").strip()
    if ":" in raw:
        country, region = [x.strip() for x in raw.split(":", 1)]
    else:
        country, region = raw, ""
    is_bd = "bangladesh" in country.lower() or country.upper() in {"BD", "BGD"}
    if is_bd:
        if not region and "," in country:
            region = country.split(",")[0].strip()
        country = "Bangladesh"
    return country, region


def parse_biosamples(xml_bytes):
    """Yield one flat dict per <BioSample> element."""
    root = ET.fromstring(xml_bytes)
    for bs in root.iter("BioSample"):
        attrs = {}
        for a in bs.findall("Attributes/Attribute"):
            key = a.get("harmonized_name") or a.get("attribute_name")
            attrs.setdefault(key, (a.text or "").strip())

        org = bs.find("Description/Organism")
        loc_raw = attrs.get("geo_loc_name", "")
        country, region = split_location(loc_raw)

        row = {
            "biosample": bs.get("accession"),
            "sra_sample": next(
                (i.text for i in bs.findall("Ids/Id") if i.get("db") == "SRA"), ""),
            "bioproject": ";".join(
                l.get("label") or "" for l in bs.findall("Links/Link")
                if l.get("target") == "bioproject"),
            "organism": org.get("taxonomy_name") if org is not None else "",
            "taxid": org.get("taxonomy_id") if org is not None else "",
            "title": bs.findtext("Description/Title", default=""),
            "owner": bs.findtext("Owner/Name", default=""),
            "submission_date": bs.get("submission_date", ""),
            "publication_date": bs.get("publication_date", ""),
            "last_update": bs.get("last_update", ""),
            "geo_loc_name_raw": loc_raw,
            "country": country,
            "region": region,
        }
        for k in KEEP_ATTRS:
            row[k] = attrs.get(k, "")
        row["all_attributes"] = json.dumps(attrs, ensure_ascii=False)
        yield row


def load_biosample_ids():
    path = CLEAN / "biosample_bangladesh.tsv"
    if not path.exists():
        sys.exit(f"{path} not found. Run --stage biosample first.")
    return pd.read_csv(path, sep="\t")["biosample"].dropna().unique().tolist()


def write_log(stage, info):
    path = ROOT / "run_log.json"
    log = json.loads(path.read_text()) if path.exists() else {}
    log[stage] = {
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "biopython": BIOPYTHON_VERSION,
        **info,
    }
    path.write_text(json.dumps(log, indent=2))


# ----------------------------------------------------------------------------
# Stage 1: BioSample
# ----------------------------------------------------------------------------
def stage_biosample(dry_run=False):
    rows, summary = [], []
    for name, taxid in ORGANISMS.items():
        term = f"txid{taxid}[Organism:exp] AND Bangladesh[All Fields]"
        count, web, qk = search_history("biosample", term)
        print(f"{name:26} txid{taxid:<6} {count:>7} raw hits")
        summary.append({"query_name": name, "taxid": taxid,
                        "query": term, "raw_hits": count})
        if dry_run or count == 0:
            continue
        for start in range(0, count, BATCH):
            with Entrez.efetch(db="biosample", query_key=qk, WebEnv=web,
                               retstart=start, retmax=BATCH, retmode="xml") as h:
                xml = h.read()
            (RAW / "biosample" / f"{taxid}_{start:06d}.xml").write_bytes(xml)
            for r in parse_biosamples(xml):
                r["query_name"], r["query_taxid"] = name, taxid
                rows.append(r)
            print(f"   fetched {min(start + BATCH, count)}/{count}")

    if dry_run:
        return
    if not rows:
        sys.exit("No records fetched.")

    df = pd.DataFrame(rows)
    # Audit file: every hit, including ones that only mention Bangladesh in text
    df.to_csv(CLEAN / "biosample_all_hits_unfiltered.tsv", sep="\t", index=False)

    # Keep only records whose geo_loc_name really says Bangladesh
    bd = df[df["country"] == "Bangladesh"].copy()

    # Deduplicate (Typhi appears under both 90370 and 28901)
    matched = (bd.groupby("biosample")["query_name"]
                 .agg(lambda s: "; ".join(sorted(set(s))))
                 .rename("matched_queries"))
    bd = (bd.drop_duplicates("biosample")
            .drop(columns=["query_name", "query_taxid"])
            .merge(matched, on="biosample"))
    bd.to_csv(CLEAN / "biosample_bangladesh.tsv", sep="\t", index=False)

    s = pd.DataFrame(summary)
    s["verified_bangladesh"] = s["query_name"].map(
        df[df["country"] == "Bangladesh"].groupby("query_name")["biosample"].nunique()
    ).fillna(0).astype(int)
    s.to_csv(CLEAN / "summary_counts.tsv", sep="\t", index=False)

    print("\nVerified, deduplicated Bangladesh BioSamples:", len(bd))
    print(s[["query_name", "raw_hits", "verified_bangladesh"]].to_string(index=False))
    print("\nTop regions:")
    print(bd["region"].replace("", "(none)").value_counts().head(15).to_string())
    write_log("biosample", {"queries": summary, "unique_verified": len(bd)})


# ----------------------------------------------------------------------------
# Stage 2: SRA run tables
# ----------------------------------------------------------------------------
def stage_sra():
    ids = load_biosample_ids()
    frames = []
    for i, ch in enumerate(chunks(ids, 100)):
        term = " OR ".join(f"{x}[BioSample]" for x in ch)
        count, web, qk = search_history("sra", term)
        print(f"SRA chunk {i + 1}: {count} experiments")
        if count == 0:
            continue
        with Entrez.efetch(db="sra", query_key=qk, WebEnv=web, rettype="runinfo",
                           retmode="text", retmax=10000) as h:
            txt = h.read()
        if isinstance(txt, bytes):
            txt = txt.decode("utf-8", errors="replace")
        (RAW / "sra" / f"runinfo_{i:04d}.csv").write_text(txt)
        df = pd.read_csv(io.StringIO(txt), low_memory=False)
        df = df[df["Run"].notna() & (df["Run"] != "Run")]
        frames.append(df)

    if not frames:
        print("No SRA runs found for these BioSamples.")
        return
    out = pd.concat(frames, ignore_index=True).drop_duplicates("Run")
    out = out[out["BioSample"].isin(ids)]  # keep only our verified samples
    out.to_csv(CLEAN / "sra_runs_bangladesh.tsv", sep="\t", index=False)
    print(f"SRA runs: {len(out)} (from {out['BioSample'].nunique()} BioSamples)")
    write_log("sra", {"runs": len(out)})


# ----------------------------------------------------------------------------
# Stage 3: Assemblies
# ----------------------------------------------------------------------------
ASM_FIELDS = [
    "AssemblyAccession", "AssemblyName", "Organism", "SpeciesName", "Taxid",
    "AssemblyStatus", "AssemblyType", "BioSampleAccn", "SubmitterOrganization",
    "SeqReleaseDate", "ContigN50", "Coverage", "FtpPath_GenBank", "FtpPath_RefSeq",
]


def stage_assembly():
    ids = load_biosample_ids()
    rows = []
    for i, ch in enumerate(chunks(ids, 100)):
        term = " OR ".join(f"{x}[BioSample]" for x in ch)
        with Entrez.esearch(db="assembly", term=term, retmax=10000) as h:
            uids = Entrez.read(h)["IdList"]
        print(f"Assembly chunk {i + 1}: {len(uids)} assemblies")
        if not uids:
            continue
        with Entrez.esummary(db="assembly", id=",".join(uids), report="full") as h:
            res = Entrez.read(h, validate=False)
        for d in res["DocumentSummarySet"]["DocumentSummary"]:
            rows.append({k: str(d.get(k, "")) for k in ASM_FIELDS})

    if not rows:
        print("No assemblies found for these BioSamples.")
        return
    df = pd.DataFrame(rows).drop_duplicates("AssemblyAccession")
    df.to_csv(CLEAN / "assemblies_bangladesh.tsv", sep="\t", index=False)
    (CLEAN / "assembly_accessions.txt").write_text(
        "\n".join(df["AssemblyAccession"]) + "\n")
    print(f"Assemblies: {len(df)}")
    print("To download genomes (NCBI Datasets CLI):")
    print("  datasets download genome accession --inputfile "
          f"{CLEAN / 'assembly_accessions.txt'} --include genome,gff3 --filename genomes.zip")
    write_log("assembly", {"assemblies": len(df)})


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=["biosample", "sra", "assembly", "all"],
                   default="biosample")
    p.add_argument("--dry-run", action="store_true",
                   help="biosample stage only: print counts, download nothing")
    args = p.parse_args()

    setup()
    if args.stage in ("biosample", "all"):
        stage_biosample(dry_run=args.dry_run)
    if args.stage in ("sra", "all") and not args.dry_run:
        stage_sra()
    if args.stage in ("assembly", "all") and not args.dry_run:
        stage_assembly()



# pip install biopython pandas
# setx NCBI_EMAIL "sajirhtml@gmail.com"
# setx NCBI_API_KEY "6dd72a8488426759cd7a22fdad3a39047508"
# python ncbi_bangladesh.py --stage biosample --dry-run   # counts only, downloads nothing
# python ncbi_bangladesh.py --stage all                   # full pipeline
