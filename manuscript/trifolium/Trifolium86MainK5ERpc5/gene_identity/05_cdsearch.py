#!/usr/bin/env python3
"""05_cdsearch.py -- conserved domains of the 353 Trifolium proteins with the NCBI Batch CD-Search web service.

Run with `python3 -I` (stdlib only).

Usage:
  python3 -I 05_cdsearch.py <data_dir> [poll_seconds=10] [max_wait_seconds=3600]
  python3 -I 05_cdsearch.py <data_dir> --reparse    # re-parse the saved service output, no new search

One batch is submitted (POST to bwrpsb.cgi; the service accepts up to 1,000 proteins per request) with the
documented defaults written out explicitly:
  db=cdd, smode=auto, useid1=true, compbasedadj=1 (composition-based statistics), filter=false (no
  low-complexity filter), evalue=0.01, maxhit=500, dmode=rep (concise: the top-scoring hit per region of the
  query), tdata=hits.
(evalue and maxhit only act in live mode; novel sequences are searched live under smode=auto.)
The job is polled every <poll_seconds> s until it reports status 0 (success); any other final status, or the
time cap, stops the script with an error (no other domain method is substituted). Warning lines that the
service prints next to status 0 are recorded in cdsearch_meta.tsv.
The service output states no CDD version. For the record only, the script also saves NCBI's FTP release file
(pub/mmdb/cdd/cdd.info); that file describes the public CDD release, not necessarily the service's database.

Reads  <data_dir>/work/queries.faa
Writes <data_dir>/cdsearch/{submit_response.txt, cdsearch_raw.txt, cdsearch_hits.tsv, cdsearch_meta.tsv,
       ncbi_ftp_cdd.info}
"""
import csv
import datetime
import os
import re
import sys
import time
import urllib.parse
import urllib.request

URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
CDD_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.info"
PARAMS = {"db": "cdd", "smode": "auto", "useid1": "true", "compbasedadj": "1", "filter": "false",
          "evalue": "0.01", "maxhit": "500", "dmode": "rep", "tdata": "hits"}
STATUS = {"0": "success", "1": "invalid search ID", "2": "no effective input", "3": "search ongoing",
          "4": "queue manager service error", "5": "data corrupted or no longer available"}


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def http(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"User-Agent": "trifolium-gene-identity/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read().decode("utf-8", "replace")


def status_of(text):
    m = re.search(r"^#status\t(\d+)", text, re.M)
    return m.group(1) if m else None


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    data = sys.argv[1]
    reparse = len(sys.argv) > 2 and sys.argv[2] == "--reparse"
    poll = int(sys.argv[2]) if len(sys.argv) > 2 and not reparse else 10
    max_wait = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
    out = os.path.join(data, "cdsearch")      # service output only (untrusted text, parsed as data)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(data, "work", "queries.faa")) as fh:
        fasta = fh.read()
    n_q = fasta.count(">")
    if reparse:
        old = dict(ln.rstrip("\n").split("\t", 1) for ln in open(os.path.join(out, "cdsearch_meta.tsv")))
        with open(os.path.join(out, "cdsearch_raw.txt")) as fh:
            text = fh.read()
        parse_and_write(out, text, old["cdsid"], old["submitted_utc"], old["finished_utc"], n_q)
        return

    submitted = now()
    body = urllib.parse.urlencode(dict(PARAMS, queries=fasta)).encode()
    resp = http(URL, body)
    with open(os.path.join(out, "submit_response.txt"), "w") as fh:
        fh.write(resp)
    m = re.search(r"^#cdsid\t(\S+)", resp, re.M)
    st = status_of(resp)
    if not m or st not in ("0", "3"):
        sys.exit(f"[05] ERROR: submission failed (status {st}: {STATUS.get(st, '?')}); see submit_response.txt")
    cdsid = m.group(1)
    print(f"[05] submitted {n_q} proteins, cdsid {cdsid}")

    waited, text, http_errors = 0, resp, 0
    query = urllib.parse.urlencode({"cdsid": cdsid, "tdata": PARAMS["tdata"], "dmode": PARAMS["dmode"]})
    while True:
        time.sleep(poll)
        waited += poll
        try:
            text = http(f"{URL}?{query}")
            http_errors = 0
        except Exception as exc:  # transient network/HTTP error: retry a few times, then stop
            http_errors += 1
            print(f"[05] poll error {http_errors}: {exc}")
            if http_errors >= 5:
                sys.exit(f"[05] ERROR: 5 consecutive polling errors for {cdsid}")
            continue
        st = status_of(text)
        if st == "0":
            break
        if st != "3":
            sys.exit(f"[05] ERROR: CD-Search status {st} ({STATUS.get(st, 'unknown')}) for {cdsid}")
        if waited >= max_wait:
            sys.exit(f"[05] ERROR: CD-Search did not finish within {max_wait} s ({cdsid})")
    finished = now()
    with open(os.path.join(out, "cdsearch_raw.txt"), "w") as fh:
        fh.write(text)

    parse_and_write(out, text, cdsid, submitted, finished, n_q)


def parse_and_write(out, text, cdsid, submitted, finished, n_q):
    # parse the hits table (header line starts with "Query\t")
    lines = text.splitlines()
    try:
        h = next(i for i, ln in enumerate(lines) if ln.startswith("Query\t"))
    except StopIteration:
        sys.exit("[05] ERROR: no hits table in the CD-Search output")
    header = lines[h].split("\t")
    rows = []
    for ln in lines[h + 1:]:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split("\t")
        rec = dict(zip(header, f))
        qm = re.search(r">(\S+)", rec["Query"])
        rec["query_id"] = qm.group(1) if qm else rec["Query"]
        rows.append(rec)
    cols = ["query_id", "Hit type", "PSSM-ID", "From", "To", "E-Value", "Bitscore", "Accession", "Short name",
            "Incomplete", "Superfamily"]
    with open(os.path.join(out, "cdsearch_hits.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([c.lower().replace(" ", "_").replace("-", "_") for c in cols])
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])

    # version information: only if the service output states one
    ver = [ln for ln in lines if re.search(r"(?i)\bcdd\b.*\bv?\d+\.\d+|version|data ?source", ln)]
    ver_str = " | ".join(v.strip("#").strip() for v in ver) if ver else "not reported by service"
    warnings = [ln.strip("#").strip() for ln in lines
                if re.match(r"^#status\t(?!\d)", ln) and ln.strip() != "#status\tsuccess"]
    try:
        cdd_info = http(CDD_INFO_URL)
        with open(os.path.join(out, "ncbi_ftp_cdd.info"), "w") as fh:
            fh.write(cdd_info)
        ftp_release = cdd_info.splitlines()[0].strip()
    except Exception as exc:  # recorded, not fatal: it is context, not the service output
        ftp_release = f"could not fetch {CDD_INFO_URL}: {exc}"
    n_with = len({r["query_id"] for r in rows})
    n_spec = len({r["query_id"] for r in rows if r.get("Hit type", "").lower() == "specific"})
    meta = [("service", URL), ("cdsid", cdsid), ("submitted_utc", submitted), ("finished_utc", finished),
            ("n_queries", n_q), ("n_queries_with_any_hit", n_with), ("n_queries_with_specific_hit", n_spec),
            ("n_hit_rows", len(rows)), ("cdd_version_in_service_output", ver_str),
            ("service_warnings", " | ".join(warnings) if warnings else "none"),
            ("ncbi_ftp_cdd_info_first_line", ftp_release)]
    meta += [(f"param_{k}", v) for k, v in PARAMS.items()]
    with open(os.path.join(out, "cdsearch_meta.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "value"])
        w.writerows(meta)
    for k, v in meta:
        print(f"[05] {k}\t{v}")


if __name__ == "__main__":
    main()
