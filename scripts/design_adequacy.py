#!/usr/bin/env python3
"""design_adequacy.py — report what the SAMPLING DESIGN can support, before any GEA runs.

Every GEA/offset method in this pipeline is fitted on sample rows, but the
environment varies only between SITES. Three properties of the design therefore
bound everything downstream:

  1. Environmental degrees of freedom.  n_samples is not the unit of replication
     for a climate association; n_sites is. 58 samples at 11 sites carry 10
     environmental d.f., not 57 — the extra rows are replicates of an
     environmental value, not new observations of it.

  2. Predictor-block conditioning.  p predictors measured at s sites give a
     correlation matrix of rank at most min(p, s-1). With p >= s the block is
     rank-deficient BY CONSTRUCTION, and any method that inverts (or implicitly
     inverts) the environmental covariance — geometric offset's quadratic form,
     RDA's constrained ordination, varpart's adjusted R^2 — is then operating on
     directions that carry no information.

  3. Site balance.  Singleton sites and wildly unequal per-site n change what a
     site-level statistic means without changing any config value.

What is and is NOT new here (be precise, the pipeline already reports some of it):

  * NEW — environmental d.f., site-level residual d.f., site balance/singletons,
    samples-per-environmental-point, within-site variance share. Nothing in the
    pipeline computes any of these.

  * ALREADY REPORTED, but only downstream and only post-hoc — rank and condition
    number. `geometric_offset.R` (~line 376) writes `env_cov_rank_numeric` and
    `env_cov_condition_number` into geometric_offset_diagnostics.tsv, and
    `pregea_rda_setup.R` (~line 130) screens site-level max|r| per RDA rung.
    Both run after a model has been fitted, and both live in optional modes.
    The contribution of this script is TIMING (mode=climate, before any fit) and
    that it does not depend on which downstream method a project happens to run.

  * The numbers will NOT match geometric_offset.R's. That script eigen-decomposes
    the environmental COVARIANCE of the per-sample block; this one decomposes the
    site-level CORRELATION matrix (scale-free, one row per site, invariant
    predictors excluded). Different matrices, different floors, different rank
    thresholds. For the geometric offset's own numerical behaviour the covariance
    figure is authoritative — that is the matrix LEA actually inverts. This table
    is the DESIGN statement: whether the predictor block could ever have been
    well-conditioned at this number of sites, regardless of method.

Reads outputs the pipeline already produces; writes one tidy TSV.

Usage:
  python3 scripts/design_adequacy.py METADATA_TSV CLIMATE_SITE_TSV PREDICTORS OUT_TSV

  METADATA_TSV      Processing/tables/metadata.tsv (columns: site, sample, latitude, longitude, ...)
  CLIMATE_SITE_TSV  climate/tables/present/climate_present_site.tsv (columns: sample, bio_*)
  PREDICTORS        comma-separated predictor names, or 'all'
  OUT_TSV           output table (metric, value, flag, note)
"""
import math
import sys

# --------------------------------------------------------------------- linalg
def jacobi_eigenvalues(a, max_sweeps=100, tol=1e-12):
    """Eigenvalues of a symmetric matrix (list of lists), cyclic Jacobi. Descending."""
    n = len(a)
    m = [row[:] for row in a]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(m[i][j] ** 2 for i in range(n) for j in range(n) if i != j))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(m[p][q]) < 1e-300:
                    continue
                theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q])
                t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    mkp, mkq = m[k][p], m[k][q]
                    m[k][p] = c * mkp - s * mkq
                    m[k][q] = s * mkp + c * mkq
                for k in range(n):
                    mpk, mqk = m[p][k], m[q][k]
                    m[p][k] = c * mpk - s * mqk
                    m[q][k] = s * mpk + c * mqk
    return sorted((m[i][i] for i in range(n)), reverse=True)


def corr_matrix(cols):
    """Pearson correlation matrix of a list of equal-length numeric columns."""
    p = len(cols)
    n = len(cols[0])
    mu = [sum(c) / n for c in cols]
    sd = [math.sqrt(sum((x - mu[j]) ** 2 for x in cols[j]) / (n - 1)) if n > 1 else 0.0
          for j in range(p)]
    out = [[0.0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if sd[i] == 0 or sd[j] == 0:
                out[i][j] = 1.0 if i == j else 0.0
            else:
                cov = sum((cols[i][k] - mu[i]) * (cols[j][k] - mu[j]) for k in range(n)) / (n - 1)
                out[i][j] = cov / (sd[i] * sd[j])
    return out, sd


# ------------------------------------------------------------------------- io
def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    return header, rows


def main(argv):
    if len(argv) != 5:
        sys.stderr.write(__doc__)
        return 2
    meta_path, clim_path, predictors_arg, out_path = argv[1:]

    mhdr, mrows = read_tsv(meta_path)
    chdr, crows = read_tsv(clim_path)

    if "sample" not in mhdr or "site" not in mhdr:
        sys.stderr.write("FATAL: metadata must carry 'site' and 'sample' columns\n")
        return 1
    if "sample" not in chdr:
        sys.stderr.write("FATAL: climate site table must carry a 'sample' column\n")
        return 1

    m_i = {k: i for i, k in enumerate(mhdr)}
    site_of = {r[m_i["sample"]]: r[m_i["site"]] for r in mrows}

    if predictors_arg == "all":
        predictors = [c for c in chdr if c != "sample"]
    else:
        predictors = [p for p in predictors_arg.split(",") if p]
    missing = [p for p in predictors if p not in chdr]
    if missing:
        sys.stderr.write(f"FATAL: predictors absent from climate table: {','.join(missing)}\n")
        return 1

    c_i = {k: i for i, k in enumerate(chdr)}
    # Collapse to ONE row per site. Predictors are expected to be constant within
    # a site (every sample at a site queries the same raster cell), but that is
    # an assumption about how the climate table was built, not a guarantee, so
    # the site value is the mean and the within-site share of variance is
    # measured and reported rather than assumed to be zero.
    per_site_rows, site_n = {}, {}
    unmapped = 0
    unparseable = 0
    for r in crows:
        smp = r[c_i["sample"]]
        st = site_of.get(smp)
        if st is None:
            unmapped += 1
            continue
        try:
            vals = [float(r[c_i[p]]) for p in predictors]
        except (ValueError, IndexError):
            # Empty / non-numeric climate cell. R's as.numeric("") gives NA and
            # is silent; float("") raises. Count it — a row that contributed no
            # predictor values must NOT be counted towards that site's n, or the
            # per-site counts and samples_per_environmental_point are inflated by
            # rows that carry no environment at all.
            unparseable += 1
            continue
        site_n[st] = site_n.get(st, 0) + 1
        per_site_rows.setdefault(st, []).append(vals)

    site_vals = {st: [sum(row[j] for row in rows_) / len(rows_) for j in range(len(predictors))]
                 for st, rows_ in per_site_rows.items() if rows_}
    sites = sorted(site_vals)
    s = len(sites)
    p = len(predictors)
    n_samples = sum(site_n.get(k, 0) for k in sites)

    rows = []

    def add(metric, value, flag="", note=""):
        rows.append((metric, value, flag, note))

    add("n_samples", n_samples)
    add("n_sites", s)
    add("n_predictors", p)

    # n_samples above counts climate rows that PARSED. The metadata count is the
    # design as sampled. They differ exactly when the pipeline dropped rows
    # (climate_na_excluded.tsv) — i.e. when the design is at its messiest — so
    # report both rather than let the denominator quietly shrink.
    n_meta = len(mrows)
    n_meta_sites = len(set(site_of.values()))
    add("n_samples_in_metadata", n_meta)
    add("n_sites_in_metadata", n_meta_sites)
    if n_meta != n_samples or n_meta_sites != s:
        add("samples_missing_from_climate", n_meta - n_samples, "WARN",
            f"metadata has {n_meta} samples at {n_meta_sites} sites; the climate table "
            f"contributes {n_samples} at {s}. Every metric below describes the SMALLER "
            "set — see climate_na_excluded.tsv for what was dropped and why")
    if unmapped:
        add("samples_unmapped_to_site", unmapped, "WARN",
            "climate rows whose sample is absent from metadata — check the sample sets")
    if unparseable:
        add("samples_unparseable_predictors", unparseable, "WARN",
            "climate rows with an empty or non-numeric predictor cell; excluded from "
            "every count and every statistic below")

    counts = sorted(site_n.get(k, 0) for k in sites)
    add("samples_per_site_min", counts[0] if counts else 0)
    add("samples_per_site_median", counts[len(counts) // 2] if counts else 0)
    add("samples_per_site_max", counts[-1] if counts else 0)
    n_single = sum(1 for c in counts if c == 1)
    add("n_singleton_sites", n_single,
        "WARN" if n_single else "",
        "a site with n=1 contributes an environmental point with no within-site "
        "replication — any per-site variance estimate is undefined there")

    ratio = n_samples / s if s else float("nan")
    env_df = s - 1
    # Replication per site is only a problem when there are few sites to begin
    # with: 5 samples at each of 60 sites is a fine design, 5 at each of 11 is
    # 58 rows carrying 10 environmental observations.
    add("samples_per_environmental_point", round(ratio, 3),
        "WARN" if (ratio >= 3.0 and env_df < 20) else "",
        "sample-level GEA treats these as independent observations of the "
        "environment; they are not. Inflation of the effective test count scales "
        "with this ratio (design pseudoreplication, invisible to lambda_GC)")

    # Within-site share of predictor variance. Expected to be ~0 (one raster cell
    # per site); anything else means climate was extracted per SAMPLE coordinate,
    # so the environmental d.f. below understates the design.
    ss_within = ss_total = 0.0
    for j in range(p):
        allv = [row[j] for st in sites for row in per_site_rows[st]]
        if len(allv) < 2:
            continue
        gm = sum(allv) / len(allv)
        ss_total += sum((v - gm) ** 2 for v in allv)
        for st in sites:
            col = [row[j] for row in per_site_rows[st]]
            sm = sum(col) / len(col)
            ss_within += sum((v - sm) ** 2 for v in col)
    wfrac = (ss_within / ss_total) if ss_total > 0 else 0.0
    add("within_site_variance_frac", round(wfrac, 6),
        "WARN" if wfrac > 0.01 else "",
        "share of predictor variance sitting WITHIN sites; ~0 means one climate "
        "value per site (the usual case) and the environmental d.f. below is exact")

    add("environmental_df", env_df,
        "FAIL" if env_df < 3 else ("WARN" if env_df < 8 else ""),
        "n_sites - 1. This, not n_samples - 1, bounds any site-level "
        "climate/geography model (varpart, dbMEM, RDA constrained rank)")

    resid_df = s - p - 1
    add("site_level_residual_df", resid_df,
        "FAIL" if resid_df <= 0 else ("WARN" if resid_df < 3 else ""),
        "n_sites - n_predictors - 1. <= 0 means the predictor block saturates the "
        "site-level design: a site-level model fits it exactly and every adjusted "
        "R^2 / variance partition on it is meaningless")

    # Screen out zero-variance predictors BEFORE the spectrum. corr_matrix gives an
    # invariant column a synthetic identity row (sd == 0 => 1 on the diagonal, 0
    # off it), which contributes an eigenvalue of exactly 1 — so an invariant
    # predictor would silently ADD to numerical_rank and to effective
    # dimensionality, the two numbers most likely to be quoted. They are reported
    # separately instead, and the spectrum describes only the predictors that vary.
    cols_all = [[site_vals[k][j] for k in sites] for j in range(p)]
    _, sd_all = corr_matrix(cols_all) if (s >= 2 and p >= 1) else ([], [0.0] * p)
    invariant = [predictors[j] for j in range(p) if sd_all[j] == 0.0]
    add("n_invariant_predictors", len(invariant),
        "FAIL" if invariant else "",
        "zero variance across sites: " + (",".join(invariant) if invariant else "none"))

    active_j = [j for j in range(p) if sd_all[j] != 0.0]
    p_act = len(active_j)
    if p_act != p:
        add("n_predictors_in_spectrum", p_act, "",
            "invariant predictors are excluded from every spectral metric below")

    if s >= 3 and p_act >= 2:
        cols = [cols_all[j] for j in active_j]
        R, sd = corr_matrix(cols)
        p = p_act

        offd = [abs(R[i][j]) for i in range(p) for j in range(i + 1, p)]
        add("max_abs_pairwise_r", round(max(offd), 4) if offd else 0.0,
            "WARN" if offd and max(offd) > 0.9 else "",
            "predictor collinearity at SITE level; a varpart climate/geography "
            "split is not comparable across projects unless this block was "
            "decorrelated first")
        add("mean_abs_pairwise_r", round(sum(offd) / len(offd), 4) if offd else 0.0)

        ev = jacobi_eigenvalues(R)
        tot = sum(ev)
        add("env_pc1_variance_frac", round(ev[0] / tot, 4) if tot else float("nan"),
            "WARN" if tot and ev[0] / tot > 0.6 else "",
            "one dominant environmental axis leaves a multivariate method "
            "(RDA, geometric offset) little to composite")

        # A site-level correlation matrix of p predictors over s sites has rank at
        # most min(p, s - 1) — the centring costs one d.f. So a rank below p is NOT
        # evidence of aliased columns whenever p >= s; it is arithmetic, and saying
        # "aliased" there would be a false diagnosis. Distinguish the two.
        rank_ceiling = min(p, s - 1)
        pos = [e for e in ev if e > 1e-10]
        if len(pos) < rank_ceiling:
            _rank_note = (f"{len(pos)} < the structural ceiling min(p={p}, n_sites-1={s - 1}) "
                          f"= {rank_ceiling}: columns are aliased beyond what the design "
                          "already forces")
            _rank_flag = "WARN"
        elif rank_ceiling < p:
            _rank_note = (f"capped at min(p={p}, n_sites-1={s - 1}) = {rank_ceiling} BY THE "
                          "DESIGN, not by the data: with p >= n_sites the block is "
                          "rank-deficient for any predictor set whatsoever. This is the same "
                          "fact as site_level_residual_df <= 0, not a second one")
            _rank_flag = "WARN"
        else:
            _rank_note = f"full rank at {p} predictors over {s} sites"
            _rank_flag = ""
        add("numerical_rank", len(pos), _rank_flag, _rank_note)
        add("numerical_rank_ceiling", rank_ceiling, "",
            "min(n_predictors, n_sites - 1); no predictor set can exceed this here")

        # Participation ratio: (sum l)^2 / sum l^2. Equals p for a spherical
        # block, 1 for a perfectly collinear one.
        pr = (tot ** 2) / sum(e * e for e in ev) if ev else float("nan")
        add("effective_dimensionality", round(pr, 3),
            "WARN" if pr < max(2.0, p / 3.0) else "",
            "participation ratio of the correlation eigenvalues; the number of "
            "predictors the design actually resolves")

        # Jacobi on a singular correlation matrix returns lambda_min at the noise
        # floor, and it can be slightly NEGATIVE. Comparing that raw to 1e-12 leaves
        # a negative condition number that passes both the FAIL and the WARN test —
        # the worst-conditioned block in the project would come out unflagged.
        lam_min = ev[-1]
        cond = (ev[0] / lam_min) if lam_min > 1e-12 else float("inf")
        add("condition_number", ("inf" if cond == float("inf") else round(cond, 1)),
            "FAIL" if cond > 1e4 else ("WARN" if cond > 1e2 else ""),
            "lambda_max / lambda_min of the site-level CORRELATION matrix (invariant "
            "predictors excluded). Above ~1e4 the smallest eigen-directions carry "
            "numerical noise. NOTE: geometric_offset.R reports its own condition "
            "number on the per-sample COVARIANCE matrix, which is the one LEA "
            "actually inverts — the two are not expected to agree")

        cum, k95 = 0.0, 0
        for e in ev:
            cum += e
            k95 += 1
            if tot and cum / tot >= 0.95:
                break
        add("n_axes_to_95pct_variance", k95)
    else:
        add("spectral_diagnostics", "skipped", "WARN",
            f"needs >= 3 sites and >= 2 VARYING predictors (have {s} sites, "
            f"{p_act} varying of {p} predictors)")

    fails = [r for r in rows if r[2] == "FAIL"]
    warns = [r for r in rows if r[2] == "WARN"]
    add("n_flags_fail", len(fails))
    add("n_flags_warn", len(warns))

    with open(out_path, "w") as fh:
        fh.write("metric\tvalue\tflag\tnote\n")
        for metric, value, flag, note in rows:
            fh.write(f"{metric}\t{value}\t{flag}\t{note}\n")

    for metric, value, flag, note in rows:
        if flag:
            sys.stderr.write(f"{flag}: {metric} = {value} — {note}\n")
    sys.stderr.write(
        f"INFO: design adequacy written to {out_path} "
        f"({len(fails)} FAIL, {len(warns)} WARN)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
