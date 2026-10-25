#!/usr/bin/env python3
"""Annotation-only guesses for the +-500 kb neighbourhood genes of the 15 EMMAX+RDA regions of Trifolium86MainK5ERpc5
(plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of ANALYSIS/Trifolium86MainK5/regionsER/
write_neighbourhood_guess.py (PC3 run; user 2026-10-14: "just your guess of the adaptive genes there"): the PC3 guesses
are kept except g31011 / g31002 (now genes of the new 4_101.3M region, curated) and the five of the two dropped chr-6
regions (g44891, g44898, g44912, g44915, g45852); guesses for the four new regions are added.

UNVERIFIED. Read from the eggNOG description / preferred name / PFAM / KEGG name only: no literature, no SNP analysis.
Neighbourhood = transcripts within region midpoint +- 500,000 bp that are not genes of that region
(tables/neighbourhood.tsv). Genes that are genes of ANOTHER of the 15 regions (4_100.9M / 4_101.3M overlap) are
skipped here: they are curated in tables/curation.tsv.

Output: tables/neighbourhood_guess.tsv (only the genes worth a look; priority A = a named class with an established
abiotic or phenology role, B = family-level link only)
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/write_neighbourhood_guess.py
"""
import csv
import os

B = os.path.dirname(os.path.abspath(__file__))
G = {  # gene -> (priority, guess)
    # 1_8649433-8911055
    "g740.t1": ("B", "drought: BURP cluster continues left of the region (15-copy array 1:557)"),
    "g755.t1": ("A", "drought: RD22-like BURP 'dehydration-responsive' (left of the region; was inside the old region)"),
    "g757.t1": ("A", "drought: RD22-like BURP 'dehydration-responsive' (left of the region; was inside the old region)"),
    "g793.t1": ("A", "light / temperature: phytochrome E (PHYE), tandem pair with g794"),
    "g794.t1": ("A", "light / temperature: phytochrome E (PHYE), tandem pair with g793"),
    # 3_147045758-147307380
    "g25499.t1": ("B", "drought?: double Clp-N P-loop NTPase = SMAX1 / SMXL-like (karrikin / strigolactone repressor)"),
    "g25541.t1": ("B", "ABA?: bZIP array (g25541, g25551-g25554); check for ABF / AREB / ABI5 members"),
    # 4_8145960-8407582
    "g26497.t1": ("A", "ABA / drought: GBF3 G-box-binding bZIP"),
    "g26512.t1": ("A", "photoprotection: violaxanthin de-epoxidase (VDE / NPQ1), high light + heat"),
    "g26589.t1": ("A", "multi-stress: A20/AN1 stress-associated protein (SAP); old call stress"),
    "g26576.t1": ("B", "light: light-inducible bZIP pair (g26576-7), HY5 / HYH-like?"),
    # 4_11490368-11751990
    "g26915.t1": ("A", "drought: Di19 (drought-induced protein 19) zinc finger"),
    "g26923.t1": ("A", "cuticle: 3-ketoacyl-CoA synthase KCS11 (very-long-chain fatty acids for wax)"),
    "g26972.t1": ("B", "K+ / salt: KUP / HAK potassium transporter pair (g26972-3)"),
    "g26984.t1": ("B", "redox: glutaredoxin array (g26984-9, GRXC11-like)"),
    "g26991.t1": ("A", "phenology: CONSTANS-like B-box + CCT"),
    "g26924.t1": ("B", "hypoxia / sugar: sucrose synthase SUS4"),
    "g26893.t1": ("B", "universal stress protein A-like"),
    # 4_15217488-15556144
    "g27338.t1": ("A", "water transport: TIP aquaporin, 46 kb left of the region"),
    "g27342.t1": ("A", "stomata: beta-carbonic anhydrase (CO2-controlled stomatal movement); 19 kb left"),
    # 4_100931231-101192853
    "g30966.t1": ("B", "class III peroxidase"),
    # 4_107336728-107600612
    "g31313.t1": ("B", "growth restraint under stress: gibberellin 2-oxidase (GA2ox)"),
    # 4_121924771-122186393
    "g32225.t1": ("A", "osmosensing: OSCA / CSC1-like hyperosmolality-gated Ca2+ channel (TMEM63), 11 bp right of the region"),
    "g32251.t1": ("A", "drought / senescence: NAC cluster of 11 copies (g32251-g32266); check for ATAF / RD26 / NAC072-like"),
    "g32183.t1": ("B", "cold: glycine-rich RNA-binding protein (GRP7 / GRP8-like?)"),
    "g32194.t1": ("B", "expansin-A pair (g32194-5)"),
    "g32164.t1": ("B", "isoflavonoid malonyltransferase / 7-IOMT cluster (legume specialised metabolism; UV / defence)"),
    # 5_5557888-5883272
    "g32831.t1": ("A", "stomatal density: EPIDERMAL PATTERNING FACTOR-like peptide"),
    "g32889.t1": ("A", "salt / methylglyoxal detox: glyoxalase I (GLO1)"),
    "g32885.t1": ("B", "carotenoid cleavage dioxygenase CCD1 pair (g32885-6), apocarotenoids"),
    "g32926.t1": ("B", "YUCCA (YUC6-like) monooxygenase"),
    # 5_111174564-111436186
    "g38759.t1": ("A", "heat: GSNOR / HOT5 S-nitrosoglutathione reductase, thermotolerance (cited in the old deep dive)"),
    "g38808.t1": ("A", "phenology?: agamous-like MADS-box; check for SOC1 / SVP / AGL24 membership"),
    "g38760.t1": ("B", "DnaJ chaperone"),
    # 6_116698620-117089945
    "g45735.t1": ("A", "phenology: VERNALIZATION INSENSITIVE 3-like (VIN3 / VEL); old call phenology"),
    "g45784.t1": ("A", "cold signalling: diacylglycerol kinase (DGK, cold-induced PA signalling)"),
    "g45717.t1": ("B", "class III peroxidase"),
    # new at 5 PCs: 2_145944878-146206511
    "g12488.t1": ("A", "ABA biosynthesis: abscisic-aldehyde oxidase (ko AAO3) array of 4 (g12488-g12490, g12492), 265-287 kb right"),
    "g12485.t1": ("B", "SPL7 (copper homeostasis; SPL family also in flowering)"),
    # new at 5 PCs: 2_176097944-176448788
    "g15007.t1": ("A", "osmosensing: OSCA / CSC1-like hyperosmolality-gated Ca2+ channel pair (g15007, g15011), 74-101 kb left"),
    "g15071.t1": ("A", "ABA / light: G-box-binding factor bZIP (ko GBF), 290 kb right"),
    # new at 5 PCs: 5_8217606-8479228
    "g33182.t1": ("A", "phenology: CONSTANS-like B-box + CCT, 121 kb right"),
    "g33206.t1": ("B", "cuticle: 3-ketoacyl-CoA synthase (KCS, very-long-chain fatty acids for wax)"),
    "g33140.t1": ("B", "AP2/ERF pair (g33140, g33208); check for DREB members"),
    "g33197.t1": ("B", "chaperonin Cpn60"),
}

nb = {r["gene_id"]: r for r in csv.DictReader(open(os.path.join(B, "tables/neighbourhood.tsv")), delimiter="\t")}
assert set(G) <= set(nb), set(G) - set(nb)
assert not any(nb[g]["in_other_region"] == "TRUE" for g in G)
with open(os.path.join(B, "tables/neighbourhood_guess.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["region_ids", "gene_id", "dist_to_edge", "priority", "guess", "desc", "old_call"])
    for g in sorted(G, key=lambda g: (nb[g]["chr"], int(nb[g]["start"]))):
        r = nb[g]
        w.writerow([r["region_ids"], g, r["dist_to_edge"], G[g][0], G[g][1], r["desc"], r["old_call"]])
print(f"INFO: {len(G)} guesses ({sum(1 for v in G.values() if v[0] == 'A')} priority A)")
