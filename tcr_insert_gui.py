#!/usr/bin/env python3
"""
tcr_insert_gui.py  —  TCR Lentiviral Insert Generator (GUI)
════════════════════════════════════════════════════════════════════
A complete, one-window tool that goes from a TCR's V gene, J gene, and
CDR3 straight to a finished, colour-coded annotated lentiviral insert.

It does the entire pipeline automatically:
    1. Runs Stitchr (with -nl, no leader) for the alpha and beta chains
    2. Trims the constant regions
         alpha -> remainder 0 (cut before human TRAC)
         beta  -> remainder 1, keeping the terminal G (cut before TRBC)
    3. Builds the lentiviral insert (embedded lab logic)
    4. Writes a colour-coded annotated document (.html, opens in any
       browser and prints/saves to PDF) plus a plain .txt of the
       nucleotide + amino-acid sequence.

You type the six values, click one button, and the annotated document
opens. No manual trimming, no editing spreadsheets.

──────────────────────────────────────────────────────────────────
REQUIREMENTS
    • Python 3  (pre-installed on Mac)
    • Stitchr installed:   pip install stitchr
                           stitchrdl -s HUMAN
    • Tkinter (usually included; on Homebrew Python:
                           brew install python-tk)

RUN
    python3 tcr_insert_gui.py
════════════════════════════════════════════════════════════════════
"""

import sys, os, re, subprocess, webbrowser, tempfile, random, datetime

# ── Tkinter (GUI) ───────────────────────────────────────────────────
try:
    import tkinter as tk
    from tkinter import ttk, messagebox, filedialog
except ImportError:
    sys.exit("ERROR: Tkinter is not installed.\n"
             "On a Homebrew Mac, run:  brew install python-tk\n"
             "Then try again.")

# ════════════════════════════════════════════════════════════════════
# CODON TABLE + SEQUENCE UTILITIES
# ════════════════════════════════════════════════════════════════════

CODONS = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_','TGC':'C','TGT':'C','TGA':'_','TGG':'W',
}
AAS = {}
for _c, _a in CODONS.items():
    AAS.setdefault(_a, []).append(_c)

def translate(seq, frame=0):
    seq = seq.upper()[frame:]
    return ''.join(CODONS.get(seq[i:i+3], '?') for i in range(0, (len(seq)//3)*3, 3))

def find_all(seq, sub):
    i, locs = 0, []
    while True:
        p = seq.find(sub, i)
        if p == -1: break
        locs.append(p); i = p + 1
    return locs

def find_restriction_sites(sites, seq):
    return {s: find_all(seq, sites[s]) for s in sites}

def swap_codons(s, site_name, idx):
    taboo = RESTRICTION_SITES[site_name]
    first = 3 * (idx // 3)
    last  = 3 * ((idx + len(taboo)) // 3)
    ci = random.choice(list(range(first, last + 3, 3)))
    codon = s[ci:ci+3].upper()
    aa = CODONS.get(codon, '?')
    nc = codon
    for cand in sorted(AAS.get(aa, [])):
        if cand != codon:
            nc = cand; break
    return s[:ci] + nc + s[ci+3:]

# ════════════════════════════════════════════════════════════════════
# INSERT GENERATION  (embedded, identical to generate_lentiviral_insert.py)
# ════════════════════════════════════════════════════════════════════

RESTRICTION_SITES = {'MfeI':'CAATTG','AccI':'GTATAC','NotI':'GCGGCCGC',
                     'SnaBI':'TACGTA','SpeI':'ACTAGT'}
CR = [
    'tctcgagta',
    'catgggctccaggctgctctgttgggtgctgctttgtctcctgggagcaggcccagtaaaggctgga',
    'aagatc',
    'acgtgacaccacccaaagtctcactgtttgagcctagcaaggcagaaattgccaacaagcagaagg'
    'ccaccctggtgtgcctggcaagagggttctttccagatcacgtggagctgtcctggtgggtcaacg'
    'gcaaagaagtgcattctggggtctgcaccgacccccaggcttacaaggagagtaattactcatattg'
    'tctgtcaagccggctgagagtgtccgccacattctggcacaaccctaggaatcatttccgctgccag'
    'gtccagtttcacggcctgagtgaggaagataaatggccagaggggtcacctaagccagtgacacaga'
    'acatcagcgcagaagcctggggacgagcagactgtggcattactagcgcctcctatcatcagggcgt'
    'gctgagcgccactatcctgtacgagattctgctgggaaaggccaccctgtatgctgtgctggtctcc'
    'ggcctggtgctgatggccatggtcaagaaaaagaactctgggagtggagccacaaatttctctctgc'
    'tgaaacaggctggagatgtggaggaaaaccccggccctatgaagagcctgcgcgtgctgctggtcat'
    'cctgtggctg',
    'tcgtgggtctggagccaa',
    'gacattcagaacccggaaccggct',
    'cagctgaaggacccccgatctcaggatagtactctgtgcctgttcaccgactttgatagtcagatca'
    'atgtgcctaaaaccatggaatccggaacttttattaccgacaagtgcgtgctggatatgaaagccat'
    'ggacagtaagtcaaacggcgccatcgcttggagcaatcagacatccttcacttgccaggatatcttc'
    'aaggagaccaacgcaacatacccatcctctgacgtgccctgtgatgccaccctgacagagaagtctt'
    'tcgaaacagacatgaacctgaattttcagaatctgagcgtgatgggcctgagaatcctgctgctgaa'
    'ggtcgctgggtttaatctgctgatgacactgcggctgtggtcctcatgaattcggaccgtgtccaat'
    'gtagc',
    'gtcgacaatcaacctctgga',
]

def trim_by_remainder(seq, factor=3, remainder=0):
    r = len(seq) % factor
    trimlen = r - remainder
    if trimlen < 0: trimlen += 3
    return seq if trimlen == 0 else seq[:-trimlen]

def generate_insert(achain, bchain):
    achain = trim_by_remainder(achain, 3, 0)
    bchain = trim_by_remainder(bchain, 3, 1)
    exp = {}
    ins = CR[0]
    exp['NotI'] = len(ins);  ins += RESTRICTION_SITES['NotI'] + CR[1] + bchain + CR[2]
    exp['SnaBI'] = len(ins); ins += RESTRICTION_SITES['SnaBI'] + CR[3]
    exp['MfeI'] = len(ins);  ins += RESTRICTION_SITES['MfeI'] + CR[4] + achain + CR[5]
    exp['AccI'] = len(ins);  ins += RESTRICTION_SITES['AccI'] + CR[6]
    exp['SpeI'] = len(ins);  ins += RESTRICTION_SITES['SpeI'] + CR[7]
    ins = ins.replace('aaggctggaaaggcta', 'aaggctgga')
    ins_aa = translate(ins, 0)
    guard = 0
    while guard < 200:
        guard += 1
        ok = True
        locs = find_restriction_sites(RESTRICTION_SITES, ins)
        for s in exp:
            for i in locs.get(s, []):
                if i != exp[s]:
                    ok = False
                    ins = swap_codons(ins, s, i)
        if ok: break
    ins = 'a' + ins
    return ins, ins_aa

# ════════════════════════════════════════════════════════════════════
# STITCHR + TRIMMING
# ════════════════════════════════════════════════════════════════════

STITCHR_CMDS = ['stitchr', '/opt/homebrew/bin/stitchr', '/usr/local/bin/stitchr']

def _extract_nt_from_stitchr(output):
    """
    Parse Stitchr's FASTA output and return the full nucleotide sequence.
    Stitchr wraps the sequence across many 60-char lines, so we must collect
    ALL sequence lines after the '>nt' header until the next blank line or
    the '>aa' header — not just the first line.
    """
    lines = output.splitlines()
    seq_lines = []
    capturing = False
    for ln in lines:
        s = ln.strip()
        if s.startswith('>'):
            if capturing:
                break                    # reached the >aa header — stop
            if s.lower().startswith('>nt') or ('nt|' in s.lower()):
                capturing = True
                continue
            # some builds emit a single '>' header; start capturing after it
            capturing = True
            continue
        if capturing:
            if s == '':
                if seq_lines:
                    break                # blank line ends the nt block
                continue
            if re.fullmatch(r'[ACGTNacgtn]+', s):
                seq_lines.append(s)
            else:
                # a non-sequence line (e.g. a divider) — stop if we have data
                if seq_lines:
                    break
    return ''.join(seq_lines).upper() if seq_lines else None


def run_stitchr(v, j, cdr3, locus):
    """Run stitchr -nl (no leader). Returns V/CDR3/J/CONSTANT nt string."""
    last = ""
    for base in STITCHR_CMDS:
        cmd = [base, '-v', v, '-j', j, '-cdr3', cdr3, '-s', 'HUMAN', '-l', locus, '-nl']
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            last = str(e); continue
        seq = _extract_nt_from_stitchr(r.stdout)
        if seq and len(seq) > 30:
            return seq
        last = (r.stdout + r.stderr)[:250]
    raise RuntimeError(
        f"Stitchr failed for {locus} ({v} / {j} / {cdr3}).\n\n"
        f"Message: {last}\n\n"
        "Check that Stitchr is installed (pip install stitchr) and that "
        "'stitchrdl -s HUMAN' has been run. If the gene name has a * or /, "
        "the app already passes it safely — but confirm the gene exists in IMGT.")

ALPHA_MARKERS = ['DIQNPDP','DIQNPEP','IQNPDPAV','IQNPEPAV','DIQNP']
# Human TRBC (beta constant) is highly conserved right after the opening E.
# These anchors are long and constant-region-specific so they cannot match
# inside the V/CDR3/J region by accident. We search for the anchor, then step
# back to the E that opens the constant region (keeping the terminal G in nt).
BETA_C_ANCHORS = ['DLNKVFPPEVAVF','DLRNVTPPKVSLF','DLKNVFPPEVAVF',
                  'DLNKVFPPEVAV','DLRNVTPPKVSL','DLKNVFPPEVAV',
                  'FPPEVAVFEPSE','PKVSLFEPSKAE']
# Fallbacks (shorter) only used if the long anchors all miss:
BETA_C_FALLBACK = ['EDLNKV','EDLRNV','EDLKNV']

def _aa_to_nt_cut(prot_index):
    """Convert an amino-acid index to the nucleotide index that starts it."""
    return prot_index * 3

def trim_alpha(nt):
    prot = translate(nt)
    for m in ALPHA_MARKERS:
        i = prot.find(m)
        if i != -1:
            return nt[:i*3]              # remainder 0
    p = nt.upper().find('ATATCCAGAA')
    if p == -1: p = nt.upper().find('ATAT')
    return nt[:p] if p != -1 else nt

def trim_beta(nt):
    """
    Cut off the human TRBC constant region, keeping the terminal G (first nt of
    the E codon that opens the constant region). Uses long, constant-specific
    anchors searched near the 3' end so a coincidental 'EDL' in the V/CDR3 region
    cannot cause an early cut.
    """
    prot = translate(nt)
    # 1) Try the long anchors. Each anchor begins with 'DL...' — the E is the
    #    residue immediately before it. Use the LAST occurrence (rightmost).
    for anchor in BETA_C_ANCHORS:
        i = prot.rfind(anchor)
        if i != -1:
            e_index = i - 1            # the E sits just before 'DL...'
            if e_index >= 0 and prot[e_index] == 'E':
                return nt[:e_index*3 + 1]   # keep the G (first nt of E codon)
            # if the residue before isn't E, still cut at the anchor's D, +1 for G
            return nt[:i*3 + 1]
    # 2) Fallback: the short 'EDLxxx' markers, but taken as the LAST occurrence
    #    and only if they appear in the last ~40% of the protein (constant region
    #    is always near the 3' end) to avoid a spurious early hit.
    min_pos = int(len(prot) * 0.55)
    best = -1
    for m in BETA_C_FALLBACK:
        j = prot.rfind(m)
        if j >= min_pos and j > best:
            best = j
    if best != -1:
        return nt[:best*3 + 1]         # 'E' is first residue of marker; keep its G
    # 3) Last resort: return unchanged (caller will warn on remainder/G check)
    return nt

# ════════════════════════════════════════════════════════════════════
# ANNOTATION -> HTML (colour-coded, opens in browser, prints to PDF)
# ════════════════════════════════════════════════════════════════════

COLORS = {
    'notI':'#B8860B','beta':'#1B5E20','linker':'#37474F','mfei':'#6A1B9A',
    'alpha':'#BF360C','bstz17i':'#B71C1C','spei':'#00695C',
    'bracket':'#B71C1C','gray':'#9E9E9E',
}

def annotate_segments(insert_nt):
    s = insert_nt
    segs = []
    # The insert begins with a single 'a' (prepended by generate_insert) then CR[0]
    # then NotI. The bracket should open BEFORE that 'a' so the 'a' is inside the
    # synthesis region: [a tctcgagta GCGGCCGC ...
    bracket_start = s.find('GCGGCCGC') - len(CR[0]) - 1
    segs.append((s[:bracket_start], COLORS['gray']))   # anything before the bracket (usually empty)
    segs.append(('[', COLORS['bracket']))
    rest = s[bracket_start:]                            # includes the leading 'a'
    def take(n, col):
        nonlocal rest
        segs.append((rest[:n], col)); rest = rest[n:]
    take(1 + len(CR[0]), COLORS['linker'])             # the 'a' + CR[0] (both inside bracket)
    take(8, COLORS['notI'])
    take(len(CR[1]), COLORS['linker'])
    take(rest.find(CR[2]), COLORS['beta'])
    take(len(CR[2] + RESTRICTION_SITES['SnaBI'] + CR[3]), COLORS['linker'])
    take(6, COLORS['mfei'])
    take(len(CR[4]), COLORS['linker'])
    take(rest.find(CR[5]), COLORS['alpha'])
    take(len(CR[5]), COLORS['linker'])
    take(6, COLORS['bstz17i'])
    take(10, COLORS['linker'])          # first 10 nt of muTrac (inside bracket)
    segs.append((']', COLORS['bracket']))   # close synthesis region here (matches reference)
    take(len(CR[6]) - 10, COLORS['gray'])   # rest of muTrac (outside bracket)
    take(6, COLORS['gray'])                 # SpeI (outside)
    segs.append((rest, COLORS['gray']))     # CR[7] tail (outside)
    return segs

def seam_check(insert_aa, cdr3_b):
    i = insert_aa.find(cdr3_b)
    if i == -1:
        return "CDR3 not found", False
    seam = insert_aa[i+len(cdr3_b):i+len(cdr3_b)+16]
    # Correct seam has the glutamate (E) opening the constant region: ...E DLRNV / DLKNV
    ok = ('EDLR' in seam or 'EDLK' in seam or 'EDLN' in seam)
    return seam, ok

def build_html(name, desc, system,
               trav, traj, cdr3_a, trbv, trbj, cdr3_b,
               insert_nt, insert_aa):
    segs = annotate_segments(insert_nt)
    seam, ok = seam_check(insert_aa, cdr3_b)

    # coloured nucleotide HTML (monospace, wraps naturally)
    nt_html = ''
    for text, col in segs:
        if not text: continue
        safe = text.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
        weight = 'bold' if col in (COLORS['notI'],COLORS['bstz17i'],
                                   COLORS['spei'],COLORS['bracket'],COLORS['mfei']) else 'normal'
        nt_html += f"<span style='color:{col};font-weight:{weight}'>{safe}</span>"

    aa_main = insert_aa.split('_')[0]
    aa_tail = insert_aa[len(aa_main):]

    seam_color = '#1B5E20' if ok else '#B71C1C'
    seam_word  = 'PASSED' if ok else 'CHECK'

    legend = [
        (COLORS['bracket'], '[ ]', 'Genscript synthesis region'),
        (COLORS['notI'],   'NotI', "5' cloning site"),
        (COLORS['beta'],   'Beta (TRB)', 'Clone-specific beta chain'),
        (COLORS['linker'], 'Scaffold', 'Signal peptides, muTrbc, P2A, muTrac'),
        (COLORS['mfei'],   'MfeI', 'Internal site'),
        (COLORS['alpha'],  'Alpha (TRA)', 'Clone-specific alpha chain'),
        (COLORS['bstz17i'],'BstZ17I', "3' cloning site"),
        (COLORS['spei'],   'SpeI', 'Internal site'),
        (COLORS['gray'],   'Gray', 'Vector context (not synthesised)'),
    ]
    legend_html = ''
    for col, lab, d in legend:
        legend_html += (f"<div class='key'><span class='sw' style='background:{col}'></span>"
                        f"<b style='color:{col}'>{lab}</b> &nbsp;{d}</div>")

    today = datetime.date.today().isoformat()
    return f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>{name} — Annotated Insert</title>
<style>
  body {{ font-family: Arial, Helvetica, sans-serif; margin: 40px; color:#12293d; }}
  h1 {{ font-size: 20px; margin-bottom: 2px; }}
  .meta {{ font-size: 13px; color:#444; margin: 2px 0; }}
  .desc {{ font-style: italic; color:#666; font-size:13px; }}
  .seam {{ font-weight:bold; font-size:12.5px; margin:8px 0; color:{seam_color}; }}
  hr {{ border:none; border-top:1px solid #c0cdd8; margin:12px 0; }}
  .seq {{ font-family: 'Courier New', monospace; font-size: 12px;
          line-height: 1.55; word-break: break-all; }}
  .sectitle {{ font-weight:bold; color:#1B3A5C; margin:14px 0 4px; font-size:13px; }}
  .key {{ font-size:12px; margin:2px 0; }}
  .sw {{ display:inline-block; width:11px; height:11px; border-radius:2px;
         margin-right:5px; vertical-align:middle; }}
  .legendbox {{ columns:2; margin:10px 0; }}
  .tail {{ color:#bbb; }}
  @media print {{ body {{ margin: 18px; }} }}
</style></head><body>
<h1>Annotated Lentiviral Insert — {name}</h1>
<div class='desc'>{desc}</div>
<div class='meta'><b>TRA:</b> {trav} / {traj} / <span style="font-family:monospace">{cdr3_a}</span>
&nbsp;&nbsp; <b>TRB:</b> {trbv} / {trbj} / <span style="font-family:monospace">{cdr3_b}</span></div>
<div class='meta'>{system} &nbsp;·&nbsp; generated {today}</div>
<div class='seam'>Seam check {seam_word}: beta J/constant junction reads
…<span style="font-family:monospace">{cdr3_b}·{seam}</span>…</div>
<hr>
<div class='legendbox'>{legend_html}</div>
<hr>
<div class='sectitle'>Annotated nucleotide sequence</div>
<div class='seq'>{nt_html}</div>
<div class='sectitle'>Translated amino acid sequence</div>
<div class='seq'>{aa_main}</div>
<hr>
<div style='font-size:11px;color:#888;font-style:italic'>
Generated by tcr_insert_gui.py · Insert length {len(insert_nt)} nt ·
To save as PDF: File → Print → Save as PDF</div>
</body></html>"""

# ════════════════════════════════════════════════════════════════════
# PIPELINE (used by the GUI button)
# ════════════════════════════════════════════════════════════════════

def run_pipeline(name, desc, system, trav, traj, cdr3_a, trbv, trbj, cdr3_b, out_dir):
    # 1. Stitchr
    a_full = run_stitchr(trav, traj, cdr3_a, 'TRA')
    b_full = run_stitchr(trbv, trbj, cdr3_b, 'TRB')
    # 2. Trim
    a_trim = trim_alpha(a_full)
    b_trim = trim_beta(b_full)
    warns = []
    if len(a_trim) % 3 != 0:
        warns.append(f"alpha remainder is {len(a_trim)%3} (expected 0)")
    if len(b_trim) % 3 != 1:
        warns.append(f"beta remainder is {len(b_trim)%3} (expected 1)")
    if not b_trim.upper().endswith('G'):
        warns.append("beta does not end in G")
    # 3. Insert
    ins_nt, ins_aa = generate_insert(a_trim, b_trim)
    # sanity
    if cdr3_a not in ins_aa: warns.append(f"alpha CDR3 {cdr3_a} not found in insert")
    if cdr3_b not in ins_aa: warns.append(f"beta CDR3 {cdr3_b} not found in insert")
    # 4. Outputs
    safe = re.sub(r'[^A-Za-z0-9_.-]', '_', name) or "insert"
    html = build_html(name, desc, system, trav, traj, cdr3_a,
                      trbv, trbj, cdr3_b, ins_nt, ins_aa)
    html_path = os.path.join(out_dir, f"{safe}_annotated_insert.html")
    txt_path  = os.path.join(out_dir, f"{safe}_output.txt")
    with open(html_path, 'w') as f:
        f.write(html)
    with open(txt_path, 'w') as f:
        f.write(f"Sample: {name}\nSystem: {system}\n")
        f.write(f"TRA: {trav} / {traj} / {cdr3_a}\nTRB: {trbv} / {trbj} / {cdr3_b}\n\n")
        f.write("Insert nucleotide sequence:\n" + ins_nt + "\n\n")
        f.write("Insert amino acid sequence:\n" + ins_aa + "\n")
    seam, ok = seam_check(ins_aa, cdr3_b)
    return html_path, txt_path, seam, ok, warns

# ════════════════════════════════════════════════════════════════════
# GUI
# ════════════════════════════════════════════════════════════════════

class App:
    def __init__(self, root):
        self.root = root
        root.title("TCR Lentiviral Insert Generator")
        root.geometry("760x620")
        root.configure(bg="#f4f6f9")

        pad = {'padx': 8, 'pady': 4}
        head = tk.Label(root, text="TCR Lentiviral Insert Generator",
                        font=("Helvetica", 17, "bold"), bg="#f4f6f9", fg="#12293d")
        head.pack(pady=(14, 2))
        sub = tk.Label(root, text="V / J / CDR3  →  Stitchr  →  trim  →  insert  →  annotated document",
                       font=("Helvetica", 10), bg="#f4f6f9", fg="#555")
        sub.pack(pady=(0, 10))

        form = tk.Frame(root, bg="#f4f6f9")
        form.pack(fill="x", padx=24)

        # Sample + system
        self.e = {}
        def row(label, key, default="", width=46):
            fr = tk.Frame(form, bg="#f4f6f9")
            fr.pack(fill="x", **pad)
            tk.Label(fr, text=label, width=14, anchor="e", bg="#f4f6f9",
                     font=("Helvetica", 10)).pack(side="left")
            ent = tk.Entry(fr, width=width, font=("Courier", 11))
            ent.insert(0, default)
            ent.pack(side="left", fill="x", expand=True)
            self.e[key] = ent

        row("Sample name", "name", "MyClone_Rank1")
        row("Description", "desc", "")
        row("System label", "system", "DR4 / clone 461")

        tk.Frame(form, height=8, bg="#f4f6f9").pack()
        alpha_hdr = tk.Label(form, text="Alpha chain (TRA)", font=("Helvetica", 11, "bold"),
                             bg="#f4f6f9", fg="#BF360C", anchor="w")
        alpha_hdr.pack(fill="x", padx=8)
        row("TRAV gene", "trav", "TRAV41")
        row("TRAJ gene", "traj", "TRAJ49")
        row("TRA CDR3", "cdr3_a", "CAAAGNQFYF")

        tk.Frame(form, height=8, bg="#f4f6f9").pack()
        beta_hdr = tk.Label(form, text="Beta chain (TRB)", font=("Helvetica", 11, "bold"),
                            bg="#f4f6f9", fg="#1B5E20", anchor="w")
        beta_hdr.pack(fill="x", padx=8)
        row("TRBV gene", "trbv", "TRBV19")
        row("TRBJ gene", "trbj", "TRBJ2-3")
        row("TRB CDR3", "cdr3_b", "CASGRSHGTDTQYF")

        # Output folder
        of = tk.Frame(form, bg="#f4f6f9")
        of.pack(fill="x", **pad)
        tk.Label(of, text="Save to", width=14, anchor="e", bg="#f4f6f9",
                 font=("Helvetica", 10)).pack(side="left")
        self.out_var = tk.StringVar(value=os.path.expanduser("~/Desktop"))
        tk.Entry(of, textvariable=self.out_var, font=("Courier", 10)).pack(
            side="left", fill="x", expand=True)
        tk.Button(of, text="Browse…", command=self.browse).pack(side="left", padx=4)

        # Run button
        self.btn = tk.Button(root, text="Generate Insert",
                             font=("Helvetica", 13, "bold"),
                             bg="#1a3d5c", fg="white", activebackground="#24557e",
                             relief="flat", padx=20, pady=8, command=self.go)
        self.btn.pack(pady=14)

        # Status
        self.status = tk.Text(root, height=6, font=("Courier", 9),
                              bg="#0f1a24", fg="#d0e0f0", relief="flat", wrap="word")
        self.status.pack(fill="both", expand=True, padx=24, pady=(0, 16))
        self.log("Ready. Fill in the fields and click Generate Insert.\n"
                 "The alpha and beta chains are built with Stitchr automatically — "
                 "no manual trimming needed.")

    def log(self, msg, clear=False):
        if clear:
            self.status.delete("1.0", "end")
        self.status.insert("end", msg + "\n")
        self.status.see("end")
        self.root.update()

    def browse(self):
        d = filedialog.askdirectory(initialdir=self.out_var.get())
        if d: self.out_var.set(d)

    def go(self):
        vals = {k: self.e[k].get().strip() for k in self.e}
        out_dir = self.out_var.get().strip()
        required = ['name','trav','traj','cdr3_a','trbv','trbj','cdr3_b']
        missing = [k for k in required if not vals[k]]
        if missing:
            messagebox.showerror("Missing fields",
                                 "Please fill in: " + ", ".join(missing))
            return
        if not os.path.isdir(out_dir):
            messagebox.showerror("Bad folder", "The output folder does not exist.")
            return

        self.btn.config(state="disabled")
        self.log("", clear=True)
        self.log(f"► {vals['name']}")
        self.log(f"  TRA: {vals['trav']} / {vals['traj']} / {vals['cdr3_a']}")
        self.log(f"  TRB: {vals['trbv']} / {vals['trbj']} / {vals['cdr3_b']}")
        self.log("  Running Stitchr and building insert…")
        try:
            html_path, txt_path, seam, ok, warns = run_pipeline(
                vals['name'], vals['desc'], vals['system'],
                vals['trav'], vals['traj'], vals['cdr3_a'].upper(),
                vals['trbv'], vals['trbj'], vals['cdr3_b'].upper(), out_dir)
            self.log(f"  Beta seam: …{seam}  [{'OK' if ok else 'CHECK THIS'}]")
            if warns:
                for w in warns:
                    self.log("  ⚠ " + w)
            self.log(f"  ✓ Annotated document: {os.path.basename(html_path)}")
            self.log(f"  ✓ Sequence text file: {os.path.basename(txt_path)}")
            self.log("  Opening the annotated document in your browser…")
            webbrowser.open('file://' + os.path.abspath(html_path))
            if not ok:
                messagebox.showwarning(
                    "Seam check",
                    "The insert was generated, but the beta junction did not match the "
                    "expected pattern. Open the document and check the seam before using it.")
            else:
                messagebox.showinfo(
                    "Done",
                    f"Insert generated successfully.\n\nSaved to:\n{out_dir}\n\n"
                    "The annotated document opened in your browser. "
                    "To save it as a PDF: File → Print → Save as PDF.")
        except Exception as ex:
            self.log("  ✗ ERROR: " + str(ex))
            messagebox.showerror("Error", str(ex))
        finally:
            self.btn.config(state="normal")


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()

if __name__ == '__main__':
    main()
