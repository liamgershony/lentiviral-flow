# lentiviral-flow
Simple end to end tool to generate annotated lentiviral TCR inserts from V/J/CDR3 input


You type the six values into a window, click a button, and the tool runs [Stitchr](https://github.com/JamieHeather/stitchr), trims the chains correctly, builds the lentiviral insert, and opens an annotated document in your browser.

No manual sequence trimming. No spreadsheets. No editing DNA by hand.

---

## What it does

For each TCR, the tool:

1. Runs Stitchr (with `-nl`, no leader) to reconstruct the alpha and beta chain nucleotide sequences.
2. Trims off the constant regions automatically — the alpha chain to remainder 0, and the beta chain keeping the terminal G so the junction glutamate (E) forms correctly.
3. Builds the lentiviral insert by inserting the chains into the vector scaffold (mouse constant regions, P2A, signal peptides, and the NotI / BstZ17I cloning sites).
4. Produces a colour-coded annotated document (opens in any browser, prints to PDF) plus a plain-text sequence file.

It also runs a **seam check** on the beta junction and flags anything that looks wrong before you use the sequence.

---

## Setup (one time)

You need **Python 3** (pre-installed on Mac) and a few things installed once.

**1. Install Stitchr and download the human gene data:**

```
pip install stitchr
stitchrdl -s human
```

**2. If you are on a Mac using Homebrew Python and the app says Tkinter is missing:**

```
brew install python-tk
```

That's it. You only do this once.

---

## Running the tool

**1. Open Terminal and point it at the folder** (type `cd `, then drag the folder onto the Terminal window, then press Enter):

```
cd ~/Desktop/lentiviral-insert-generator
```

**2. Launch the app:**

```
python3 tcr_insert_gui.py
```

**3. Fill in the fields** — a sample name, a system label, and the six TCR values:

| Field | Example |
|-------|---------|
| TRAV gene | `TRAV41` |
| TRAJ gene | `TRAJ49` |
| TRA CDR3 | `CAAAGNQFYF` |
| TRBV gene | `TRBV19` |
| TRBJ gene | `TRBJ2-3` |
| TRB CDR3 | `CASGRSHGTDTQYF` |

**4. Click "Generate Insert."** The annotated document opens automatically in your browser. Two files are saved to your chosen folder: the annotated document (`.html`) and a plain-text sequence file (`.txt`).

To save the annotated document as a PDF: in the browser, **File → Print → Save as PDF**.

---

## Notes on gene names

- You can enter gene names with or without the allele suffix — Stitchr defaults to `*01`. So `TRBV19` and `TRBV19*01` give the same result.
- For a specific non-`*01` allele, type it in full (e.g. `TRBV7-9*03`).
- Gene names with a slash (e.g. `TRAV29/DV5`) are fine to type as-is.

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| "Tkinter is not installed" | Run `brew install python-tk`, then relaunch. |
| "Stitchr failed" | Run `pip install stitchr` and `stitchrdl -s human`. |
| "No species data detected" | Run `stitchrdl -s human` (Stitchr's gene data isn't downloaded yet). |
| Seam check says CHECK | Open the document and verify the beta junction shows an **E** before `DLRNV`. Re-check the CDR3 and J gene. |
| Alpha or beta CDR3 not found | A gene or CDR3 was likely mistyped — double-check the six values. |

---

## How the insert is built

The generated insert has this structure (5′ → 3′), with the region between the brackets being what gets synthesised:

```
[ NotI — signal peptide — BETA chain — muTrbc — P2A — signal peptide — ALPHA chain — muTrac start ] — rest of muTrac — SpeI
```

The alpha and beta chains use the **human variable region** (V/CDR3/J) paired with **mouse constant regions**, delivered via a lentiviral vector.

---

## Credits

Built on [Stitchr](https://github.com/JamieHeather/stitchr) by James Heather for TCR sequence reconstruction. Insert-generation logic adapted from the Benaroya Research Institute TCR lentiviral pipeline.
