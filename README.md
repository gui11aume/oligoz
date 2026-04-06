# Oligoz

Oligoz is a single-file Python tool (`oligoz.py`) for PCR primer design.
It can be used:

1. as a command-line program to design primers from input sequences, and
2. as an importable Python module for sequence and Tm calculations.


## What the software does

Oligoz designs primer pairs using a **product-centric search**:

1. It builds candidate forward and reverse primers from a template sequence.

2. It validates each candidate primer independently (**consistency checks**):
   - primer length,
   - melting temperature (Tm),
   - self-annealing / self-dimer potential.
3. It validates primer pairs together (**compatibility checks**):
   - difference between primer Tm values,
   - hetero-dimer and non-specific annealing between primers.
4. It ranks valid pairs and prints the best one for each target sequence.

Tm is calculated from nearest-neighbor thermodynamics with salt corrections.
Internally, temperatures are handled in Kelvin and displayed in Celsius.


## Default design rules

Standard mode defaults:

- `min_length`: 17
- `max_length`: 30
- `min_Tm`: 58 C
- `max_Tm`: 65 C
- `max_self_anneal`: 3
- `max_delta_Tm`: 3
- `max_hetero_anneal`: 3
- `max_lanneal`: 4

Default hybridization conditions:

- oligo concentration (`--conc`): 1 uM
- monovalent ions (`--mono`): 50 mM
- divalent ions (`--di`): 1.5 mM


## Input format

You can provide input as:

- FASTA (multi-entry supported), or
- one raw DNA sequence per line (no header).

Non-DNA letters (outside DNA/degenerate handling expected by the code path)
can cause a sequence to be skipped during search.


## Installation (get the code)

Oligoz is a single script, so installation means getting the repository
onto your machine.

### Option A: Clone with Git

Using HTTPS:

```bash
git clone https://github.com/gui11aume/oligoz.git
cd oligoz
```

Using SSH (if your GitHub SSH keys are already configured):

```bash
git clone git@github.com:gui11aume/oligoz
cd oligoz
```

### Option B: Download without Git

1. Open [https://github.com/gui11aume/oligoz](https://github.com/gui11aume/oligoz).
2. Click **Code** -> **Download ZIP**.
3. Extract the ZIP archive.
4. Open a terminal (or PowerShell) in the extracted `oligoz` folder.

After either option, continue with the Python run instructions below.


## How to run on Linux, macOS, and Windows

Oligoz does not require installation as a package. You run the script directly.

### 1) Check Python is available

Linux/macOS:

```bash
python3 --version
```

Windows (PowerShell or cmd):

```powershell
py -3 --version
```

If Python is missing, install Python 3 from [python.org](https://www.python.org/downloads/).

### 2) Run Oligoz

From the folder containing `oligoz.py`:

Linux/macOS:

```bash
python3 oligoz.py sequences.fasta
```

Windows:

```powershell
py -3 oligoz.py sequences.fasta
```

You can also pipe FASTA text directly through standard input.


## Command line usage

```text
oligoz.py [options] [infile]
```

Main options:

- `--qPCR` : qPCR-oriented search (targets between 110 and 180 bp with stricter defaults)
- `--force` : if no hit is found, retry with relaxed rules
- `--approx d` : allow approximate product boundaries by trimming up to `d` nt from ends
- `--extraL L` : prepend extra nucleotides to left primer
- `--extraR R` : prepend extra nucleotides to right primer
- `--conc C` : oligo concentration in uM
- `--mono C` : monovalent ion concentration in mM
- `--di C` : divalent ion concentration in mM
- Rule overrides: `--min_length`, `--max_length`, `--min_Tm`, `--max_Tm`,
  `--max_self_anneal`, `--max_lanneal`, `--max_delta_Tm`, `--max_het_anneal`

Show full help:

Linux/macOS:

```bash
python3 oligoz.py --help
```

Windows:

```powershell
py -3 oligoz.py --help
```


## Example 1: Standard PCR search

Input (`example1.fasta`):

```text
>example1
TGGTAGTAGCATCTGACCCGAAGTCTGTTGCTATTACACGGTTGAGCTATGAGCGAGTTATCCAGTTGGGGTCAAGCTTAGCACGGAGCCAGTGAATAATTCATCATACGAAAATGATGCCCCCGACAAAAAGAGGAGCCCACTTTCGTGGATAAACAGTGTCTTTTGTATACAGAGATCGGTTAGATGCCGGTGGCTAG
```

Run:

Linux/macOS:

```bash
python3 oligoz.py example1.fasta
```

Windows:

```powershell
py -3 oligoz.py example1.fasta
```

Output:

```text
example1

(63.1 deg) TGGTAGTAGCATCTGACCC
(63.5 deg) CTAGCCACCGGCATCTA
---
```


## Example 2: qPCR search

Input (`example_qpcr.fasta`):

```text
>example_qpcr
GTCAAAGATAACCATACAATACATATGCTAGTATGAAATCCGCCGGTTTAGCGCGAACTGTCCAAGGGCCAAGCAAGCGCGACGTACATGTTCCCATCCGTCGCGCTTATTTTACTGGATCGGTACCCCCACCATATCTAGAAATAGAATCTGGGCAACCCCGATAGGCTATGTAGAGGTGTGTTTCTTCGAAGCGTGCGGATATCTGTCACGAACTCTCGACCCATTTATCTCGTTAAATCTTAGTGTGGCTAGCCTTACTATTTCAAGCAATTGAACAACGTGCCTGTCTCACCGATCATGATGTGTCTACCTTTCCG
```

Run:

Linux/macOS:

```bash
python3 oligoz.py --qPCR example_qpcr.fasta
```

Windows:

```powershell
py -3 oligoz.py --qPCR example_qpcr.fasta
```

Output:

```text
example_qpcr

(60.6 deg) CGCTTATTTTACTGGATCGG
(60.1 deg) AGATAAATGGGTCGAGAGTT
---
```


## Example 3: Force a result with relaxed rules

If strict defaults find no pair, use `--force`:

Linux/macOS:

```bash
python3 oligoz.py --force difficult_targets.fasta
```

Windows:

```powershell
py -3 oligoz.py --force difficult_targets.fasta
```

When a sequence only yields a relaxed-rule result, the title is marked
with `(forced)` in the output.


## Using as a Python module

You can import `oligoz.py` and use its classes directly:

```python
import oligoz

sol = oligoz.OligoSol("GATCGATCGATCGATCGTAC")
tm_kelvin = sol.Tm
tm_celsius = tm_kelvin + oligoz.ABS0

print(tm_kelvin)
print(tm_celsius)
```

Key classes:

- `DNAseq`, `RNAseq`: sequence helpers and reverse-complement support
- `Oligo`: oligo sequence object with annealing-related methods
- `OligoSol`: oligo plus solution conditions; computes Tm and occupancy


## Notes

- The script prints the single best pair per input target.
- For very difficult templates, `--force` may still report `** pathological **`.
- Input and rule tuning matter; if needed, relax constraints progressively.
