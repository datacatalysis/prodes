# What Prodes calls one residue

A PDB file names a residue by three things: the chain, the residue sequence number, and the insertion code in column 27. Prodes read the first two and ignored the third, so two residues that differed only by their insertion code were read as one residue carrying both sets of atoms. Separately, it had no notion of a `MODEL` record, so an NMR ensemble holding the same molecule 20 times was read as one structure holding all 20 copies at once.

From version 8.0 a residue is `(chain, number, insertion code)`, and a file holding several models is described by its first.

The two are one change because the damage was one mechanism. Prodes could not tell two residues apart when only the insertion code separated them, and it could not tell one model's copy of a residue from the next model's, so in both cases the surplus atoms landed on whichever residue happened to have been made most recently. What differs is where the surplus came from: from residues the file genuinely holds, or from further copies of the ones it holds already.

## Insertion codes

An insertion code lets a depositor add a residue between numbers 100 and 101 without renumbering the chain after it. It is not a corner of the format for the structures Prodes is aimed at: Kabat and Chothia antibody numbering write **every** CDR insertion this way, so `H100A`, `H100B` and so on are the normal appearance of an antibody heavy chain.

### What the merge cost

4NZU, a 1.20 Å Fab, is the shipped example. Three of its residue numbers carry insertion codes, and between them 12 residues were lost:

| read as | atoms | distinct atom names | really |
|---|---|---|---|
| `H52` SER | 30 | 21 | `52` SER, `52A` TYR |
| `H82` MET | 49 | 20 | `82` MET, `82A` ASN, `82B` SER, `82C` LEU |
| `H100` ASP | 111 | 35 | `100` ASP, `100A` ALA, `100B` ASP, `100C` CYS, `100D` HIS, `100E` LYS, `100F` GLY, `100G` ALA, `100H` PHE |

`H100` was nine residues read as one. The file holds 434 residues and Prodes reported 422.

**The mass was wrong.** `Residue.mass` looks up the merged residue's single name, so the other eleven residues contributed nothing: the residues summed to 45 084 Da against the file's 46 392 Da, an error of 1.3 kDa or 2.8 per cent. Every residue-composition and surface-fraction column was attributed to the wrong residue type along with it.

**The charge was wrong, by whole units.** An atom keeps the residue name from its own record, but takes its pKa from the `Residue` object it was merged into, so an atom was recognised as chargeable and then titrated against another residue's group. Inside `H100`, read as an aspartate:

| atom | pKa used | charge at pH 7 | should be |
|---|---|---|---|
| `LYS NZ` | ASP, 3.86 | 0.00 | +1 |
| `CYS SG` | ASP, 3.86 | **−1.00** | ~0 |
| `TYR OH` (in `H52`) | SER, none | 0.00 | ~0 |

A free thiol has a pKa of 8.33 and is neutral at pH 7; titrated against aspartate's 3.86 it contributed a full −1.

**It destroyed a disulfide bond, reintroducing the defect version 6.0 fixed.** `cysteine_sulfurs` skips any residue not named `CYS`. `H100C` is a cysteine, but it had been merged into a residue named `ASP`, so it was not a cysteine as far as detection was concerned, and 4NZU's fourth `SSBOND` record was discarded:

```
an SSBOND record joins H98 to H100, which are not both cysteines in this file; ignoring it
```

Four of five disulfides were found, nine cysteines were seen where the file has ten, and both halves of the missing bond titrated as free thiols. See `docs/disulfide_bond_detection.md`.

### What changed

- `build_structure` keys residues on `(chain, residue number, insertion code)`. `Residue` and `Atom` both carry the code.
- `Residue.label` formats a residue as it is named in a warning: `A22`, or `H100C`. Every message that names a residue reads it from there, because a message saying `H100` would be describing nine of them.
- `read_ssbond_line` reads columns 22 and 36, the insertion codes of an `SSBOND` record, and disulfide records are resolved on the full identity. The codes are read by slice rather than by index, so a record trimmed after the residue number keeps its bond rather than losing it to a missing column.
- `write_pdb` writes column 27, on `ATOM`, `HETATM` and `TER` records alike. A writer that dropped it would hand back a file whose two residues are one again.

### What 4NZU now reports

| | before | now |
|---|---|---|
| residues | 422 | 434 |
| molecular weight | 45 102 Da | 46 410 Da |
| formal charge at pH 7 | −4 | −2 |
| disulfide bonds | 4 | 5 |

Measured over the whole pipeline on 4NZU, at `--full-features --ionic-strength 0`, 69 of the 105 feature columns move and 36 do not. Nothing geometric moves: no atom moves, and no atom's own residue name changes, so `Area`, `NSurfPoints`, both `Shape` columns and every hydrophobicity column are identical. What moves is everything that reads a residue's type, mass or charge.

No structure without an insertion code moves at all: the other seven structures in `tests/data` parse to the same atoms, the same residues, the same molecular weight and the same charge, atom for atom and coordinate for coordinate.

## Models

An NMR structure is an ensemble. `MODEL` and `ENDMDL` wrap each member, and the file holds the same molecule many times over at the same chain names and residue numbers.

Prodes read every one of them into a single structure. Because a residue key that has already been seen on a chain does not start a new residue, every model after the first had its atoms appended to whichever residue was created most recently. On 1PIT, bovine pancreatic trypsin inhibitor, 20 models of 58 residues:

| | read whole, before | first model, now |
|---|---|---|
| atoms | 17 780 | 889 |
| residues | 58 | 58 |
| largest residue | 16 901 atoms | 24 atoms |
| formal charge at pH 7 | −1096 | +6 |

The real charge of BPTI at pH 7 is about +6. Nothing warned about any of this until version 7.1, and what it warned was that the answer was not usable.

Issue #13 reports the charge as -1362. That figure was taken before the disulfide work of version 6.0 and the alternate-conformation work of 7.0, both of which moved it; the version this change replaces gives -1096 on the same file.

`tests/data/1PIT.pdb.zip` is the RCSB deposition unchanged, all 20 models and all 17 780 coordinate records. It is not trimmed the way `1AO6.pdb.zip` was, because its models are the point of it and its header is 1 KB of 1.4 MB.

### What changed

`keep_first_model` in `prodes.io.parser` keeps the records of the first model and logs how many the file held and how many records were dropped. The structure that comes out is identical, atom for atom and coordinate for coordinate, to parsing a file holding that model alone.

"First" is the lowest model number among the records of the type being built, rather than model 1 outright, so a file whose `HETATM` records begin in a later model is described rather than refused for holding no records at all.

The selection happens **before** the conformer election, which settles two things that would otherwise be decided across models. The election keys residues on chain, number and insertion code and carries no model, so one model's alternate conformations could decide another model's residue; and the residue grouping is by key already seen, so a repeated residue number is what caused the pile-up above. Both become unreachable once one model is chosen.

`prodes_run.json` records `models_in_file`. The bundle ships the input file unchanged, so a bundle whose features describe model 1 of 20 is otherwise indistinguishable from one describing a crystal structure.

### What this is not

It is not conformational averaging. Prodes describes one member of the ensemble, chosen by position in the file rather than by anything about the structure, which is what most tools do with an NMR file and what a user handing one over almost always means. Averaging descriptors over the members is a different piece of work and would answer a different question; see the README's comparison section on what the commercial packages do instead.

## What is still narrower than a PDB file

**A pKa file names a residue by chain and number, but not by insertion code.** From version 9.0 `Structure.redo_pkas` and the converters in `prodes.io.pka_converter` key on `(chain, residue number)`, so a value predicted for one chain is no longer offered to another chain's same-numbered residue (see the README's [pKa values and protonation states](../README.rst) section). Insertion codes are a different, still-open gap: a chain and number alone do not distinguish two residues that differ only by insertion code.

The predictor does the same thing on its way in. PROPKA 3.5.1 run on 4NZU writes five different residues into its summary as number 100 of chain H, the insertion code nowhere in the line:

```
   ASP 100 H     1.76       3.80        (H100)
   ASP 100 H     1.76       3.80        (H100B)
   HIS 100 H     3.31       6.50        (H100D)
   CYS 100 H    49.99       9.00        (H100C)
   LYS 100 H     5.81      10.50        (H100E)
```

so `convert_propka` keys all five under `(H, 100)` and nothing in the file says which residue each belongs to. What saves it is that entries are applied by **group** name: the histidine value lands on the histidine, the lysine value on the lysine, and the cysteine's 99.99-style marker is dropped because that cysteine is in a disulfide and does not titrate, which PROPKA agrees with. Measured on this file, every value lands on the residue it was predicted for.

Two residues of the same type sharing a number in one chain are the case that stays wrong: `H100` and `H100B` are both aspartates, and both take 1.76 whichever of them it was predicted for. Before chains were tracked at all, the position was worse rather than better, since all nine insertions were one residue named `ASP` and four of the five predictions were dropped for naming a group it did not have.

Fixing this properly means adding an insertion code to the pKa file format too, a further change to the same documented public format.

**A residue written in two places is still merged.** The grouping test is membership in the keys already seen on the chain, so a file that writes one residue's atoms, then another residue, then the first residue again puts the late atoms on whichever residue was made most recently. Selecting one model put the case that mattered out of reach rather than changing this rule.

## Where the code is

`prodes.io.parser`: `build_structure` groups residues, `keep_first_model` selects the model, `read_ssbond_line` reads the record columns and `write_pdb` writes them back. `prodes.core.residue.Residue.label` is how a residue is named in a message. The insertion code and the model number each arrive as their own field on an `AtomRecord`, read by `prodes.io.pdb_reader`; see `docs/pdb_reading.md` for why they were not there before version 7.1.
