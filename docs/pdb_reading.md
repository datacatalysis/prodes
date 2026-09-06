# How prodes reads a PDB file

From version 7.1 the coordinate records of a PDB file are read by `Bio.PDB.PDBParser`, and not by prodes' own column arithmetic. Biopython is used as a **record reader**, not as an object model: `PDBParser` is driven with a `StructureBuilder` of prodes' own, `prodes.io.pdb_reader.PdbRecordBuilder`, which collects one flat `AtomRecord` per coordinate line instead of assembling Biopython's `Structure > Model > Chain > Residue > Atom` tree.

That last sentence is the whole of this document. Everything below says why, because it is not the obvious choice and the obvious one does not work.

## The shape of the parse

```
file text
   │
   ├─ prodes.io.pdb_reader.read_atom_records   Bio.PDB.PDBParser + PdbRecordBuilder
   │      one AtomRecord per ATOM/HETATM line, in file order
   │
   ├─ prodes.io.parser.check_records_were_all_read
   │
   ├─ prodes.io.parser.keep_first_model        one model per structure, from 8.0
   │
   ├─ prodes.io.conformers.elect_conformers    one conformation per residue, the 7.0 rules
   │
   └─ prodes.io.parser.build_structure         Structure, Chain, Residue, Atom
```

`AtomRecord` is the piece that matters. It carries every field separately: the atom name, the alternate location, the residue name, the chain, the residue number, **the insertion code**, the coordinates, the occupancy, the segment id, the element, **the model number**, and whether the line was an `ATOM` or a `HETATM`. The reader that came before it had no insertion code on an atom and no notion of a model at all, which is why the two parser defects open at 7.1 — merging residues that differ only by an insertion code, and merging the models of an NMR ensemble — were parser rewrites rather than changes to a grouping key. Both were fixed in version 8.0 as exactly that: the insertion code joined the residue key in `build_structure`, and the model number became the filter `keep_first_model` applies before the election. See `docs/residue_identity.md`.

## Why not Biopython's entity tree

Biopython's `Structure` tree is the natural thing to build on, and it was tried first. Four things stop it, three of which are only visible if you run it against the cases prodes' own test suite already pins.

**It silently drops residues prodes reads today.** `StructureBuilder.init_residue` raises when two conformations of one residue share an atom written with a blank alternate location. In the permissive mode everybody uses, `PDBParser` catches that, sets the current residue to `None`, and every following atom of that residue is discarded — with a warning and no error.

**An atom with no alternate location that arrives after a disordered residue has formed is discarded** as "defined twice". That atom is exactly the one the 7.0 rules always keep: an alternate location marks a conformation, and an atom without one belongs to every conformation. It is the reason the ordinary case, an ordered backbone with a disordered side chain, comes out whole.

**A disordered atom whose occupancy column is blank crashes.** `DisorderedAtom.disordered_add` compares the new occupancy against the best seen so far, and a blank column is `None`. prodes treats a blank occupancy as a real case — zero is a real occupancy, and a file that stops before column 60 has not said the atom is absent — so this is not an obscure corner for it.

**Coordinates are `float32`.** `PDBParser` stores `np.array((x, y, z), "f")`, which perturbs every coordinate in the seventh significant digit.

Biopython's disorder handling is also weaker than the rules prodes settled on in version 7.0, which is a separate objection and the one that is easiest to state: `DisorderedAtom` chooses per atom by highest occupancy, which assembles one residue out of two rotamers, and `DisorderedResidue` chooses by record order, which never consults occupancy at all. See `docs/alternate_conformations.md`.

Collecting records avoids all of it, because no `DisorderedAtom` and no `DisorderedResidue` is ever constructed. Nothing is dropped, nothing is relabelled, nothing crashes on a blank column, and prodes keeps every decision about identity and disorder that it had before.

## Coordinates

The builder stores `round(float(value), 3)`, which recovers the file's own number **exactly** rather than approximately. A PDB coordinate column is eight characters holding three decimals, so the largest value it can carry is 9999.999; over that range the spacing of `float32` is at most 0.000977, half of which is 0.000488, which is below the 0.0005 that would round to a different three-decimal value. Checked over the whole three-decimal grid the format can express.

A file that writes a fourth decimal into that column loses it, and one that writes a value the column cannot hold is shifted. Neither is silent: `prodes.io.parser.warn_about_unwritable_coordinates` says so once per file.

## What Biopython does that prodes now relies on

- The column arithmetic for every coordinate record, including the `MODEL`/`ENDMDL` bookkeeping and the skipping of `ANISOU`, `SIGATM`, `SIGUIJ` and `TER`.
- The hetero flag, which is how an `ATOM` record is told from a `HETATM` one.
- Permissive handling of unreadable occupancy and temperature-factor columns.

It is constructed `PDBParser(PERMISSIVE=True, QUIET=False)` with the complaints captured and re-emitted through prodes' logger. `PERMISSIVE=False` would make a blank occupancy fatal; `QUIET=True` would throw the complaints away, and at the default they would go to stderr, two per atom, for a file trimmed at column 54.

## What Biopython does not do

- **It does not parse `SSBOND`.** Its header stops at the first coordinate record and carries the title, the resolution and the journal. `prodes.io.parser.ssbond_records` reads those lines from the text, as prodes always has, from anywhere in the file. A lost `SSBOND` is a charge change and not a cosmetic one; see `docs/disulfide_bond_detection.md`.
- **It does not read past an `END` or `CONECT` record.** Everything after the first one is returned as an unparsed trailer, so a file that writes coordinates after either would lose them. `check_records_were_all_read` compares the number of coordinate lines in the text against the number of records and raises, naming the cause, rather than describing a structure that is missing atoms. Note that only the six-column `END   ` of the specification counts: a bare three-character `END`, which prodes' own `write_pdb` produces, is not recognised and truncates nothing.
- **It does not infer a blank element column; `pdb_reader.infer_element` does.** Biopython's `Atom` class guesses the element from the atom name, but no `Bio.PDB.Atom` is ever constructed here, so calling it would mean building one per record just for this. `infer_element` follows the same rule directly instead: the first character of `fullname` (columns 13-16 verbatim) tells an inorganic element written from column 13, `"FE  "`, from an organic one left blank, `" CA "`, and a name starting with a digit, a stereo-specific hydrogen such as `1HB1`, is read from its second character. It only runs where the column is blank; a column that holds something, however unlikely, is kept as written. Left uninferred, a blank element reaches `data.vdw_radius("")` from the surface calculation and raises, so a file with no element column produced no features at all before this. The count of atoms inferred is logged once per file, over every record the file holds, and reaches `prodes_run.json` as `inferred_elements`, because a structure whose radii come from guesses is a different thing from one whose radii were read. The two counts are not always the same number: the log is taken at read time, before `keep_first_model`, `elect_conformers` and the `ATOM`/`HETATM` filter have thrown anything away, while `inferred_elements` is taken over only the records that survived all three and became part of the built structure. They agree for a single model file read for one record type with no alternate conformations, which is most files; a multi-model ensemble, or a file whose inferred atoms are on a losing alternate conformer or the record type not requested, logs a larger count than it records.
- **It does not write PDB files.** `write_pdb` still writes the columns itself, because it also writes the dummy surface-point records that `PDBIO` has no way to express.
- **Its header parser is not used.** `read_atom_records` blanks every line before the first coordinate record rather than handing it over, because `Bio.PDB.parse_pdb_header` looks for a date anywhere in a `HEADER` or `REVDAT` record and looks the month up in a list of English abbreviations: a structure stamped `02-Mai-99`, or any header text matching the same pattern, raises `ValueError: 'Mai' is not in list` and takes the whole parse with it. Blanked and not dropped, so the line numbers the reader quotes are still the file's own.

## Everything that changed when the reader changed

No feature value moves on any well-formed structure: the eight structures in `tests/data` come out identical atom for atom, coordinates as exact floats, and identical in all 106 feature columns. What follows is the complete list of what *does* differ, all of it on input outside the format's contract.

| what | before | now |
|---|---|---|
| a file holding no records of the requested type | `IndexError: index -1 is out of bounds` | `ValueError` naming the file and the record type |
| an unreadable coordinate, residue number or insertion code column | `ValueError`, or nothing at all for a malformed `ANISOU` | `ValueError` naming the file; a malformed `ANISOU` now stops the parse, because Biopython reads those columns before the builder is reached |
| an element symbol written in lower or mixed case | kept as written, `Fe` | uppercased, `FE`. The format specifies right-justified uppercase, and no shipped structure writes anything else |
| two spellings of one residue number, `  30` and `0030` | two separate conformer elections | one election, because the residue number is now compared as the integer the format says it is |
| a coordinate record after a six-column `END   ` or any `CONECT` | read, though the reader that produced the file may not have meant it to be | `ValueError`, because the records are invisible to the reader and losing them quietly is worse |
| a coordinate written to more than three decimals | kept | read to three decimals, with one warning per file |
| a file holding more than one `MODEL` | merged into one impossible structure, silently | merged the same way, with a warning saying so. From version 8.0 the first model is described and the rest are dropped, which is where that warning went |
| a residue key in the alternate-conformer report | the raw columns, so a blank chain gave a double space and `0022` stayed padded | the parsed chain and number, so `PRO 22` and `PRO A22` |
| a record name padded with something other than spaces, such as `ATOM\t\t` | read, because the name was compared stripped | not read, and the reader says it ignored an unrecognised record |
| a negative occupancy | silent | one log line, from the reader. The conformer chosen is unchanged |
| a coordinate column holding `nan` or `inf` | read as `nan`/`inf`, silently | read the same way, with one log line saying every distance from those atoms is meaningless |
| a blank element column | `KeyError` from the van der Waals radius lookup, naming neither the file nor the atom; every feature was lost | inferred from the atom name, the same rule `Bio.PDB.Atom` uses, with the count logged once per file and recorded as `inferred_elements` |

The alternate-conformer report row and the blank-element row both reach `prodes_run.json`. Nothing in the list changes a number for a file that conforms to the format.

One asymmetry is worth naming because it is the only case where the change turns a working number into a crash rather than the other way round: `Cr` is prodes' own key for an aromatic carbon in `data.vdw_radius`, so a file that literally writes `Cr` in the element column used to get a radius and now raises. Chromium is `CR` in the format and no real file writes `Cr`, but the asymmetry is real.

## What follows from this

`MMCIFParser` drives the same `StructureBuilder` interface and calls no method `PdbRecordBuilder` does not implement, so reading mmCIF is a second entry point onto the same `AtomRecord`s rather than a second parser.
