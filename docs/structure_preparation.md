# Preparing a structure before running Prodes

Prodes describes the structure you give it. It does not correct one. If a surface lysine has no side chain in the file, Prodes gives that residue no charge and says nothing about it, and every charge-derived feature comes out wrong.

[PDB2PQR](https://github.com/Electrostatics/pdb2pqr) fixes most of that. It is free, BSD-3-Clause, pure Python, and a single command. This page is why you should run it, what it changes, and the traps.

## Why it matters

PDB entry **4HKZ** is a Fab. Two of its lysines, A188 and B221, are modelled only as far as CB, which is ordinary for surface residues at that resolution. Prodes reads the file without comment.

Running PDB2PQR first rebuilds those two side chains. The only difference between these two columns is whether those two side chains are present:

| | side chains unmodelled | side chains rebuilt |
| --- | --- | --- |
| Formal charge at pH 7.4 | **-2** | **0** |
| Isoelectric point | **6.46** | **7.61** |
| Mean surface electrostatic potential | -0.013 V | **+0.011 V** |
| Positive surface points | 12,976 | 13,684 |

**31 of the 54 default features move.** The net charge is wrong by 2, the isoelectric point is wrong by 1.14 units, and the mean surface potential has the wrong sign.

For choosing between an anion and a cation exchanger at a given buffer pH, an isoelectric point of 6.46 and one of 7.61 are different recommendations.

This is not universal. Many deposited structures are complete, and 1AO6, 1GPB, 4NZU, 1GDW, 1CBN, 1IGT and 1BEY all are. It is a property of resolution and of surface disorder, and the residues most often left unmodelled are the long charged ones, which are exactly the residues that carry the charge.

## What PDB2PQR does that Prodes does not

| PDB2PQR | Prodes on its own |
| --- | --- |
| **Rebuilds missing heavy atoms** from the residue's reference geometry | nothing; the atom stays absent and everything computed from it is wrong |
| **Detects backbone chain breaks** and names them | nothing; a nine residue gap in 1BEY passes without a message |
| **Refuses to continue** when a gap is too large to rebuild, naming the residue | nothing |
| **Deletes atoms that do not belong to a residue**, logging each one | keeps them, and assigns a radius from the first letter of the atom name |
| **Adds every hydrogen**, including polar hydrogens whose position depends on the hydrogen bond network | reads hydrogens if present, then filters them out; every calculation is heavy-atom only |
| **Optimises the hydrogen bond network**, including flipping Asn, Gln and His where the deposited assignment is worse than the flip | nothing |
| **Resolves steric clashes** created by rebuilt atoms, by rotating side chain dihedrals | nothing |
| **Adds a missing `OXT`** and treats termini as proper patches, with neutral variants available | assigns the terminus but does not add `OXT`, and places the C-terminal charge on the backbone carbonyl carbon |
| **Assigns a partial charge and a radius to every atom** from a named force field | formal charges on ionisable groups only, and one radius per element from a five-entry table |
| **Reads 57 PDB record types**, including `MODEL`, `MODRES` and `LINK` | reads coordinates and `SSBOND` |
| **Fetches a structure by PDB ID** | reads a local `.pdb` or `.pdb.zip` |

Prodes keeps hold of two of these itself, and does them well: it **recognises disulfide bonds**, so a cystine is not titrated as a free thiol, and it **resolves alternate conformations** by keeping the best occupied one per residue.

## What Prodes actually uses from the prepared file

Only the **geometry**. Prodes reads `ATOM` records, filters out hydrogens, and re-titrates every residue from its own pKa table at the pH you ask for.

So the hydrogens PDB2PQR adds, the protonation states it assigns and the partial charges it writes are all discarded. That is deliberate: it is what lets you run Prodes at several pH values against one prepared structure. What you gain from the preparation step is that the atoms are there and in sensible positions.

## What happens to a predicted structure: the terminal oxygen

Run PDB2PQR over a set of predicted structures and you will usually see it add exactly one atom to each, and nothing else. That atom is `OXT`, and it is worth understanding because of what Prodes currently does with it.

**What `OXT` is.** The C terminus of a protein is a carboxylate, `-C(=O)-O`. The PDB format names the first oxygen `O`, which every residue has, and the second one `OXT`, which only the last residue of a chain has. It is an ordinary heavy atom of the real molecule, not an artefact.

**Structure sources disagree about whether to write it.** Measured over the structures at hand:

| source | structures with `OXT` |
| --- | --- |
| AlphaFold database | 200 of 200 sampled |
| Boltz2 predicted monomers | 0 of 200 sampled |
| crystal structures | 21 of 29 |

Boltz2 leaves it off, AlphaFold always writes it, and a crystal structure has it or not depending on whether the terminal residue was resolved. So the same protein arrives with a different atom count depending on where its structure came from.

**Why that matters here.** Prodes assigns hydrophobicity per residue, giving every heavy atom the value of the residue it belongs to. `OXT` is the one exception: `Property_point.set_lipo` gives it a fixed value of 1.0, which on the default `mj_scaled` scale is the maximum, the same as phenylalanine. A carboxylate oxygen is one of the most hydrophilic groups in a protein, so the value is the wrong way round.

This behaviour is inherited from the original [tneijenhuis/prodes](https://github.com/tneijenhuis/prodes) repository, which this fork preserves.

**The practical consequence is inconsistency, and that is what the repair fixes.** Adding one `OXT` to a Boltz2 model moves a median of 9 of the 54 default features. Most of that is the ordinary effect of one more atom on the surface, but two hydrophobicity features move in one direction: `NSurfPosMhp` rises on 33 of the 39 structures where it moves at all, and `SurfMhpMean` rises on 12 of the 13.

The size of that shift is small. The problem is that today it is applied to every AlphaFold structure, no Boltz2 structure, and most but not all crystal structures. A model trained across sources is then partly learning which predictor produced each file.

Repairing every structure with PDB2PQR removes that, because every structure ends up with exactly one `OXT` per chain. The artefact becomes a constant rather than a per-source bias, which is the thing that actually damages a model. Until the constant itself is changed, consistency is what matters, and this is the cheapest way to get it.

## The commands

Install PDB2PQR into an environment of its own, not alongside Prodes:

```
conda create -n pdb2pqr -c conda-forge pdb2pqr      # once
```

Prodes never imports PDB2PQR. It only reads the repaired PDB file that PDB2PQR writes, so the two never have to be importable at the same time, and installing them together buys nothing. What it costs is that each one's pins constrain the other's for as long as both are installed. They are independently maintained projects on separate release cycles, and the first time one of them moves a shared dependency the other has not caught up with, a single environment stops solving and takes a working Prodes install down with it.

A `pip install pdb2pqr` into the Prodes environment usually works if you would rather have one environment. Note that the PyPI package pins `docutils<0.18`, which collides with Sphinx and several other common packages; the conda-forge build does not carry that pin.

Then:

```
mkdir prepared

conda activate pdb2pqr
pdb2pqr --ff=PARSE --keep-chain --pdb-output=prepared/1GDW.pdb 1GDW.pdb prepared/1GDW.pqr

conda activate prodes
propka3 prepared/1GDW.pdb                    # the prepared file, not the download
python -m prodes.io.pka_converter 1GDW.pka propka -o 1GDW_pka.json
python -m prodes prepared/1GDW.pdb 1GDW.zip --ph 7.4 --pka 1GDW_pka.json
```

Inside a script, `conda run -n pdb2pqr pdb2pqr ...` calls PDB2PQR in its own environment without switching by hand, so a loop can run start to finish with the Prodes environment active.

`--pdb-output` is the flag that matters. The `.pqr` file is PDB2PQR's normal output and Prodes cannot read it; `--pdb-output` writes the prepared structure as a PDB, which Prodes can.

**Keep the original file name and change the directory.** Prodes takes the `ID` column from the file name, so a prepared file called `1GDW_prep.pdb` labels that row `1GDW_prep`. Over a dataset that is tedious to undo.

**`--keep-chain` is in the command on purpose.** It does not change the file Prodes reads, but it does change PDB2PQR's other output. See the section below.

**PROPKA writes its `.pka` into the directory you are standing in**, named after the input's base name, not next to the input file. That is why the third line above reads `1GDW.pka` and not `prepared/1GDW.pka`.

Prodes will print `Ignoring unrecognized record 'TER'` on a PDB2PQR output file. That is harmless: `TER` is a chain separator and carries no coordinates.

Install from **conda-forge**, not pip. The PyPI package pins `docutils<0.18`, which collides with Sphinx and with several other common packages. The conda-forge build does not carry that pin.

## Run PROPKA on the prepared file, not on the download

PDB2PQR can run PROPKA itself, with `--titration-state-method=propka`. **That does not replace the `propka3` step**, for two reasons:

* PDB2PQR does not write a `.pka` file. It uses the values internally to pick protonation states and then discards them, so `prodes.io.pka_converter` has nothing to read.
* PDB2PQR bakes in one pH. Prodes separates prediction from calculation so that you can predict once and then run at as many pH values as you like.

What you should take from PDB2PQR is the **ordering**. It runs PROPKA after repairing the structure, and that is the right way round. On 4HKZ:

| | titratable groups predicted |
| --- | --- |
| PROPKA on the deposited file | 145 |
| PROPKA on the prepared file | 149 |

The four extra groups are the two rebuilt lysine side chains and the two C-termini that only exist once `OXT` has been added. Among the 145 groups both runs share, six shift by more than 0.5 pKa units.

## Chain identifiers

Prodes needs chain identifiers. It groups residues by chain, decides which cysteines are bonded into a disulfide from the chain and residue number together, and applies per-residue pKa values by residue number. Flattening the chains would change all three.

PDB2PQR writes two files, and they behave differently.

| output | chain identifier |
| --- | --- |
| `--pdb-output`, the file Prodes reads | **always written**, with or without `--keep-chain` |
| the `.pqr`, PDB2PQR's normal output | **blank** unless `--keep-chain` is given |

The reason is in `io.print_biomolecule_atoms`: the PDB branch calls `atom.get_pdb_string()`, which takes no chain flag, while the PQR branch calls `atom.get_pqr_string(chainflag=chainflag)`.

Verified rather than assumed. On structures of two, four, six and eight chains, the `--pdb-output` file carries every chain identifier, and running PDB2PQR with and without `--keep-chain` produces byte-for-byte identical PDB output. On the six-chain structure 1F6R, Prodes reads the same six chains, the same 727 residues and the same 24 disulfide bonds from the repaired file as from the original.

So the flag is not required for this pipeline. It is in the documented command anyway, for two reasons: the `.pqr` is what APBS and the tools built on it consume, and nobody should have to remember which of two output files preserves what.

## Traps

**Do not use `--ffout`.** It renames residues into the force field's own scheme, putting `ASH`, `GLH`, `LYN` and `HID` into the file. Prodes raises `KeyError` on any of those. Without `--ffout` both PARSE and AMBER keep canonical residue names, verified at pH 4.0 and pH 7.4.

**PDB2PQR applies terminal patches at chain breaks.** It puts a real +1 and -1 at a position where the chain is actually continuous and merely unmodelled. Prodes assigns termini only at the first and last residue of a chain, which is the better answer, and it keeps doing so on the prepared file because it re-titrates from scratch. The charges PDB2PQR invented at the gap are discarded along with the rest of its charge assignment.

**PDB2PQR can fail on sequence microheterogeneity.** Crambin (1CBN) has residue 22 modelled as both PRO and SER. PDB2PQR 3.6.1 picks the highest occupancy per *atom* without keeping the residue identity consistent, ends up with a proline carrying a serine hydroxyl, and exits with `Unable to debump biomolecule`. Prodes handles this case correctly on its own, so if PDB2PQR refuses a structure, running Prodes on the original file is a reasonable fallback.

**mmCIF input does not work reliably.** PDB2PQR 3.6.1 accepts a `.cif` but produced an empty `.pqr` from a valid one during testing. Convert to PDB first.

**Be consistent within a dataset.** The same rule as for PROPKA: prepare all of your structures or none of them. Mixing prepared and unprepared structures puts two different kinds of number in the same feature column.

## If you would rather not add the dependency

The single most valuable thing PDB2PQR gives you here is knowing that a residue is incomplete. Until Prodes reports that itself, you can check without installing anything by comparing each residue's atom names against the standard set for its type, and treating any structure with missing side-chain atoms on charged residues as suspect.

## References

PDB2PQR asks to be cited as both of these, and prints both on every run:

Dolinsky TJ, Czodrowski P, Li H, Nielsen JE, Jensen JH, Klebe G, Baker NA. PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations. *Nucleic Acids Research* **35**(Web Server issue), W522-W525 (2007). https://doi.org/10.1093/nar/gkm276

Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, et al. Improvements to the APBS biomolecular solvation software suite. *Protein Science* **27**(1), 112-128 (2018). https://doi.org/10.1002/pro.3280
