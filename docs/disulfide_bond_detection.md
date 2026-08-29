# Cysteines and disulfide bonds

How Prodes decides which cysteines are bonded into a disulfide, why it matters to the charge features, and what it still gets wrong. Behaviour from version 6.0.

## Why it matters

A cysteine whose `SG` is bonded to another cysteine's has no thiol proton, so it does not titrate at all. Prodes finds those bonds and gives the two cysteines no side-chain pKa, which keeps them neutral at every pH.

Before version 6.0 it had no concept of a disulfide and gave every `CYS` the free-thiol pKa of 8.33. On lysozyme, which has four disulfides and no free cysteine, that put the formal charge at pH 8.5 at **-1 when it is +7**, and the isoelectric point at 8.9 when it is 10.3.

How much of your output moves depends on the pH you run at, but it is never nothing:

* **At pH 7**, one of the 54 default features changes on lysozyme: the isoelectric point, which is a property of the whole titration curve and so shifts at any pH you ask for it. With `--full-features` 20 of the 105 change, because the `*Average` columns use fractional charges and a free thiol at pKa 8.33 is already 4 per cent ionised at pH 7.
* **At pH 8.5**, 21 of the 54 change, including the formal charge, the dipole and the whole `SurfEp` block.

Because this changes the value of features that were already being reported, feature tables calculated with an earlier version are not comparable, and a model fitted on them needs refitting. That is why it is a major version bump.

## How the bonds are found

* An `SSBOND` record is authoritative **for the two cysteines it names**. Records that carry a crystallographic symmetry operator, that join a residue to itself, or that name a residue not in the coordinates are skipped, and the rest of the cysteines still go to the geometry. A record whose two sulfurs are not at bonding distance is honoured, with a warning, because that usually means a strained bond and occasionally means a reduced structure carrying stale records.
* Every cysteine no record claims is paired **by distance**, with a cutoff of **2.5 Å** between the `SG` atoms, each sulfur taking at most one partner and the shortest bond winning a contested one.

Records are applied per cysteine rather than per file, so a file whose record set is incomplete keeps the bonds its records forgot, and a file whose only records are crystallographic still gets its real bonds by geometry.

AlphaFold models carry no `SSBOND` records but do place the sulfurs at bonding distance, so they are handled by the geometric route and need nothing special.

### Why 2.5 Å

It is the same cutoff PROPKA and PDB2PQR use, so the distance criterion is the one your pKa predictor already applies. The two can still reach different answers, since Prodes also honours `SSBOND` records and PROPKA does not read them.

It is not tighter because real bonds run past 2.35 Å. The SARS-CoV-2 spike trimer, 6VXX, a 2.8 Å cryo-EM structure, models `Cys391`-`Cys525` at 2.361 Å in all three protomers, while every other bond in that entry sits between 2.018 and 2.058 Å.

It is not wider because there is no room above. The AlphaFold model of metallothionein-2, a protein with twenty metal-binding cysteines and no disulfides at all, places two of them 2.050 Å apart, with the rest of the metal cluster at 3.3 Å and beyond. Widening the cutoff toward those distances would turn a single false positive into a dozen.

A dihedral filter was considered and rejected. Some tools screen on the `CB-SG-SG-CB` angle at 90 ± 10°, but measured across 142 annotated disulfides that angle spans 37° to 174°, with 54 per cent of real bonds outside 80° to 100°. It would reject more real bonds than it saves.

## Checking that it worked

The count found is printed at the start of the run and recorded in `prodes_run.json` inside the output bundle. It is worth a glance: a structure you expect to have disulfides that reports none is being titrated as though every cysteine were free.

> **This matters most when you are not using PROPKA.** PROPKA detects disulfides itself and reports a bridged cysteine as `99.99`, its marker for a group that does not titrate, and Prodes has always passed that through. So a run with `--pka` was already close to right, and a run without it was not.

## What is still not handled

In every case because the evidence is not an `SG`-`SG` distance:

* **A cysteine bound to a metal**, as in a zinc finger, is coordinated and deprotonated rather than protonated, and Prodes titrates it as a free thiol. The metal is in a `HETATM` record, which the parser does not read.
* **A thioether link**, such as the two cysteines that bond to the haem of cytochrome c, is nowhere near `SG`-`SG` bonding distance.
* **A real bond a model has stretched.** A low-confidence AlphaFold region can place a genuine disulfide well past 2.5 Å, in which case both cysteines are titrated. This fails in the safe direction, back to the pre-6.0 behaviour.
* **A bond that is not there.** The reverse also happens, as the metallothionein model above shows. No distance cutoff can tell a 2.05 Å metal-cluster contact apart from a real bond.
* **An inter-chain bond in a single-chain model.** An antibody heavy chain modelled on its own has nothing to bond to, so the cysteines that would join the light chain look free.
* **A cysteine with alternate locations.** Prodes reads the `altLoc` column as part of the atom name, so a disordered `SG` is not recognised as an `SG` at all. Such a cysteine is invisible both to this detection and to the charge calculation, which means it carries no charge either way.

One thing deliberately unchanged: a cystine is more hydrophobic than a free thiol, but `CYS` keeps a single hydrophobicity value, so the `Mhp` features and the `CYSSurfFrac` column do not distinguish the two.

## Where the code is

`prodes.calculations.disulfides`, called by `PDBparser` on every structure it reads, so a structure that has been parsed from a file always knows which of its cysteines are bonded. A bonded cysteine is marked with a `disulfide_partner`, which stops `Residue.pkas` offering a titratable side chain; the charge correction follows from that rather than being special-cased anywhere.
