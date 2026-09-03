# Alternate conformations

A crystal structure often models a residue in more than one position. Each conformation is written out in full, and every atom of it carries a letter in column 17, the `altLoc` indicator, together with an occupancy saying how much of the time that conformation is the one that is there.

Before version 7.0 Prodes read the atom name as columns 13 to 17, which is the name **plus** that letter. A disordered aspartate's oxygens arrived called `OD1A`, `OD1B`, `OD2A` and `OD2B`, and a disordered cysteine's sulfur as `SG A`.

From 7.0 the name and the indicator are read as the fields they are, and one conformation per residue is kept.

## What the old behaviour cost

A name with the letter attached matches nothing the rest of the package looks for, and the package looks up a great deal by atom name.

| what is looked up by name | consequence of the mangled name |
|---|---|
| the charged atoms of a residue | the side chain carries no charge at any pH |
| the terminus, matched as `N` or `C` | a disordered terminal residue loses its terminal charge |
| the aromatic carbons, which carry their own radius | an aromatic carbon gets the plain carbon radius, 2.0 Å instead of 1.85 Å |
| the `SG` a disulfide is found from | a disordered cysteine is invisible to the geometric search |

At the same time every conformation survived as a separate atom, so a disordered residue put two overlapping side chains into the surface calculation.

Nothing warned. The output was a plausible-looking CSV row.

On 1CBN, a 0.83 Å crambin structure where 13 of 46 residues are disordered, that came to 272 of 772 atoms parsed under a wrong name, 400 heavy atoms where the structure has 327, ASP 43 reading a formal charge of 0 where it is −1, and 6 of the 18 aromatic carbons getting the right radius.

## Which numbers move

Measured by running the pipeline before and after on three structures spanning the range:

| | 1CBN, 28% of residues disordered | 1F8N, 2.1%, six of them titratable | 7BZ5, 1.0%, none titratable |
|---|---|---|---|
| default features changed | 44 of 54 | 38 of 54 | 12 of 54 |
| Formal charge | 0.0 → 0.0 | −12.0 → −11.0 | 7.0 → 7.0 |
| Isoelectric point | 6.170 → 6.898 | 6.043 → 6.137 | unchanged |
| Dipole | 70.8 → 106.8 | +13% | +0.3% |
| Area | −3.0% | −0.17% | ~0 |

**The size of the change tracks how many *titratable* residues are disordered, not how many residues are disordered.** 7BZ5 has alternate conformations throughout and barely moves, because none of them is on a residue that carries charge. This is worth stating plainly because the intuition runs the other way: the duplicated side chains look like the problem, and the surface is in fact the smallest of the effects.

AlphaFold models contain no alternate conformations, so the standard AlphaFold workflow produces byte-identical output to 6.0.

## How the conformation is chosen

Per residue, not per atom.

The distinction is not academic. A residue can be modelled as two **different** amino acids, which is called microheterogeneity: in 1CBN residue 22 is a serine at occupancy 0.20 and a proline at 0.60, and residue 25 an isoleucine and a leucine. Choosing the best occupied atom for each atom name independently would assemble the serine's `OG` beside the proline's `CG` and `CD`, a side chain no protein has. Choosing per residue keeps one whole amino acid.

The rules, in order:

1. Residues are grouped by chain, residue sequence number **and insertion code**. The insertion code is in the key even though Prodes otherwise ignores it and merges residues that differ only by it, which is a separate defect. Leaving it out here would elect a single conformation across two different residues, which would be a new one.
2. Each conformation is ranked by the **median** occupancy of its atoms. Occupancy is usually constant across a conformation, and the median recovers that constant while ignoring an atom refined away from the rest; occupancy is not always constant, and 1CBN residue 34 carries both 0.80 and 1.00 under one letter. A mean is dragged by such strays, and a sum would prefer whichever conformation has more atoms rather than the one that is there more of the time.
3. Ties go to the letter appearing first in the file, which in conventionally written files is `A`. This is not a corner case: a large share of disordered residues are written at exactly equal occupancy. The comparison holds a small tolerance, because a median of two-decimal occupancies is not always the same float for the same nominal value: 0.04 and 0.37 give exactly 0.205 while 0.01 and 0.40 give 0.20500000000000002, and 125 reachable median values have that property. Comparing exactly would decide those on float noise and the tie-break would never run.
4. Atoms with no letter belong to every conformation and are always kept. This is the ordinary case, an ordered backbone with a disordered side chain, and it is what most disordered residues look like.
5. **Nothing is carried over from a losing conformation.** A winning conformation that is missing an atom stays missing it. An earlier version of this change completed the winner from the runners-up, on the theory that a conformation written only as far as it is ordered would otherwise leave the residue truncated. Measured over thousands of disordered residues, that rule never once supplied a heavy atom: a conformation is written short precisely when it is the minor one, so the best occupied conformation is never the less complete of the two. What it did supply was hydrogens carrying the other rotamer's coordinates, several landing a fraction of an Angstrom from an atom they are not bonded to — on 4NZU it placed an `HG23` 0.86 Å from an `OG1` — and a duplicate of any atom the residue also wrote without a letter. It was removed, and disabling it changed no feature on any structure tested.
6. The residue takes the name of the winning conformation, not of whichever atom arrived first. In 1CBN residue 22 the serine conformers are written before the proline that wins, so naming the residue after its first atom would report a serine holding a proline's atoms.

## What is reported

One warning per structure naming how many residues were disordered and how many atoms were dropped, and a further line for each residue whose winning conformation carries a different amino acid, since that changes the sequence rather than only the geometry.

The same counts go into `prodes_run.json` as `alternate_conformers`. That matters more than it looks: the output bundle ships the **input** file unchanged, so the structure a viewer draws still holds every conformation while the features describe only one.

## Caveats

**A renamed residue changes the reported sequence and mass.** On 1CBN, residues 22 and 25 become proline and leucine where the old parser reported serine and isoleucine, which moves the molecular weight by +10.04 Da. Nothing else reports the substitution, so the warning and the run record are the only notice a user gets.

**A pKa file no longer agrees with a renamed residue.** `redo_pkas` matches predicted pKas to residues by residue name. A pKa file is generated by propka3 from the original coordinates, which still hold every conformation, so for a residue whose name changed the entry is dropped with a warning. Harmless for the serine and proline of 1CBN residue 22, since neither titrates. A real error for a pair such as aspartate and asparagine, where one titrates and the other does not.

**A tie between two amino acids is decided by record order.** Where a residue is modelled as two different amino acids at equal occupancy, nothing but the order of the records separates them, and the reported sequence follows that order. 1ALX chain A residue 11 is a tyrosine and a tryptophan at a median occupancy of 0.50 each. Reordering the file would change the residue Prodes reports and the molecular weight with it. The choice is deterministic for a given file, not physically determined.

**The choice is between letters, not between chemical species.** Each letter is ranked on its own. Where one amino acid is spread over two letters and its rival over one, the rival is compared against each letter separately rather than against their total. In 1EJG residue 22 the serine occupies letters B and C at 0.21 and 0.22 while the proline holds letter A at 0.57, so the proline wins either way; a structure where the split mattered would be decided differently by a rule that aggregated first.

**Insertion codes are still ignored.** Residues 30 and 30A are still read as one residue. This is very visible on Kabat-numbered antibodies: in 4NZU, H100 through H100H arrive as a single residue of 111 atoms holding 35 distinct atom names. Separate defect, separate issue.

**No conformational averaging.** Prodes describes one conformation, the best occupied one. It does not average descriptors over the ensemble, which is what some commercial packages do. A structure whose alternates are close to evenly occupied is genuinely ambiguous, and the run record is what tells you it was.

## Where the code is

`prodes.io.conformers`: `elect_conformers` makes the choice, `rank_letters` and `conformer_rank` order the conformations, and `residue_key` says which records belong to one residue. The alternate location and the occupancy each arrive as their own field on an `AtomRecord`, read by `prodes.io.pdb_reader`. The parse reads every record into memory first, because which conformation to keep cannot be decided until all of them have been seen.

From version 7.1 the records come from `Bio.PDB.PDBParser`, driven with a record collecting `StructureBuilder` of prodes' own rather than through Biopython's `Structure` entity tree. That is deliberate: Biopython represents disorder rather than resolving it, and its defaults for resolving it would undo the rules above. See `docs/pdb_reading.md`.
