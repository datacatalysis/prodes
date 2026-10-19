"""Looks up what each Prodes feature is called and what it means.

This module calculates nothing. It is a lookup table over the features that
``prodes.run.calculate`` produces, so that pipelines and plotting code can turn a
column heading into a human-readable label without hard-coding strings.

A **feature code** is the short identifier Prodes writes as the CSV column
heading, such as ``SurfEpMeanFormal``. Every lookup here takes a feature code and
returns something more readable. The direction is one-way by design: there is no
route back from a label or description to a feature code.

    from prodes.feature_dictionary import FeatureDictionary

    fd = FeatureDictionary()

    label = fd.get_plot_name("SurfEpMeanFormal")        # 'Surface EP mean (formal)'
    title = f"Impact of {fd.get_description(code)} on R2"
    fd.get_long_description(code)                       # how it is calculated
    fd.get_unit(code)                                   # 'Da', 'v', ... or None
    fd.get_original_explanation(code)                   # wording from the 2024 paper
    fd.get_reason_dropped("SurfEpMeanAverage")          # None when the feature is kept

Prodes writes a reduced, non-redundant feature set by default and the full legacy
set under ``--full-features``. Both code lists are available, so a pipeline can
line up datasets calculated under either setting without re-running Prodes or
diffing CSV headers by hand:

    fd.get_feature_codes()             # the 65 written by default
    fd.get_feature_codes(full=True)    # the 116: the 105 legacy columns and the 11 patch ones
    fd.get_dropped_feature_codes()     # the 51 the default leaves out
    fd.get_feature_codes(residue_sets=["CDR"])  # plus the columns of a residue set

    # cut a full-feature CSV down so it lines up with a reduced-feature one
    df = fd.drop_redundant_features(pd.read_csv("full_output.csv"))

The reduced list is in the order Prodes writes those columns. The full list is
the reduced codes followed by the full-only ones, which is not the column order of
a ``--full-features`` run, since that interleaves the two; read the CSV header if
you need the exact order.

``ID`` is a row label rather than a feature, so it is in neither list; use
:data:`ID_COLUMN` when you need it.

Everything, including which feature codes exist, comes from two human-readable
YAML files in ``prodes/data``: ``features_reduced.yaml`` for the 65 calculated by
default and ``features_full_only.yaml`` for the other 51, the latter carrying the
reason each one was dropped. Those files are the single source of truth, so the
list of features is never written down twice and cannot disagree with itself. Each
is read exactly once, at import. Edit them by hand: adding a feature means adding
an entry there, and the tests will tell you if the entry and the code disagree.

What is still worth checking is that they agree with the code, which
``tests/test_feature_dictionary.py`` does by comparing them against the real
output of ``calculate()``.

Why the other 51 are dropped by default is set out in
``docs/redundant_feature_analysis.md``.

The last 11 of the default features are the surface patch features, marked
``per_residue_set`` in ``features_reduced.yaml``. A run given a residue set NAME
writes each of them a second time as ``<code>_NAME``, preceded by
``SurfArea_NAME``. ``get_patch_feature_codes`` returns the 11, every getter also
accepts a patch code with a residue set suffix, such as ``SurfMhpPatchArea_CDR``,
and ``get_feature_codes(residue_sets=...)`` includes the suffixed columns. The
method is in ``docs/surface_patches.md``.
"""

import importlib.resources

import yaml

#: The row-label column Prodes writes before the features. Not a feature itself.
ID_COLUMN = "ID"

REDUCED_DICTIONARY_FILE = "features_reduced.yaml"
FULL_ONLY_DICTIONARY_FILE = "features_full_only.yaml"


def _read_features(filename: str, key: str = "features") -> dict:
    """Returns one mapping, ``features`` by default, of one YAML dictionary file."""

    with importlib.resources.open_text("prodes.data", filename, encoding="utf-8") as handle:
        return yaml.safe_load(handle)[key]


# Each file is read exactly once, at import. The two are kept as separate dicts
# because that is already the split callers ask for: the reduced set, the
# full-only set, and the union of the two. The files are expected to partition the
# feature set; tests/test_feature_dictionary.py checks that they do.
_REDUCED_ENTRIES = _read_features(REDUCED_DICTIONARY_FILE)
_FULL_ONLY_ENTRIES = _read_features(FULL_ONLY_DICTIONARY_FILE)
_ALL_ENTRIES = {**_REDUCED_ENTRIES, **_FULL_ONLY_ENTRIES}

#: The non-redundant subset written by default, in the order Prodes writes it.
REDUCED_FEATURE_CODES = tuple(_REDUCED_ENTRIES)

#: Written only under ``--full-features``, in the order Prodes writes them.
DROPPED_FEATURE_CODES = tuple(_FULL_ONLY_ENTRIES)

#: Every feature code the original Prodes algorithm produces: the reduced set
#: followed by the full-only set. Note that this is not the column order of a
#: ``--full-features`` run, which interleaves the two; read the CSV header if you
#: need that.
FULL_FEATURE_CODES = REDUCED_FEATURE_CODES + DROPPED_FEATURE_CODES

# The patch features are part of the reduced set; they are the entries marked
# per_residue_set, which a run repeats for every residue set it is given.
_PATCH_ENTRIES = {code: entry for code, entry in _REDUCED_ENTRIES.items() if entry.get("per_residue_set")}
_RESIDUE_SET_ONLY_ENTRIES = _read_features(REDUCED_DICTIONARY_FILE, "residue_set_only_features")

#: The surface patch features, in the order Prodes writes them. With a residue
#: set NAME, each is written again as ``<code>_NAME``.
PATCH_FEATURE_CODES = tuple(_PATCH_ENTRIES)


def residue_set_feature_codes(name: str) -> list[str]:
    """Returns the columns Prodes writes for the residue set NAME, in output order."""

    return [f"{code}_{name}" for code in (*_RESIDUE_SET_ONLY_ENTRIES, *PATCH_FEATURE_CODES)]


def _patch_entry(feature_code: str) -> dict | None:
    """Returns the entry of a patch code, either whole-surface or ``<code>_NAME``, or None.

    Only patch codes take a residue set suffix, so ``Area_CDR`` is not resolved
    to ``Area``.
    """

    code, separator, name = feature_code.rpartition("_")
    # The rule PatchSettings applies to a set name: ASCII letters and digits.
    if separator and name.isascii() and name.isalnum():
        return _PATCH_ENTRIES.get(code) or _RESIDUE_SET_ONLY_ENTRIES.get(code)

    return None


class FeatureDictionary:
    """Translates a Prodes feature code into a plot label or an explanation.

    Every getter takes a feature code, the short identifier Prodes writes as the
    CSV column heading (``SurfEpMeanFormal``), and returns something more
    readable. Construction is cheap; the underlying files are read once per
    process and shared between instances.

        fd = FeatureDictionary()
        label = fd.get_plot_name("SurfEpMeanFormal")

    Raises:
        KeyError: from any getter given a string that is not a feature code.
    """

    def __init__(self) -> None:
        self._entries = _ALL_ENTRIES

    def __repr__(self) -> str:
        return f"{type(self).__name__}({len(REDUCED_FEATURE_CODES)} reduced, {len(FULL_FEATURE_CODES)} full feature codes)"

    def __contains__(self, feature_code: str) -> bool:
        return feature_code in self._entries or _patch_entry(feature_code) is not None

    # --- the feature code lists -------------------------------------------------

    def get_feature_codes(self, full: bool = False, residue_sets=()) -> list[str]:
        """Returns the feature codes, excluding ``ID``.

        The reduced set comes back in the order Prodes writes those columns. The
        full set is the reduced codes followed by the full-only ones, which is not
        the column order of a ``--full-features`` run: that interleaves the two.
        Read the CSV header if you need the exact column order. The columns of
        the residue sets come last, in the order Prodes writes them.

        Args:
            full: return the full set instead of the reduced default.
            residue_sets: the names of the run's residue sets, in the order
                they were given.

        Returns:
            A new list of feature codes.
        """

        codes = list(FULL_FEATURE_CODES if full else REDUCED_FEATURE_CODES)
        for name in residue_sets:
            codes.extend(residue_set_feature_codes(name))

        return codes

    def get_patch_feature_codes(self) -> list[str]:
        """Returns the 11 surface patch feature codes, which a run repeats for every residue set."""

        return list(PATCH_FEATURE_CODES)

    def get_dropped_feature_codes(self) -> list[str]:
        """Returns the feature codes in the full set that the reduced set leaves out."""

        return list(DROPPED_FEATURE_CODES)

    # --- code to human-readable text --------------------------------------------

    def get_plot_name(self, feature_code: str) -> str:
        """Returns the short axis or legend label for a feature code.

        >>> FeatureDictionary().get_plot_name("SurfEpMeanFormal")
        'Surface EP mean (formal)'
        """

        return self.get_entry(feature_code)["plot_name"]

    def get_description(self, feature_code: str) -> str:
        """Returns the one-line description of a feature code."""

        return self.get_entry(feature_code)["short_explanation"]

    def get_long_description(self, feature_code: str) -> str:
        """Returns the full explanation of a feature code, including how it is calculated."""

        return self.get_entry(feature_code)["long_explanation"]

    def get_unit(self, feature_code: str) -> str | None:
        """Returns the unit of measurement for a feature code.

        Returns:
            The unit as given by Neijenhuis et al. 2024, or None when the feature
            is dimensionless or the paper does not document it.
        """

        return self.get_entry(feature_code).get("unit")

    def get_original_explanation(self, feature_code: str) -> str | None:
        """Returns the description of a feature code from the original paper.

        Kept alongside :meth:`get_description` rather than replacing it, because
        the two write-ups emphasise different things: the paper is terser, while
        the longer text says more about how the value is actually computed.

        Returns:
            The wording from Supplemental Table 1 of Neijenhuis et al. 2024, or
            None for the few features the paper does not document.
        """

        return self.get_entry(feature_code).get("original_explanation")

    def get_reason_dropped(self, feature_code: str) -> str | None:
        """Returns why the default feature set leaves a feature out.

        Args:
            feature_code: a Prodes feature code.

        Returns:
            The explanation, or None when the feature is in the reduced set and so
            was never dropped.
        """

        return self.get_entry(feature_code).get("drop_note")

    # --- whole entries ----------------------------------------------------------

    def get_entry(self, feature_code: str) -> dict:
        """Returns everything known about one feature code.

        Returns:
            A new dict, safe to mutate. Full-only features carry three extra keys
            (``dropped_because``, ``redundant_with``, ``drop_note``).
        """

        if feature_code in self._entries:
            return dict(self._entries[feature_code])

        entry = _patch_entry(feature_code)
        if entry is not None:
            return dict(entry)

        raise KeyError(
            f"{feature_code!r} is not a Prodes feature code. Use get_feature_codes(full=True) for the {len(FULL_FEATURE_CODES)} "
            "valid codes; a patch code may also carry a residue set suffix, such as SurfMhpPatchArea_CDR."
        )

    def get_dictionary(self) -> dict:
        """Returns every entry of the reduced and full sets, keyed by feature code.

        The patch features are included under their plain codes; the residue
        set columns, such as ``SurfMhpPatchArea_CDR``, are not, but get_entry
        looks them up.

        Returns:
            A new dict of new dicts, safe to mutate.
        """

        return {code: dict(entry) for code, entry in self._entries.items()}

    # --- acting on a dataframe --------------------------------------------------

    def drop_redundant_features(self, dataframe, keep_id: bool = True):
        """Cuts a full-feature Prodes output down to the reduced feature set.

        Not a getter: this performs an action on the frame you pass in. Use it to
        make a dataset calculated with ``--full-features`` directly comparable
        with one calculated under the default.

        Args:
            dataframe: a DataFrame whose columns are Prodes output.
            keep_id: keep the ``ID`` column in front of the features when present.

        Returns:
            A new DataFrame with only the reduced features, in output order.

        Raises:
            KeyError: if any reduced feature is missing, listing what was absent.
        """

        missing = [code for code in REDUCED_FEATURE_CODES if code not in dataframe.columns]
        if missing:
            raise KeyError(f"input is missing {len(missing)} of the {len(REDUCED_FEATURE_CODES)} reduced features: {missing}")

        columns = list(REDUCED_FEATURE_CODES)
        if keep_id and ID_COLUMN in dataframe.columns:
            columns.insert(0, ID_COLUMN)

        return dataframe[columns]
