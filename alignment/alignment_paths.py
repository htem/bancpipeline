"""alignment_paths — D.6.2 self-describing filename builder (Python).

Mirror of `alignment_path()` / `legacy_alignment_path()` in
`alignment/alignment-data-sources.R`. See `alignment/RENAME_NOTES.md`
§D.6.2 for the full grammar.

Reads:
  - env var BANC_ALIGNMENT_NAMING: "v2" (default) or "legacy". Lets a
    producer dual-write or scripts honour the same toggle the R side
    reads.

CLI:
  None — this is a library.

Used by:
  Future migration of alignment/banc-alignment-run.py,
  alignment/banc-alignment-sweep.py, etc. to opt into the v2 names.

Notes:
  - Only the stages used by the paper-run set are wired into
    legacy_alignment_path(); unknown stages raise to flag the gap.
"""
from __future__ import annotations

import datetime as _dt
import os as _os
from typing import Optional as _Optional, Sequence as _Sequence


_FIELD_SEP = "__"
_TAG_SEP = "-"


def _date_str(date) -> str:
    if date is None:
        return _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%d")
    if isinstance(date, (_dt.datetime, _dt.date)):
        return date.strftime("%Y%m%d")
    s = str(date).replace("-", "")
    if not (len(s) == 8 and s.isdigit()):
        raise ValueError(
            f"alignment_path(): date must be YYYY-MM-DD or YYYYMMDD; got {date!r}"
        )
    return s


def _tag(prefix: str, value) -> _Optional[str]:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    return f"{prefix}{_TAG_SEP}{s}"


def alignment_path(
    stage: str,
    query: str,
    target: str,
    *,
    region: _Optional[str] = None,
    side: _Optional[str] = None,
    vq: _Optional[str] = None,
    vt: _Optional[str] = None,
    syn: _Optional[str] = None,
    extra: _Optional[_Sequence[str]] = None,
    date=None,
    ext: str = "csv",
    dir: str = ".",
    scheme: str = "auto",
) -> str:
    """Build a self-describing alignment-output filename.

    Schema:
        {stage}__{query}_{target}[__region-…][__side-…][__vq-…][__vt-…]
              [__syn-…][__<extra>]__{YYYYMMDD}.{ext}

    Args:
        stage:  pipeline stage tag (e.g. "prep-banc-meta",
                "align-cosine-tier", "validate-holdout-accuracy"). Must
                not contain "__".
        query:  query dataset name, lowercase-with-dashes (e.g. "banc").
        target: target dataset name (e.g. "fafb", "manc").
        region: optional region tag (e.g. "whole-brain", "optic-lobe").
                Emitted as "region-<value>".
        side:   optional "right" | "left" | "both". "side-<value>".
        vq:     optional query version pin. "vq-<value>".
        vt:     optional target version pin. "vt-<value>".
        syn:    optional synapse-table version ("v2" | "v3"). "syn-<value>".
        extra:  optional list of preset-defined extra tags (each
                lowercase-with-dashes, no "__").
        date:   ISO date, datetime, or YYYYMMDD string. Default today UTC.
        ext:    file extension (no leading dot). Default "csv".
        dir:    output directory. Default "." (current).
        scheme: "auto" (read $BANC_ALIGNMENT_NAMING, default "v2"),
                "v2", or "legacy".

    Returns:
        File path string. If `dir` is "." or None, just the filename.

    Raises:
        ValueError: invalid stage / date / extra tag format.
    """
    if scheme == "auto":
        scheme = _os.environ.get("BANC_ALIGNMENT_NAMING", "v2").lower()
    if scheme == "legacy":
        return legacy_alignment_path(
            stage=stage, query=query, target=target,
            region=region, side=side, vq=vq, vt=vt, syn=syn,
            extra=extra, date=date, ext=ext, dir=dir,
        )
    if scheme != "v2":
        raise ValueError(f"alignment_path(): scheme must be 'v2' or 'legacy'; got {scheme!r}")
    if _FIELD_SEP in stage:
        raise ValueError(f"alignment_path(): stage may not contain '__'; got {stage!r}")

    parts = [stage, f"{query}_{target}"]
    for prefix, value in (
        ("region", region), ("side", side),
        ("vq", vq), ("vt", vt), ("syn", syn),
    ):
        t = _tag(prefix, value)
        if t is not None:
            parts.append(t)
    if extra:
        for e in extra:
            if e is None:
                continue
            es = str(e).strip()
            if not es:
                continue
            if _FIELD_SEP in es:
                raise ValueError(f"alignment_path(): extra tag may not contain '__'; got {es!r}")
            parts.append(es)
    parts.append(_date_str(date))
    fname = _FIELD_SEP.join(parts) + "." + ext
    if dir in (".", None):
        return fname
    return _os.path.join(dir, fname)


def legacy_alignment_path(
    stage: str,
    query: str,
    target: str,
    *,
    region: _Optional[str] = None,
    side: _Optional[str] = None,
    vq: _Optional[str] = None,
    vt: _Optional[str] = None,
    syn: _Optional[str] = None,
    extra: _Optional[_Sequence[str]] = None,
    date=None,
    ext: str = "csv",
    dir: str = ".",
) -> str:
    """Map the v2 parameters to the legacy filename convention.

    Only stages used by the paper-run set are wired up. Producers may
    opt into dual-write by setting BANC_ALIGNMENT_NAMING=legacy and
    writing the same content to both paths during the D.6.2 migration.
    """
    _LEGACY_INFIX = {"optic-lobe": "optic", "whole-brain": "brain"}
    infix = ""
    if region is not None:
        if region not in _LEGACY_INFIX:
            raise ValueError(
                f"legacy_alignment_path(): no legacy infix for region={region!r}"
            )
        infix = _LEGACY_INFIX[region]

    side_l = side or ""

    def _q_t_infix_side() -> str:
        return "_".join(p for p in (query, infix, side_l) if p)

    def _q_target_infix_side() -> str:
        return "_".join(p for p in (query, target, infix, side_l) if p)

    if stage == "prep-banc-meta":
        fname = f"{_q_t_infix_side()}_meta.csv"
    elif stage == "prep-target-meta":
        bits = "_".join(p for p in (infix, side_l) if p)
        fname = f"{target}_{bits}_meta.csv"
    elif stage == "prep-banc-edges":
        fname = f"{_q_t_infix_side()}_edgelist.feather"
    elif stage == "prep-target-edges":
        bits = "_".join(p for p in (infix, side_l) if p)
        fname = f"{target}_{bits}_edgelist.feather"
    elif stage == "prep-seeds":
        fname = f"{_q_t_infix_side()}_seeds.csv"
    elif stage == "prep-seeds-production":
        fname = f"{_q_t_infix_side()}_seeds_production.csv"
    elif stage == "prep-nblast":
        fname = f"{_q_target_infix_side()}_nblast.csv"
    elif stage == "prep-capacity":
        fname = f"{target}_type_capacity.csv"
    elif stage == "prep-forbidden":
        fname = "forbidden-matches.csv"
    elif stage in (
        "validate-holdout-accuracy",
        "validate-holdout-confusions",
        "validate-nt-mismatches",
        "validate-mismatches",
    ):
        suffix = stage[len("validate-"):]
        fname = f"{_q_t_infix_side()}_{suffix.replace('-', '_')}.csv"
    elif stage == "validate-discrepancies":
        extra_tag = ""
        if extra:
            tags = [str(e) for e in extra if e and str(e).strip()]
            if tags:
                extra_tag = "_" + "_".join(tags)
        fname = f"{_q_t_infix_side()}_discrepancies{extra_tag}.csv"
    elif stage == "align":
        fname = f"{_q_t_infix_side()}_alignment.{ext}"
    elif stage.startswith("align-"):
        rest = stage[len("align-"):].replace("-", "_")
        fname = f"{_q_t_infix_side()}_alignment_{rest}.{ext}"
    elif stage.startswith("ntac"):
        rest = stage.replace("-", "_")
        fname = f"{_q_t_infix_side()}_{rest}.{ext}"
    else:
        raise ValueError(f"legacy_alignment_path(): no legacy mapping for stage={stage!r}")

    if dir in (".", None):
        return fname
    return _os.path.join(dir, fname)


__all__ = ["alignment_path", "legacy_alignment_path"]
