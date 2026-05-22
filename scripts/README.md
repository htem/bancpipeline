# `scripts/` — miscellaneous helpers

One-shot shell helpers that don't belong in a pipeline subdirectory. Not part of any recurring chain.

## Files

- `seeds_to_production.sh` — Convert whole-brain alignment `seeds.csv` → `seeds_production.csv` by flipping every `is_holdout=TRUE → FALSE`. Used to promote a held-out alignment run to the production seeding pass (no holdout — every GT cell becomes an anchor).

## Executable for users?

Yes — `bash scripts/seeds_to_production.sh <seeds.csv>` from the repo root. Operates on a local CSV; no HMS access needed.

See [`alignment/README.md`](../alignment/README.md) for the alignment-pipeline context.
