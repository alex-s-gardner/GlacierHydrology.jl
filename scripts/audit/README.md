# Audit scripts

Read-only scripts that reproduce findings of the GEMB pipeline review
(`notes/gemb_pipeline_review_2026-10-13.md`) from the production data. They write nothing.
Each script's header gives the command to run and the output it printed when written.
They need the JPL data directories listed in `src/local_paths.jl`.

| Script | Findings | What it reports |
|---|---|---|
| `discharge_units_and_surface_masks.jl` | 6, 7, 8 | Antarctic discharge stored as Gt/yr and converted with δice (110.6 → 101.4 Gt/yr). The geotile discharge column before the fixes, with only the double conversion, and corrected. The RGI7/RGI6 glacier-area ratio per geotile. |

Run from the repository root:

    ~/.juliaup/bin/julia +release --project=. scripts/audit/discharge_units_and_surface_masks.jl
