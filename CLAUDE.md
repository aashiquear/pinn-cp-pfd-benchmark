# Agent entry point

1. `PROJECT_GUIDELINES.md` binds every conversation (language, citations,
   provenance, accuracy tiers, compute, reproducibility, framework changes).
2. `.claude/FINDINGS.md` holds the author's decisions, measured facts and
   traps that are not derivable from the code. Read it before running anything.
3. `docs/run_ledger.md` lists every run with its deck, scripts, revision,
   host, cores and outcome. Add a row for every run you make.
4. Results go in a numbered `docs/0N_*.md`; workflow and reproduction steps are
   separate documents (`docs/workflow_*.md`, `docs/reproduce_*.md`).
5. Ask the author before any multi-run campaign; check one case first.
