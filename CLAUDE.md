# CLAUDE.md

--START OF CODING GUIDELINES--

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

## 5. Environment

Never use `pip`, `pip-tools`, `poetry`, or `conda` directly.

### Running Python Code

- Run a script: `uv run <script>.py`
- Run a tool: `uv run pytest` / `uv run ruff` / `uv run mypy`
- Launch a REPL: `uv run python`


When adding or modifying code in `src/dust_hbm/`, add or update tests in `tests/` to cover the new behaviour. Keep tests fast (no MCMC, CPU-only) — trace numpyro models with `seed`+`trace` for a single forward pass instead.

---
**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

--END OF CODING GUIDELINES--

## 6. Main CLI

```bash
python scripts/make_simlib.py <opsim.db> [options]
# Key options:
#   --Nfields 50000        number of fields to sample
#   --ddf_tags             separate WFD/DDF fields
#   --host_file hosts.fits host galaxy library
#   --wgtmap wgtmap.txt    SNANA weight map for host resampling
#   --author "Name"        header metadata
```

## 7. Architecture

**Data flow:**
```
OpSim SQLite → OpSimSurvey → formatObs() → SNANA_Simlib / SNSIM_obsfile
```

**`opsimsummaryv2/summary_opsim.py` — `OpSimSurvey`**
- Reads OpSim SQLite via SQLAlchemy, filters/processes observation table
- Builds `sklearn.BallTree` spatial index for fast RA/Dec queries
- `compute_hp_rep()`: aggregates observations into HEALPix pixels (nside configurable)
- `sample_survey()`: draws N random HEALPix fields, returns observation sets
- `formatObs()`: converts raw OpSim columns → PSF, ZPT, sky noise (per-band calibration)
- `get_survey_hosts()`: spatial join via GeoPandas polygons + multiprocessing (`host_joiner`)

**`opsimsummaryv2/sim_io.py` — output writers**
- `SNANA_Simlib`: writes SIMLIB format (observation library) + HOSTLIB (host galaxy table)
- `SNSIM_obsfile`: writes Parquet observations + YAML survey config (alternative sim code)

**`opsimsummaryv2/utils.py` — helpers**
- `host_joiner()`: parallelizable spatial join (used by `get_survey_hosts`)
- `host_resampler()`: resamples host catalog by SNANA weight map
- `download_rubinlsst_baseline_dbfile()`: fetches OpSim DB from SLAC

## 8. Key Data Formats

- **Input**: OpSim SQLite (Rubin LSST simulation output), optional FITS host library, optional SNANA WGTMAP
- **Output SIMLIB**: text format required by SNANA; one `LIBID` block per pointing
- **Output HOSTLIB**: tab-separated table appended to SIMLIB or standalone
- **Output Parquet**: columnar observation table + sidecar YAML config for snsim

## 9.  Docstring Style

Google-style (Parameters / Returns / Notes sections).
