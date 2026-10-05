---
description: Implement a day-N step from implementation-steps/ following PLAN.md
argument-hint: <step-file> (e.g. day-1.md or implementation-steps/day-1.md)
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
---

# Implement step: $ARGUMENTS

You are implementing **exactly one** step of the TRISO-CFD project. The step file is: `$ARGUMENTS`

## Scope rule (read first, do not violate)

Implement **only** what the step file at `$ARGUMENTS` explicitly requires. Nothing else.

- Do NOT implement work from other day files, even if they appear "next" or "prerequisite-ish."
- Do NOT add functions, files, tests, figures, or YAML entries that this step does not name.
- Do NOT preemptively scaffold for future days (no stub modules, no placeholder configs, no "we'll need this later").
- PLAN.md and CLAUDE.md are **context and constraints**, not a to-do list. Use them to shape *how* you implement this step, never to expand *what* you implement.
- If the step depends on artifacts from earlier days that do not yet exist, stop and report the gap — do not silently build them.
- If something in the step is ambiguous, pick the narrowest reasonable interpretation and note the assumption in the final report rather than expanding scope to cover both readings.

Every file you touch and every symbol you add must trace back to a specific sentence in `$ARGUMENTS`. If you cannot point to that sentence, do not write it.

## 1. Load context (do these reads in parallel)

- `PLAN.md` — project goal, scope, cases, repo layout (§2.2), deliverables, resources.
- `.claude/CLAUDE.md` — mandatory engineering rules for this repo.
- The step file itself. Resolve `$ARGUMENTS` as follows:
  - If it already contains a `/`, use it as-is.
  - Otherwise prepend `implementation-steps/` (so `day-1.md` → `implementation-steps/day-1.md`).
- Any prior day step files in `implementation-steps/` numbered lower than this one — they establish what should already exist.
- Current repo state: `ls`/`Glob` the directories touched by this step (e.g. `src/`, `params/`, `tests/`, `results/`, `notes/`) to see what is already there. Do NOT reinvent files that exist; extend them.

## 2. Plan before writing code

Write a short plan (5–10 bullets) that lists:
- Every file you will create or modify, mapped to the repo layout in PLAN.md §2.2.
- Every function/class signature you will add, with units on all physical arguments.
- Every unit test you will add and what published or hand-worked value it checks against (per CLAUDE.md rule 4).
- The "Done when" acceptance criteria copied verbatim from the step file, and how each one will be verified.
- Any parameter values needed — these must go in `params/*.yaml` with a source comment (paper + page, or `assumed`), never inlined as magic numbers (CLAUDE.md rule 2).

Do not start editing until the plan is written out in the response.

## 3. Implement

Follow these rules strictly:

- **SI units everywhere.** Docstrings state units on every physical quantity.
- **No invented MFiX keywords.** If this step touches `.mfx` files, cross-check every keyword against the installed MFiX version's keyword reference. If uncertain, say so and stop rather than guess (CLAUDE.md rule 1).
- **Parameters in YAML, not code.** New physical values go in `params/particles.yaml`, `params/gases.yaml`, or `params/geometry.yaml` with a `# source:` comment. Scripts read from YAML.
- **Correlations get unit tests.** Every function added to `src/correlations.py` gets a test in `tests/` that pins it to a hand-worked or published number, with the source cited in a comment.
- **Reuse over rewrite.** If `src/gasprops.py` or `src/correlations.py` already exists, extend it — don't create a parallel module.
- **Notebooks last.** If the step touches `analytics.ipynb`, implement and test the underlying Python first, then have the notebook import from `src/` and call it. Notebook cells should be thin.
- **Figures go to `results/`** with the exact filename the step specifies. Every figure caption states the averaging window if it shows time-averaged data (CLAUDE.md rule 6).
- **Raw simulation output is gitignored.** Only commit scripts, YAML inputs, small CSVs, and figures (rule 7). If `.gitignore` is missing entries this step would produce, add them.

## 4. Verify acceptance criteria

For each bullet in the step's "Done when:" line:

1. Run the check (execute the tests, regenerate the figure, print the hand-calc comparison, etc.) via `Bash`.
2. Report the result plainly. If a check fails, do NOT paper over it — per CLAUDE.md rule 9, first write a short note explaining the likely physical or numerical cause, then propose (and only then apply) a fix.
3. Do not claim the step is done unless every criterion passes.

## 5. Log the run

Append a row to `results/run_log.csv` (create the file with a header row if it doesn't exist) recording: date, step file, key parameters touched, files changed, tests passed/failed, wall-clock time for any simulation runs, and status (`done` / `partial` / `blocked`). This is CLAUDE.md rule 8 — it applies to every step, not just MFiX runs.

## 6. Write the day summary

For every step implemented, produce a markdown summary at `notes/dayN_summary.md` (where `N` matches the step file). Model it on `notes/day1_summary.md` — same structure, same tone, same depth. It is the human-readable record of what this day added to the project, and later days will read it as context.

The file MUST contain these three sections, in this order:

1. **`## What was implemented`** — a one-paragraph framing of the step's deliverable, followed by a bulleted list of the key artifacts (source files, YAML inputs, tests, notebooks, figures) with a one-line purpose for each. Use the same paths that appear in §3 of your report.

2. **`## How it works`** — a numbered walkthrough of the actual mechanics: each function/class or pipeline stage, what equation or algorithm it implements, what inputs and outputs it has (with units), and how the pieces compose. Include the verification step: which hand calc or published number the tests pin to, and the numeric agreement achieved.

3. **`## What it models (physically)`** — explain the physics the step captures and, more importantly, *why it matters* for the coater. Tie each modeled quantity to a downstream operating decision (spouting velocity, fountain height, residence time, sensitivity to gas temperature or particle density, etc.). This is where the step earns its place in PLAN.md — make the connection explicit.

Style rules for the summary:

- Use SI units on every physical quantity, and state them inline (e.g. `µ [Pa·s]`, `ρ_p [kg/m³]`).
- Reference specific published values or hand-calc results with the numeric agreement (e.g. "hand calc gives 6.44 m/s; solver returns 6.442 m/s").
- Cite PLAN.md sections when a modeling choice traces back to a scope decision (`PLAN §1.3`, etc.).
- Do not restate CLAUDE.md rules — assume the reader knows them.
- Do not include a changelog, task list, or run-log data — that lives in `results/run_log.csv`.
- Length: roughly the same as `notes/day1_summary.md` (≈40–60 lines of prose). Longer is fine if the step genuinely warrants it; shorter is a signal the step was under-documented.

The summary is a deliverable of the step, not an afterthought. Do not mark the step `done` in `results/run_log.csv` until `notes/dayN_summary.md` exists and covers all four sections.

## 7. Report back

End with a concise summary in the chat response:
- Files created/modified (with paths), including `notes/dayN_summary.md`.
- Tests added and their pass/fail status.
- Which "Done when" criteria are satisfied and which (if any) are not.
- Anything the user should decide before moving to the next day.

Keep this chat summary tight — two short paragraphs or a short bulleted list, not a wall of text. The long-form write-up belongs in `notes/dayN_summary.md`, not here.
