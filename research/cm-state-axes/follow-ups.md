# Follow-ups

## FU-1: Understand gaps that block any direction decision or Deepen (from C001 review F4, 2026-11-07)

Still unread:

- Li 2025 Supplementary Data 2 and 5, and the second halves of the Li 2025 and Honkoop 2019 full
  texts
- Wu 2016 full text, and GSE74652
- Ben-Yair 2019 and Sallin 2015 full texts
- Beisaw 2020
- A forward-citation search from Honkoop 2019 and Li 2025
- A sourced zebrafish cell-cycle set (for example Tirosh 2016 orthologs)

Restart condition: these must be done before a direction decision is saved.

## FU-2: original label-bearing files (2026-11-07)

`data/hu2022/*.h5ad` still carry the authors' labels. A human should move them out of agent reach
(for example into `holdout/`), because agents may not touch `holdout/`. Owner: the user.

## FU-3: Zotero exports (2026-11-07)

The four CSVs in `data/zotero/` are identical. Owner: the user, to supply the separate exports if
they exist.

## FU-5: teammates' unmerged branches hold existing work for this study (2026-11-07 15:45)

`main` was force-pushed (fea5089). The new main only strips annotation details from `docs/data.md`
and `scripts/get_data.sh`. Eleven unmerged branches add:

- `origin/ad/markers-li2025-spatial`: Li 2025 Supplementary Data 2 parsed into
  `markers/sources/li2025_bzm.tsv`, with Ensembl IDs and a filter for neighbouring non-CM cell
  types. This closes part of FU-1, and its Ensembl IDs can settle F1 (the *vmhc* symbol).
- `origin/hk/docs-research-plan`: the unit of replication is the **heart**, not the library
  (Hr2a/Hr2b, H8a/H8v, Hr6a/Hr6v and Hr7a/Hr7v are one heart each). It proposes a binomial
  depth-thinning check, decoys matched on detection rate, and a dissociation control axis.
- `origin/ad/guardrails`, `ad/loop-config` and `ad/write-guards`: `config/axes.yaml` and
  `config/checks.yaml` (decoys, a dissociation axis, partial-correlation checks), plus changes to
  `AGENTS.md` and `.claude/settings.json`.
- `origin/ad/markers-literature`: its own `markers/*.yaml`, which will **collide** with this
  study's untracked `markers/*.yaml` on merge.
- `origin/ad/data-agent-input`: a second allowlist script, which duplicates `src/zf/io.py::sanitize_h5ad`.

Owner: the user and the team decide what to merge. This study reads the branches with
`git show` and does not check them out.

## FU-6: the trap test was hinted for *ttn.2* (2026-11-07 15:45)

*ttn.2* is named in the source paper (Hu et al. 2022), so any *ttn.2* result is not found unaided (R-E6). The agreement check is not void, because no label values were seen.

## FU-7: doublet threshold failure (C005 review F11, 2026-11-07 16:00)

Scrublet's automatic threshold called about 0 doublets in every library. Rerun with a per-library
expected doublet rate (from 10x loading) or with a second method (scDblFinder or SOLO), then test
whether the agent's state calls are enriched among doublets. Restart condition: before any claim
that the state calls are free of doublets.

## FU-8: supervisor dates ran ahead of the clock (2026-11-07)

Several supervisor timestamps were estimated rather than read from the system clock, and were
corrected in `JOURNAL.md`. Package lesson: read the clock (`date`) before stamping any record.

## FU-9: earlier results were not re-checked after a confound was found (user question, 2026-11-07 16:37)

The chamber confound (F8, C003) led to ventricle-only analysis from C004 on. But the C002
results (null, depth ties, rank against score_genes, per-gene depth flags), the C003 ambient
panel and the C005 doublet correlations were computed on all cardiomyocytes and never re-run
within ventricle cells. The user's question exposed this.

General failure: when a new confound is found, every earlier result that rests on the
confounded population must be re-checked, not only the next step. Action: C008 re-runs them
within ventricle cells. Package lesson for the caller: add "re-check earlier results under each
new confound" to the review step.

## FU-4: package lesson, a subagent tool was unavailable (2026-11-07)

`AGENTS.md` asks for one subagent per atomic step, but the researcher's runtime had no subagent
tool. This is not specific to this study. Report it to the caller.

## FU-1 status update (researcher, C006, 2026-11-07 16:10)

- Done:
  - forward-citation search from Honkoop 2019 and Li 2025 (OpenAlex, R-E19);
  - Li 2025 Supplementary Data 2 and 5 opened (structure, cluster lists, targeted rows).
- Still open:
  - Wu 2016 full text and GSE74652;
  - the second half of Honkoop 2019;
  - Ben-Yair 2019 and Sallin 2015 full texts;
  - Beisaw 2020;
  - a sourced zebrafish cell-cycle set;
  - a full read of Neufeld 2024.
- Partly met. A direction decision for A (Li SD2-based) rests on read tables. B and C still rest
  on summaries.
