# SpaceMarkersExperiment Pipeline Methods

**Date:** 2026-04-20
**Target branch:** `claude/sme-pipeline-methods` (off `claude/keen-chatelet`)
**Merges into:** `claude/keen-chatelet`

## Motivation

The `SpaceMarkersExperiment` (SME) class introduced on `claude/keen-chatelet`
serves as a unified container for SpaceMarkers inputs and outputs, but the
pipeline functions used in the step-by-step vignette
(`find_all_hotspots`, `get_pairwise_interacting_genes`, `get_im_scores`, the
directed-pipeline functions, and the LR-scoring chain) operate on raw
matrices/data.frames. To build up an SME incrementally, the current
`SpaceMarkersStepByStep_vignette.Rmd` reaches into internal slots with
`sme@spacemarkers$results$... <- ...` and
`S4Vectors::metadata(sme)$hotspots$... <- ...`. The goal of this change is
to close that gap: every pipeline step gets an SME-in/SME-out method, and
every result slot gets a setter generic symmetric to its existing reader.

End-user result (per-step chain becomes):

```r
sme <- sme |>
    find_all_hotspots() |>
    get_pairwise_interacting_genes(mode = "DE", analysis = "enrichment") |>
    get_im_scores() |>
    calculate_overlap_undirected()
```

## Scope

In scope (added in this PR):

- Setter generics and `SpaceMarkersExperiment` methods for every result slot
  already exposed via a read accessor.
- Conversion of the pipeline step functions to S4 generics with two methods
  each: the existing body as the default (matrix/data.frame inputs), plus an
  `"SpaceMarkersExperiment"` method.
- Refactor of `.undirected_SpaceMarkers_sme` and `.directed_SpaceMarkers_sme`
  in `R/OneSpaceMarkers.R` to chain the new SME methods (dogfooding — no
  behavior change for end users of `SpaceMarkers()`).
- Rewrite of `vignettes/SpaceMarkersStepByStep_vignette.Rmd` to use the new
  chainable API, and re-knit to HTML.
- testthat coverage in `tests/testthat/test-SpaceMarkersExperiment-pipeline.R`.

Out of scope:

- Plot functions. `plot_spatial_data_over_image`, `plot_im_scores`, etc. keep
  their current signatures.
- `load10X`. Already SME-returning.
- The top-level `SpaceMarkers()` wrapper signature. Only its two internal
  helpers change.
- The SpatialExperiment vignette, which exercises the top-level wrapper and
  does not demonstrate per-step composition.

## Design

### Setter generics

Added to `R/AllGenerics.R`. Each mirrors a current reader.

| Setter generic | Writes to | Special |
|---|---|---|
| `hotspots<-` | `metadata(x)$hotspots[[type]]` | `type = c("undirected","pattern","influence")`, same `match.arg` as the reader |
| `interactions<-` | `metadata(x)$interactions` | assigns the full named list; per-pair update is not needed for the pipeline |
| `influence_map<-` | `metadata(x)$influence` | |
| `undirected_scores<-` | `x@spacemarkers$results$undirected_scores` | |
| `directed_scores<-` | `x@spacemarkers$results$directed_scores` | |
| `overlap_scores<-` | `x@spacemarkers$results$overlap_scores` | |
| `lr_scores<-` | `x@spacemarkers$results$lr_scores` | |
| `analysis_type<-` | `x@spacemarkers$analysis` | Value must be `%in% c("undirected","directed","both")`. Validation runs via the existing `setValidity("SpaceMarkersExperiment", ...)` since `@spacemarkers<-` on an S4 object re-validates. |

All setters are defined in `R/SpaceMarkersExperiment-methods.R` next to their
readers. Each setter preserves sibling slots untouched.

Ligand and receptor score slots under `metadata(x)$ligand_scores` /
`metadata(x)$receptor_scores` are written by the
`calculate_gene_set_score` and `calculate_gene_set_specificity` SME methods
respectively. Rather than expose setters for these intermediate slots, the
writes are performed by the SME methods directly via `S4Vectors::metadata(x)`.

### Pipeline function conversions

Each of the following becomes an S4 generic. The existing function body
becomes the default method (dispatched on the current primary argument).
A second method dispatching on `"SpaceMarkersExperiment"` is added.

**Generic argument naming.** Each generic keeps the first-argument name used
by the current function (`data` for `get_pairwise_interacting_genes` and
`calculate_gene_scores_directed`, `spPatterns` for `find_all_hotspots` and
`calculate_influence`, `df` for `find_hotspots_gmm`, `hotspots` for
`calculate_overlap_undirected`, `pat_hotspots` for
`calculate_overlap_directed`, `SpaceMarkers` for `get_im_scores`,
`ligand_scores` for `calculate_lr_scores`). Dispatch is on that first arg.
This preserves back-compat for **all existing call sites** — both positional
and named — because no arg name changes. SME method bodies may alias the
first arg locally (e.g., `sme <- data`) for readability.

| Function | Default-method dispatch | SME method reads | SME method writes |
|---|---|---|---|
| `find_all_hotspots` | `x = "data.frame"` | `.sme_spPatterns(x)`, `spatial_params(x)` | `hotspots(x, "undirected") <- hs` |
| `find_hotspots_gmm` | `df = "data.frame"` | `.sme_spPatterns(x)` or `influence_map(x)` (chosen by `type` arg on SME method) | `hotspots(x, type) <- hs` |
| `calculate_influence` | `spPatterns = "data.frame"` | `.sme_spPatterns(x)`, `spatial_params(x)` | `influence_map(x) <- inf` |
| `get_pairwise_interacting_genes` | `data = "ANY"` (matrix-like) | `.sme_expr(x)`, `.sme_spPatterns(x)`, `spatial_params(x)`, `hotspots(x, "undirected")` | `interactions(x) <- res`; updates `x@spacemarkers$params$mode`, `$analysis_method`, `$min_overlap` |
| `get_im_scores` | `SpaceMarkers = "list"` | `interactions(x)` | `undirected_scores(x) <- get_im_scores(ints)` |
| `calculate_gene_scores_directed` | `data = "ANY"` | `.sme_expr(x)`, `hotspots(x, "pattern")`, `hotspots(x, "influence")` | `directed_scores(x) <- scores` |
| `calculate_overlap_undirected` | `hotspots = "data.frame"` | `hotspots(x, "undirected")` | `overlap_scores(x) <- ov`; `analysis_type(x) <- "undirected"` (or `"both"` if directed already present) |
| `calculate_overlap_directed` | `pat_hotspots = "data.frame"` | `hotspots(x, "pattern")`, `hotspots(x, "influence")` | `overlap_scores(x) <- ov`; `analysis_type(x) <- "directed"` (or `"both"`) |
| `calculate_gene_set_score` | existing signature (`data`, `gene_set`) | `directed_scores(x)` and `params(x)$lr_pairs$ligand.symbol` | writes `metadata(x)$ligand_scores` |
| `calculate_gene_set_specificity` | existing signature | `directed_scores(x)` and `params(x)$lr_pairs$receptor.symbol` | writes `metadata(x)$receptor_scores` |
| `calculate_lr_scores` | `ligand_scores = "ANY"` | `metadata(x)$ligand_scores`, `metadata(x)$receptor_scores`, `params(x)$lr_pairs` | `lr_scores(x) <- lr` |

### Shared helpers

Two unexported helpers in `R/SpaceMarkersExperiment-methods.R` (used by every
SME method):

```r
.sme_expr <- function(sme) {
    if ("logcounts" %in% SummarizedExperiment::assayNames(sme))
        SummarizedExperiment::assay(sme, "logcounts")
    else
        SummarizedExperiment::assay(sme, 1L)
}

.sme_spPatterns <- function(sme) {
    data.frame(
        barcode = colnames(sme),
        as.data.frame(SpatialExperiment::spatialCoords(sme)),
        as.data.frame(spatial_patterns(sme)),
        check.names = FALSE,
        row.names = colnames(sme)
    )
}
```

The SME pipeline-methods call these, then delegate to the default method with
the extracted inputs, then write the result via the matching setter. No
business logic is duplicated — every SME method is purely extract-call-pack.

### Error handling

Each SME method checks the one or two prerequisite slots it depends on. When
missing, `stop()` with a message naming the step the caller should run:

- `get_pairwise_interacting_genes(sme)` when `hotspots(x, "undirected")` is
  `NULL`: `"Run find_all_hotspots(x) before get_pairwise_interacting_genes()."`
- `get_im_scores(sme)` when `interactions(x)` is `NULL`:
  `"Run get_pairwise_interacting_genes(x) before get_im_scores()."`
- `calculate_gene_scores_directed(sme)` when either `hotspots(x, "pattern")`
  or `hotspots(x, "influence")` is `NULL`: name the specific missing step
  (`find_hotspots_gmm(x, type = "pattern")` or `calculate_influence(x)` + the
  matching GMM call).
- `calculate_overlap_undirected(sme)` when `hotspots(x, "undirected")` is
  `NULL`: name `find_all_hotspots()`.
- `calculate_overlap_directed(sme)` when either pattern/influence hotspot
  slot is `NULL`: name the missing step.
- `calculate_lr_scores(sme)` when ligand_scores, receptor_scores, or
  `params(x)$lr_pairs` is missing: name the missing prerequisite.

No defensive re-validation beyond that. Internal consumers are trusted.

### Dogfooding

`.undirected_SpaceMarkers_sme` becomes:

```r
.undirected_SpaceMarkers_sme <- function(sme, cpus = 1, genes = NULL,
                                         min.gene.expr = 10, ...,
                                         returnSME = TRUE) {
    sme <- .apply_sme_filters(sme, genes = genes, min.gene.expr = min.gene.expr)
    sme <- find_all_hotspots(sme)
    sme <- get_pairwise_interacting_genes(sme, mode = "DE",
                                          analysis = "enrichment",
                                          minOverlap = 10, workers = cpus, ...)
    sme <- get_im_scores(sme)
    sme <- calculate_overlap_undirected(sme)
    if (!returnSME) return(undirected_scores(sme))
    sme
}
```

The directed helper is analogous. The filter step (gene/barcode filtering and
`get_spatial_parameters` fallback) stays in an internal `.apply_sme_filters()`
helper — factored out of the current `.process_sme_input()` — because it isn't
naturally a pipeline method (it mutates multiple slots at once).

### Data flow

Undirected:
```
sme --find_all_hotspots--> metadata$hotspots$undirected
    --get_pairwise_interacting_genes--> metadata$interactions, params
    --get_im_scores--> results$undirected_scores
    --calculate_overlap_undirected--> results$overlap_scores, analysis
```

Directed:
```
sme --calculate_influence--> metadata$influence
    --find_hotspots_gmm(type="pattern")--> metadata$hotspots$pattern
    --find_hotspots_gmm(type="influence")--> metadata$hotspots$influence
    --calculate_gene_scores_directed--> results$directed_scores
    --calculate_overlap_directed--> results$overlap_scores, analysis
```

LR:
```
sme --calculate_gene_set_score--> metadata$ligand_scores
    --calculate_gene_set_specificity--> metadata$receptor_scores
    --calculate_lr_scores--> results$lr_scores
```

## Testing

All tests written before the corresponding implementation (TDD). Added to
`tests/testthat/test-SpaceMarkersExperiment-pipeline.R`. Fixtures reuse
`system.file("extdata", "CoGAPS_result.rds", package = "SpaceMarkers")` and
whatever is already in `tests/testthat/assets/` for the existing pipeline
test — no new test data files.

1. **Setter/getter symmetry.** For each new setter, set a value then read it
   back with the matching accessor; assert equality.
2. **SME method parity.** For each pipeline step, run the matrix/data.frame
   default method and the SME method on the same inputs; assert the slot
   written by the SME method equals the default method's return.
3. **Prerequisite errors.** For each SME method with a prerequisite, call it
   on a fresh SME (no prior step run) and assert `expect_error(...,
   regexp = "<named prerequisite>")`.
4. **Undirected chain.** Build SME, chain the four undirected SME methods,
   assert final SME's `undirected_scores`, `overlap_scores`,
   `analysis_type`, and `hotspots(x, "undirected")` are non-null and
   well-formed.
5. **Directed chain.** Same, for the five directed SME methods, including LR
   when `params(x)$lr_pairs` is provided.
6. **Back-compat.** Since the default method *is* the existing function body,
   the existing test files (`test-findPatternHotspots.R`,
   `test-getInteractingGenes.R`, `test-getPairwiseInteractingGenes.R`,
   `test-getOverlapScores.R`, etc.) must continue to pass unchanged after
   the refactor. No new assertions needed beyond ensuring these pass.
7. **Dogfooding regression.** Call `SpaceMarkers(sme, directed = FALSE)` and
   `SpaceMarkers(sme, directed = TRUE, lr_pairs = ...)` and assert the
   results match what the pre-refactor `.undirected_SpaceMarkers_sme` /
   `.directed_SpaceMarkers_sme` produce. The baseline is captured in the
   same test file by running the pre-refactor code paths once at the top
   (or, equivalently, by asserting the resulting SME structure matches the
   shape documented in the class docstring).

## Vignette

`vignettes/SpaceMarkersStepByStep_vignette.Rmd` is rewritten:

- Undirected section: replaces the block that currently does
  `sme_ud@spacemarkers$results <- list(...)` and
  `S4Vectors::metadata(sme_ud)$hotspots$undirected <- hs` with a four-step
  chain of SME methods.
- Directed section: analogous rewrite.
- LR section: if present in the current vignette, also rewritten.
- Narrative around the setters and the explicit step sequence is kept — the
  point of the vignette is still to show what each step does — but the
  "assemble an SME by hand" boilerplate is deleted.
- Re-knitted to `SpaceMarkersStepByStep_vignette.html` (currently on tip is
  the knit from earlier in this session; will be replaced).

## Naming, documentation, exports

- All new generics and methods go in roxygen-annotated blocks; `NAMESPACE`
  regenerated via `devtools::document()`.
- All setters get `@rdname` matching their reader so they appear on the same
  help page.
- New generics and their SME methods inherit `@rdname` from the existing
  function's documentation, with a `@param x A SpaceMarkersExperiment or the
  corresponding matrix/data.frame` note explaining the dispatch.

## Migration / back-compat

- Every existing call site of the promoted functions keeps working because
  the default method *is* the current function body. This is the whole point
  of the two-method pattern.
- `OneSpaceMarkers.R` end-to-end behavior is unchanged (same outputs for same
  inputs).
- No deprecations introduced.

## Open items / deferred

- Could add a `setMethod("show", ...)` enhancement that prints which pipeline
  step is next runnable given which slots are populated. Not in this PR.
- Could add `inherits()`-style "fluent" aliases (`run_hotspots(sme)` etc.) —
  intentionally not added; the functions keep their current names so that
  users of the matrix API and readers of the existing manual have continuity.
