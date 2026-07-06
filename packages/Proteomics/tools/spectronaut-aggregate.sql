-- Aggregate a long-form Spectronaut Report (peptide/precursor/fragment-level rows)
-- down to one row per (protein-group × condition × replicate), in the shape
-- packages/Proteomics/src/parsers/spectronaut-parser.ts expects.
--
-- This is the bundled, documented manual fallback (D-05) for a Spectronaut report
-- too large for the in-browser streaming importer, AND the D-04 equivalence oracle
-- the Wave-0 golden test pins to.
--
-- Template tokens __IN_PATH__ and __OUT_PATH__ are replaced by the shell wrapper
-- before this script is passed to duckdb. (DuckDB's COPY ... TO does not accept
-- variables, so we can't use SET VARIABLE here.)

-- Streaming aggregation — duckdb spills to disk if it can't fit memory.
PRAGMA memory_limit='8GB';
PRAGMA threads=8;

COPY (
  SELECT
    "PG.ProteinGroups",
    -- DIVERGENCE FROM /tmp/spectronaut-aggregate.sql: the two carry-along
    -- any_value() SELECT terms for the gene-symbol and protein-accession columns
    -- (present in the /tmp original) were DROPPED here. The package's Spectronaut
    -- data (files/demo/spectronaut-hye-mix.tsv and the synthetic precursor
    -- fixture) does not carry those two columns and
    -- pivotSpectronaut/aggToPivotResult never consume them, so they are not part
    -- of the parity contract; keeping them would Binder-Error duckdb on the
    -- fixture (ignore_errors does NOT cover Binder Errors) and the entire Wave-0
    -- oracle chain (golden → sidecar → 12-03 golden test) would have no source.
    --
    -- ============================ WARNING ============================
    -- REFERENCE-FILE-ONLY: the CASE below is a one-off correction for the
    -- mislabeled reference file `2026-05-13 BP DMD WT.tsv` ONLY. Cross-tab
    -- against R.FileName showed every WT* filename tagged DMD and every DMD*
    -- filename tagged WT (24/24); this restores the convention so DE direction
    -- matches what the filenames imply. It is a STRUCTURAL NO-OP on any data
    -- whose R.Condition is not literally 'DMD'/'WT' — including the CondA/CondB
    -- synthetic fixture, which is exactly why this same script doubles as the
    -- D-04 golden oracle without perturbing the golden. The streaming TS
    -- aggregator MUST NOT port this flip (RESEARCH Pitfall 1 — the single
    -- highest-risk parity trap). Do not "generalize" or remove it either: the
    -- manual-fallback path still needs it for the mislabeled reference file.
    -- =================================================================
    CASE "R.Condition"
      WHEN 'DMD' THEN 'WT'
      WHEN 'WT'  THEN 'DMD'
      ELSE "R.Condition"
    END                              AS "R.Condition",
    "R.Replicate",
    any_value("R.FileName")          AS "R.FileName",
    -- PG.Quantity is constant per (protein-group, run) in Spectronaut output;
    -- max() collapses the precursor-level duplicates safely.
    max(TRY_CAST("PG.Quantity" AS DOUBLE)) AS "PG.Quantity",
    -- Best precursor q-value backing this protein in this run.
    min(TRY_CAST("EG.Qvalue" AS DOUBLE))   AS "EG.Qvalue"
  FROM read_csv(
    '__IN_PATH__',
    delim='\t',
    header=true,
    sample_size=-1,           -- scan everything for type inference; columns are messy
    ignore_errors=true,       -- tolerate stray malformed lines in a 2.6 GB file
    nullstr=['', 'NaN', 'NA']
  )
  WHERE
    -- Mirror the parser's precursor-level filter: drop EG rows with q-value > 0.01.
    -- Non-numeric q-values (e.g. 'Profiled') pass — same as the TS parser.
    (TRY_CAST("EG.Qvalue" AS DOUBLE) IS NULL OR TRY_CAST("EG.Qvalue" AS DOUBLE) <= 0.01)
    -- Drop decoys / contaminants up front; parser does this too but it's cheaper here.
    AND "PG.ProteinGroups" IS NOT NULL
    AND "PG.ProteinGroups" NOT LIKE 'CON\_\_%' ESCAPE '\'
    AND "PG.ProteinGroups" NOT LIKE 'REV\_\_%' ESCAPE '\'
  GROUP BY "PG.ProteinGroups", "R.Condition", "R.Replicate"
) TO '__OUT_PATH__' (FORMAT csv, DELIMITER '\t', HEADER);

-- Summary so we can sanity-check the collapse ratio.
SELECT
  count(*)                                          AS output_rows,
  count(DISTINCT "PG.ProteinGroups")                AS proteins,
  count(DISTINCT ("R.Condition" || '_' || "R.Replicate")) AS samples,
  list(DISTINCT "R.Condition")                      AS conditions
FROM read_csv('__OUT_PATH__', delim='\t', header=true);
