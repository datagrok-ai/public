# JMP Plugin for Datagrok — Specification

> **Purpose:** This spec defines the Datagrok JMP plugin: a browser-side reader for SAS JMP data tables (`.jmp`). It is written for the developer (human or AI agent) building the package. It covers scope, architecture, the file format as far as it is known, type mapping, testing, and the phased roadmap. The [Minitab plugin](../../Minitab/docs/MINITAB_PLUGIN.md) is the structural template.

**Status:** Phase 1 built, 2026-10-21; Phase 2 in progress, 2026-10-23 (v0.1.0, not yet released). The parser matches JMPReader.jl on all 19 reference files; 116 package tests and 12 BDD scenarios pass on a local stand. §11 lists what is left.

---

## Table of Contents

1. [Problem and Goals](#1-problem-and-goals)
2. [Scope](#2-scope)
3. [Package Layout](#3-package-layout)
4. [The `.jmp` File Format](#4-the-jmp-file-format)
5. [Type Mapping](#5-type-mapping)
6. [User-Facing Behaviour](#6-user-facing-behaviour)
7. [Error Handling and Robustness](#7-error-handling-and-robustness)
8. [Performance](#8-performance)
9. [Testing and Validation](#9-testing-and-validation)
10. [Licensing and Attribution](#10-licensing-and-attribution)
11. [Roadmap](#11-roadmap)
12. [Open Questions](#12-open-questions)

---

## 1. Problem and Goals

### Problem

JMP (SAS Institute) is widely used in pharma R&D, CMC, and manufacturing for exploratory statistics, DOE, and quality work. In many organisations, analyses live as `.jmp` files on scientists' local drives and team shares. They are hard to find, hard to share, and invisible to anyone who doesn't have a JMP licence. The reasoning behind decisions made with those analyses is lost with them.

Datagrok can't open `.jmp` files today. Users have to open each one in JMP and export it to CSV or Excel. That step loses types and metadata, and it doesn't scale to backfilling hundreds of historical files.

### Goals

1. **Open any `.jmp` data table in Datagrok** by drag-and-drop, File > Open, or a double-click in the file browser, with no JMP installation and no server-side component.
2. **Preserve the data faithfully**: column names, numeric/character/date/time/duration types, and missing values, matching what JMP shows.
3. **Preserve enough metadata to be useful**: JMP version, save time, and (later) column properties, value labels, and row states.
4. **Enable bulk ingestion** of historical JMP files into shared Datagrok spaces (Project Portal), with searchable metadata.
5. **Later, re-create saved JMP analyses as live Datagrok viewers** (from table scripts), so the analysis comes across, not just the data.

### Non-goals

- Writing `.jmp` files. The write path is undocumented, and JMP imports CSV/Excel/Parquet well enough for the round trip.
- Running JSL (JMP Scripting Language). We may *parse* a narrow subset of saved table scripts to build viewers (Phase 3), but we never execute JSL.
- JMP projects (`.jmpprj`), journals (`.jrn`), standalone scripts (`.jsl`), and add-ins (`.jmpaddin`). These are out of scope for v1 and listed in [Open Questions](#12-open-questions).

---

## 2. Scope

| Phase | Deliverable | Section |
|---|---|---|
| **1 — Data bridge (MVP)** | Import + preview of `.jmp` data tables; core types; missing values; compression | §4–§9 |
| **2 — Fidelity** | Row states, value labels, column properties (units, notes, spec limits), formula text, table variables | §11 |
| **3 — Analyses** | Parse saved table scripts and re-create common platforms as Datagrok viewers | §11 |
| **4 — Portal ingestion** | Bulk import of a folder or share; metadata extraction for findability | §11 |

---

## 3. Package Layout

A new standalone package, `public/packages/Jmp` (`@datagrok/jmp-reader`), mirrors `public/packages/Minitab`:

```
packages/Jmp/
├── package.json            # @datagrok/jmp-reader, category "Data Access"; runtime deps: datagrok-api, u2
├── README.md               # user-facing: what opens how, limitations, how to test
├── CREDITS.md              # JMPReader.jl (MIT) attribution + test fixture provenance
├── docs/JMP_PLUGIN.md      # this file
├── files/
│   ├── jmpreader/          # JMPReader.jl's 19 test files (MIT), unchanged
│   ├── golden/             # what JMPReader.jl reads from each (tools/jmp2json.jl), committed
│   ├── tests/              # made here: unsupported-column, truncated, not-a-jmp, saved-bar-chart, saved-fit-model; jmp15-graph-builder (saved by JMP 15, anonymized)
│   └── public/             # eight openly licensed research files saved by JMP 10–18 (§9.2, CREDITS.md §5)
├── tools/
│   ├── jmp2json.jl         # golden-output generator (Julia + JMPReader.jl), dev-only
│   ├── jmp-corpus.ts       # `npm run corpus -- <dir>`: read a folder of .jmp files, report (§9.4)
│   ├── jmp-large.ts        # `npm run large -- --mb 100`: write a large synthetic file and time reading it (§8)
│   ├── jmp-survey.ts       # bundled to one jmp-survey.mjs: survey saved analyses where the files are, shareable summary (§9.4)
│   └── make-script-fixtures.ts # writes the synthetic saved-analysis fixtures (§9.2)
├── bdd/                    # Gherkin features → Playwright specs (libraries/bdd), §9.5
│   ├── features/           # open, preview, problem-files, analyses
│   ├── fixtures/           # the files the features drop or pick from disk
│   └── generated/          # compiled specs, committed, never edited by hand
└── src/
    ├── package.ts          # fileHandler + fileViewer
    ├── package-test.ts     # browser entry + testNode() for the headless Node pass
    ├── formats/jmp/
    │   ├── jmp-types.ts        # constants (port of constants.jl), JmpTable/JmpColumn, JmpFormatError
    │   ├── jmp-reader.ts       # ByteReader: bounds-checked little-endian cursor
    │   ├── jmp-parser.ts       # header, metadata, per-column decode (port of metadata.jl + column.jl)
    │   ├── jmp-decompress.ts   # gzip sections via DecompressionStream
    │   ├── jmp-scripts.ts      # saved table scripts: where they sit in the file (§4.5)
    │   ├── jmp-properties.ts   # column properties and value labels in a column's header (§4.6)
    │   ├── jsl.ts              # a reader for declarative JSL: tokens → tree, never evaluated
    │   ├── jsl-to-viewers.ts   # a saved analysis → Datagrok viewer specs (§11, Phase 3); no platform imports
    │   ├── jmp-to-dataframe.ts # JmpTable → DG.DataFrame (types, tags, semTypes, row states); adds analyses to a view
    │   └── jmp-viewer.ts       # u2 views: the Browse preview and the can't-open state
    └── tests/
        ├── golden.ts               # compares a parsed table with JMPReader.jl's output (no platform imports)
        ├── large-jmp.ts            # makes a large .jmp from compressed.jmp for the §8 measurements
        ├── jmp-reader-tests.ts     # header, reference, tables (runtests.jl ported), robustness; Node pass
        ├── jmp-script-tests.ts     # JSL reader, saved scripts, analyses → viewers; Node pass
        └── jmp-dataframe-tests.ts  # DataFrame conversion, row states, analyses on a view, the preview; browser pass
```

**Why a standalone package, and not Crux?** Crux is a chemistry engine (molecules, SMARTS, fingerprints). A `.jmp` file is a generic data table, so a reader in Crux would share no code with the rest of it. A combined "stats tool formats" package (Minitab + JMP) is a reasonable future refactor. It doesn't change the parser, so v1 ships as its own package.

**Registration.** A `fileHandler` for `jmp` (bytes → `DataFrame[]`) and a `fileViewer` for `jmp` (a `FileInfo` → a view). The shell hands a file from disk (drop, File > Open) to a package's file viewer before its handler (`core/client/xamgle/lib/src/commands/file/file_commands.dart`, `_openFile`), so the viewer tells the two apart by the file's connection: none means the file came from disk, and it opens straight as a table; a file on a server gets the preview (§6).

Because parsing runs in the browser, no Docker container or server function is needed. Drag-and-drop, File > Open, and the file browser all work through the handler.

---

## 4. The `.jmp` File Format

The format is **binary, proprietary, and undocumented**. SAS has not published it, and users have long asked for it in a long-running "Open up JMP table format" thread on the JMP Community wish list. The only serious independent implementation is **[JMPReader.jl](https://github.com/jaakkor2/JMPReader.jl)** (Julia, MIT, Jaakko Ruohio, v0.2.4, actively maintained as of Oct 2026). It is about 540 lines of source. Its author reports it reads 605 of the JMP 17 sample files and 612 of the JMP 18 sample files, plus about 1,000 real-world files.

**Implementation strategy: a faithful TypeScript port of JMPReader.jl**, kept structurally close to the original (same function boundaries, same constant tables), so upstream fixes can be diffed and back-ported. Record the upstream commit we ported from in `CREDITS.md`.

What follows is the format as JMPReader.jl understands it. Everything is little-endian.

### 4.1 Header

| Offset | Bytes | Meaning |
|---|---|---|
| 0 | 4 | Magic `FF FF 00 00` |
| 4 | 12 | Version triple, three `UInt32` values. **Supported:** `(7, 4, 0x00010001)`. **Rejected:** `(3, 4, 3)`, an older format, and files starting `00 00 FF FF`: the same old format big-endian, as classic Mac JMP wrote it (seen in 2007–08 research files, with 31-character space-padded column names in `windows-1252`). Anything else → "Unknown JMP file format". |

### 4.2 Table metadata

The metadata block is found **heuristically**: search for the byte sequence `07 00 08 00 00 00`, then step back 38 bytes. From there:

| Field | Type | Notes |
|---|---|---|
| `nrows` | `Int64` | |
| `ncols` | `Int32` | |
| *(unknown)* | 5 × `Int16` | |
| `charset` | string, 4-byte length prefix | NUL-padded |
| *(unknown)* | 3 × `UInt16` | |
| `savetime` | `Float64` | Seconds since **1904-01-01** (JMP epoch) |
| *(unknown)* | `UInt16` | |
| `buildstring` | string, 4-byte length prefix | e.g. `"... Version 18.0.1"`; JMP version parsed from it |

Next comes the **column index**, also found heuristically. Scan for runs of `FF FF`, skip 10 bytes, and read `n_visible: UInt32` and `n_hidden: UInt32`. Accept the position when `n_visible + n_hidden == ncols`. Then read:

- `idx_visible: UInt32[n_visible]`, `idx_hidden: UInt32[n_hidden]`
- `colwidths: UInt16[ncols]` (display widths)
- 7 × `UInt32` (unknown)
- zero or more optional blocks, each introduced by `FC FF` / `FD FF` / `FE FF` / `FF FF` + `Int64` length (skipped)
- `ncols` again as `Int32` (sanity check: it must match)
- `coloffsets: Int64[ncols]`, the absolute file offsets of each column's data block

At each column offset, a 2-byte-length-prefixed string holds the column name.

> ⚠️ Both seeks are heuristics. JMPReader.jl's own source notes that the optional-block detection can misfire if a column count happens to match those byte patterns. The port must turn each assertion into a clear error (§7), not silent garbage.

### 4.3 Column blocks

Each column block starts with the column name, followed by six descriptor bytes `dt1..dt6` that set storage and display type:

| `dt1` | Meaning |
|---|---|
| `0x01` | Numeric, uncompressed |
| `0x0A` | Numeric, gzip-compressed |
| `0x02` | Character, uncompressed |
| `0x09` | Character or row state, gzip-compressed (row state when `dt2 == 0x03`) |
| `0x03` (with `dt2 == 0x03`) | 8-byte integer storage |
| `0xFF`, `0xFE`, `0xFC` | Byte integers (1/2/4-byte, width in `dt5`) |

**Compression.** When `dt1 ∈ {0x09, 0x0A}`, the block contains a section marker `EF BE FE CA` (`0xCAFEBEEF`), then `gziplen: UInt64` and `gunziplen: UInt64`, then a standard **gzip** stream. The browser's `DecompressionStream('gzip')` handles it with no dependency. Fall back to `pako` only if a target browser lacks it.

**Numeric data.** The values are the *last* `nrows × width` bytes of the (decompressed) block. `dt6` gives the width for `dt1 ∈ {0x01, 0x0A}`: 1 → `Int8`, 2 → `Int16`, 4 → `Int32`, 8 → `Float64`. `dt4` (display format) and `dt5` (input format) then classify the column:

| Classification | Rule (from JMPReader.jl `constants.jl`) |
|---|---|
| Plain number / fixed decimal | `dt4 == dt5 ∈ {0x00, 0x03, 0x42, 0x43, 0x44, 0x59, 0x60, 0x61, 0x63}` or `dt5 ∈ {0x5E, 0x63}`. `0x61` is ours, not JMPReader.jl's: the p-value format of `Prob>\|t\|` columns (values in [0, 1], JMP 15). |
| Currency | `dt4 == dt5 == 0x5F` |
| Longitude | `dt4 == dt5 ∈ {0x54, 0x55, 0x56}` |
| Latitude | `dt4 == dt5 ∈ {0x51, 0x52, 0x53}` |
| Date | `dt4 ∈ {0x65, 0x66, 0x67, 0x6E–0x72, 0x75, 0x76, 0x7A, 0x7F, 0x88, 0x8B, 0x90}`. `0x90` is ours: a date column saved by JMP 18, every value a midnight. |
| Date-time | `dt4 ∈ {0x69, 0x6A, 0x73, 0x74, 0x77–0x79, 0x7B–0x7E, 0x80–0x82, 0x86, 0x87, 0x89, 0x8A}` |
| Duration | `dt4 ∈ {0x6C, 0x6D, 0x83, 0x84, 0x85}`. The ambiguous `:day:hr:m` / `:day:hr:m:s` formats are treated as durations. |

Dates, date-times and durations are stored as `Float64` seconds from the 1904-01-01 epoch.

**Missing values.** For `Float64`, the sentinel is `NaN`. For integer widths, the sentinel is `typemin + 1` (e.g. `-127` for `Int8`, `-2147483647` for `Int32`).

**Character data.** There are three layouts:

1. **Fixed width**: `(dt3, dt4) == (0, 0) && dt5 > 0`, or `1 ≤ dt3 ≤ 7 && dt4 == 0`. Width is `dt5`, and each value is truncated at the first NUL.
2. **Variable width**: `(dt3, dt4, dt5) == (0, 0, 0)`. A width-type byte (1/2/4) is followed by per-row lengths and then the concatenated UTF-8 bytes. The header offsets differ between the compressed and uncompressed variants (see `column.jl`).
3. **Pooled ("compact" subtype)**: compressed, and the first `Int64` equals the block length. The block holds per-row indices into a string pool (index width of 1/2/4 bytes), then the pool itself (lengths of 1 or 2 bytes). Index 0 means the empty string. This maps directly onto Datagrok's categorical string storage.

**Row states** (`dt1 == 0x09 && dt2 == 0x03`): fixed-width records per row (8 bytes in every file seen). Bytes 6–7 hold the marker index (0–31 map to JMP's marker glyphs; higher values are a Unicode code point). Bytes 1–4 hold the colour, either an RGB triple when byte 4 is `0xFF`, or an index into JMP's 65-colour palette. Bytes 5 and 8 are zero in every fixture. JMPReader.jl doesn't decode the **selected / excluded / hidden / labelled** flags, and no fixture sets them. In JSL a row state is a 16-bit number whose low four bits are selected, excluded, hidden and labelled, then the marker and the colour; the record may hold them in bytes 5 and 8, which the reader keeps as `rowStateFlags` and the corpus run (§9.4) reports when set.

A row-state **column** is data the user saved, one state per row. The table's own row states (which rows are excluded or hidden now) are something else, and where the file keeps them isn't known: no fixture has any.

**Fixed decimal**: JMPReader.jl's source notes that `dt5` of `0x5E` or `0x63` means fixed decimal, with `dt3` the width and `dt4` the number of decimals. No fixture has one (every `0x63` column has `dt4 == dt5`, JMP's "Best"). The reader takes `dt4` as the decimals when `dt5` is one of those, `dt4 != dt5` and `dt4 ≤ 30`, and sets the column's format to that many decimals.

**Currency code**: a currency column's descriptor doesn't vary with the currency. The code is in the column header, as `13 00`, a `u32` length of 5, the three-letter ISO code, then two bytes that are `02 00` in every fixture (perhaps the decimals). The reader looks for that within 256 bytes after the descriptor. Confirmed on all nine columns of `currencies.jmp`.

**Unknown combinations** of `dt1..dt6`: JMPReader.jl logs an error and returns a NaN column. We do the same (§7).

### 4.4 What is *not* decoded yet (by anyone, publicly)

- What a column property *means* beyond notes, units, formulas, value order, spec limits and value labels (§4.6): the rest are kept as JSL text. `Value Colors` needs JMP's colour numbers, which aren't published: small indices (`"Yes" = 4`) and negative numbers that look like packed RGB (`-4354269`) both occur.
- Table variables and table notes
- Saved table scripts: the scan (§4.5) found all 45 scripts in the public corpus (§9.4); long scripts and run-on-open scripts are still unconfirmed
- Row-state flags (selected / excluded / hidden / labelled), and the table's own row states (§4.3)
- Column groups
- The visible/hidden column split: two lists of column positions, visible then hidden. In all 305 public files read (§9.4) the two lists together hold every column exactly once, and the hidden ones are what people hide: survey bookkeeping ("Response Type", "Start Date", "Recorded Date", "Survey duration (seconds)") and source details ("Source Table", "Plan File Name"). Display width doesn't track it. Taken as confirmed, though no file says which columns JMP shows hidden.

These are reverse-engineering targets for Phases 2–3. Each new decoding gets a fixture file and a note in §4.

### 4.5 Saved table scripts (provisional)

Not in JMPReader.jl. Inferred from the `Source` script JMP saved in the `compact_*` fixtures, the only scripts in any fixture:

```
... build string | 05 00 04 00 00 00 00 00 01 00 04 00 | u64 683 | 03 00 | u64 593 | 01 00 | 03 | u32 586 | "Source(Open(...))" | ff ff ...
```

The text is plain UTF-8 JSL, written as `Name(body)`: the script's name called on the script. The nested lengths look like a property block holding a script list, but one example can't confirm that. The reader (`jmp-scripts.ts`) doesn't rely on them. It scans the bytes between the build string and the first column for `0x03`, a `u32` length, and text of that length that reads as a JSL call, and takes `New Script("name", body)` as well as `Name(body)`. A script that has that shape but doesn't parse is kept with its error, so it can be listed. Any surprise in this scan is caught; it can't cost the table's data.

Confirmed on the public corpus (§9.4): up to 9 scripts in one table, all found. A name JSL can't write bare is saved as `Name("1a_logistic-regr")(body)`; JMP 15 writes some such names bare, `HC Met255 %Oxidation vs. Week_5C(…)`. Tables also keep hand-written JSL programs as scripts (loops, dialogs, `New Window`); these are listed as programs, not analyses. Still to confirm: long scripts (compressed?) and scripts JMP marks as run-on-open.

### 4.6 Column properties and value labels (provisional)

Not in JMPReader.jl. The container is the one [jmp-convert](https://github.com/davidmathers/jmp-convert) (MIT) reads; what the entries hold was read from the public corpus (§9.4). After a column's six descriptor bytes come 5 bytes for a numeric column (`dt1` `0x01` or `0x0A`) or 4 for any other, then a property list:

```
u16 count | u64 total length | count × (u16 id | u64 length | payload)
```

| id | Payload | Seen |
|---|---|---|
| `0x07` | JSL properties: each a `u32` length and the property as JSL text, such as `Notes("Was the offer accepted?")`; the list ends with `0xFFFFFFFF` | 91 properties in 13 kinds, in the public corpus and the fixtures: `Formula`, `Notes`, `Units`, `Value Order`, `"Value Ordering"`, `Value Colors`, `Coding`, `Design Role`, `Factor Changes`, `Response Limits`, `Response Probability`, `Intermediate`, `Predicting` |
| `0x0E` | Value labels: `u32 n`, `u32` length, `n` keys (`f64` for a numeric column, NUL-ended text otherwise), `u32 n`, `u32` length, `n` NUL-ended labels | 2 columns (`Validation`: 0, 1, 2 → Training, Validation, Test) |
| `0x04` | A coded text column's value pool (jmp-convert) | not seen; no column fails without it |

Columns with no properties have a count of 0. A header that doesn't add up has no properties (`jmp-properties.ts`); a surprise there can't cost the column's data. Every property in the corpus parses as JSL, once the reader takes `.` (missing), a leading `<<`, `Name("x")(...)` and a string called as a name (`"Value Ordering"(...)`, JMP's form for a property it has no name for).

The vendored `example1.jmp` carries `Formula(:ints * 2)` and `Units("kPa")`. The public files have no licence, so they aren't vendored; the other cases are tested on headers built to the layout above.

---

## 5. Type Mapping

| JMP type | Datagrok column | Details |
|---|---|---|
| Numeric `Float64` | `double` | `NaN` sentinel → null |
| Numeric `Int8/16/32` (byte integers) | `int` | `typemin+1` sentinel → null. Datagrok `int` is 32-bit with null at `-2^31`, so there is no collision. |
| 8-byte integer storage | `bigint` if any value exceeds the 32-bit range, otherwise `int` | |
| Fixed decimal | `double` | `format` = `0.00` with the column's decimals (§4.3; from JMPReader.jl's note, no fixture yet) |
| Currency | `double` | `format` = `#,##0.00`; `units` and the tag `jmp.currency` = the ISO code (`EUR`) |
| Longitude / Latitude | `double` | semType `DG.SEMTYPE.LONGITUDE` / `LATITUDE`, which the GIS package reads |
| Date | `datetime` | Convert the 1904 epoch to ms since 1970: `(s − 2082844800) × 1000`. Set a date-only `format` tag. |
| Date-time | `datetime` | Same conversion |
| Duration | `double` (seconds) | Set the `units` tag to `s`. Datagrok has no native duration type. |
| Character (any layout) | `string` | Decoded with `TextDecoder('utf-8')`. |
| Row state | two `string` columns, `<name> marker` and `<name> color` | The marker column holds JMP's glyph, marker-coded with the nearest Datagrok marker (18 of 32 have one); the colour column holds `#RRGGBB`, colour-coded with itself. Both carry the tag `jmp.rowState` = the JMP column's name. The flags aren't decoded yet (§4.3). |
| A column with value labels | `string` | Each value's label, or the value itself where it has none, as JMP shows it. The values and labels are kept in the `jmp.valueLabels` tag, `[[value, label], ...]`. |
| A column the reader can't decode | Not a column | Left out, listed in the preview and named in a warning balloon (§6, §7). It would hold nothing. |

**DataFrame-level tags** (prefix `jmp.`):

| Tag | Source |
|---|---|
| `jmp.version` | parsed from `buildstring` |
| `jmp.buildString` | `buildstring` |
| `jmp.saveTime` | `savetime`, ISO-8601 |
| `jmp.charset` | `charset` |
| `jmp.scripts` | the saved table scripts, as JSON `[{name, text}]`, present only when there are any (§4.5) |

**Column properties** (§4.6): `Notes` → `description`; `Units` → `units`; `Formula` → `jmp.formula`, the JSL expression (the column keeps JMP's cached values); `Value Order` with a custom order, or `"Value Ordering"` → the column's category order (through its labels); `Spec Limits` → `jmp.specLimits`, `{"lsl", "usl", "target"}`. Every property, mapped or not, is kept as JSL text in `jmp.properties`, `[{"name", "text"}]`.

**Column-level tags**: `jmp.format` (raw `dt1..dt6` as hex, for debugging and future decoding), `jmp.displayWidth`, and `jmp.hidden` = `true` on a column JMP hides (§4.4). A hidden column stays in the table, as in JMP, and is hidden in the grid of a table view the package opens.

Column order follows the file. Duplicate names can't happen in JMP, but names are passed through `DG.DataFrame`'s normal de-duplication anyway.

---

## 6. User-Facing Behaviour

All UI is built in u2 (`libraries/u2`): `appView` hosts the content, the Open button sits in the ribbon, and the grid is u2's `viewers.grid`.

### 6.1 A file in the Browse tree: preview

Clicking or double-clicking a `.jmp` file in the Browse tree shows the preview:

- a summary line: a `JMP <version>` badge, `<rows> rows × <columns> columns` (the columns that come across), and when the file was saved
- a folded **Columns not imported** section, present only when there are any, with a hint (`1 of 12`) and one line per column: its name and why
- an open **Saved analyses** section, present only when the table has saved analyses, with a hint (`1 of 2 open as viewers`) and one line per analysis: its name, the viewers it opens as or why it doesn't, what doesn't come across (a `By()` split, pages, a smoother), and a **Script** link that shows the script's text
- the table in a grid
- **Open** in the ribbon, which opens the table in the workspace with the converted analyses as viewers

Double-click previews rather than imports, as for other formats with a package viewer (Curves' `.pzfx`).

### 6.2 A file from disk: straight to a table

A `.jmp` file dropped onto the window, or picked through File > Open > File... or the Browse toolbar's Open local file icon, opens as one table view named after the file, as CSV or Excel do, with each converted analysis as a viewer titled after its script. Columns and analyses that don't come across are named in warning balloons (`JMP: 1 column not imported: formula`, `JMP: 1 saved analysis not converted: Fit A`).

### 6.3 A file that can't be read

A file that isn't a readable JMP table (not a JMP file, an older format, truncated, empty) opens a view that says so: "This file can't be opened as a JMP table", the reason, and a hint to re-save it in JMP. It opens no table. The same view appears in the Browse preview.

### 6.4 Programmatic

`JmpReader:importJmp(bytes)` is the file handler, callable from scripts and from the Portal ingestion job (Phase 4) without the UI. It throws a `JmpFormatError` for an unreadable file. It returns the table only: it has no view to add analyses to, but the scripts are kept in the `jmp.scripts` tag.

---

## 7. Error Handling and Robustness

The parser relies on heuristics, so failures must be **loud, local, and recoverable**:

| Failure | Behaviour |
|---|---|
| Empty file | Throw "The file is empty" |
| Bad magic | Throw "Not a JMP data table" |
| Old format `(3,4,3)`, little- or big-endian | Throw "This JMP file uses an older format that isn't supported yet. Re-save it in a current JMP version." |
| Unknown version triple | Throw "Unrecognised JMP file version (x.y.z)". Include the triple so we can collect samples. |
| Metadata seek fails / `ncols` mismatch | Throw with the offsets tried. Never return a partially-wrong table. |
| Unknown column descriptor | **Column-local**: the column is left out, listed in the preview, named in a warning. The rest of the table still loads. |
| Decompression error | Column-local, same as above |
| Truncated file (offsets past EOF) | Throw "File appears truncated" |

Every read goes through `ByteReader`, which checks the real buffer length. A row count beyond what the file could hold (more than 1,100× its size, past gzip's ceiling) and a negative block length are rejected up front, so a corrupt header can't make the reader allocate or loop without end. The `JMP: robustness` fuzz test mutates fixtures (byte flips, truncations) and asserts the reader either throws a `JmpFormatError` or returns a table; a 5,100-mutation run during development hit neither a hang nor any other error.

---

## 8. Performance

- **Typed arrays, end to end.** The reader keeps each column the way a Datagrok column stores it (`JmpData` in `jmp-types.ts`): doubles, dates and durations in a `Float64Array` (NaN missing), integers in an `Int32Array` (missing = −2³¹, which is `DG.INT_NULL`), text as distinct values plus an `Int32Array` of per-row indexes. The conversion hands these over with `fromFloat64Array`, `fromInt32Array`, `Column.fromIndexes` and, for dates, `setRawData` on a date-time column (microseconds since 1970). Nothing is boxed per value. An 8-byte numeric block that is 8-byte aligned is copied as one `Float64Array`. Tests and tools read single values through `jmpValue`.
- **Text.** De-duplicating text in a JS `Map` while reading and handing it over with `fromIndexes` is far faster than one string per row: in the browser, for 1.9 M unique strings, 0.4–0.6 s against 7–22 s for `Column.fromList`, and 0.05 s against 0.2 s for 6,000 distinct values. `fromIndexes` also doesn't guess types, so JMP text that reads as numbers stays text (the earlier `fromStrings` would have made it numeric).
- **Lazy, parallel decompression.** Gzip each column independently with `DecompressionStream`, using `Promise.all` across columns. This is the counterpart of JMPReader.jl's threaded column loop.
- **Large files.** `fileHandler` receives the whole file as bytes, so practical limits are browser memory. Target: a 100 MB `.jmp` file (several million rows) opens in under 10 s on a standard laptop. Above roughly 1 GB, show an error suggesting a server-side path (a Phase 4 concern).
- **Web worker.** Run in a worker only if the main-thread parse of the 100 MB target blocks the UI for more than 2 s. Measure before adding it.
- **Yielding.** For tables of 100,000 rows or more, the reader yields to the event loop between columns, so the UI waits at most one column at a time. Text that is all ASCII is decoded once per column and cut, not decoded string by string.
- **Measured on the reference files (Node 22, M-series laptop):** 5–8 M cells/s (10,000 × 7 compressed: 9 ms; 65,537 × 1 pooled strings: 11 ms).
- **Measured at 100 MB (2026-10-22).** The file is synthetic: `compressed.jmp` refilled with 1,886,792 rows of varied values in its seven columns (number and date-time about 95% distinct, short text with 6,000 values, one-letter text, duration, date, and text of 20–130 characters, all unique), 5% missing (`src/tests/large-jmp.ts`, `npm run large`). It's valid for this reader; JMP never opened it.

| Where | Read | DataFrame | Table view | Total | Longest main-thread block |
|---|---|---|---|---|---|
| Node 22 (`npm run large`), 3 runs | 1.7–2.5 s; the table holds 397 MB | | | | |
| Chrome, local stand (`grok test --benchmark`), 4 runs | 3.6–4.4 s | 0.08–0.12 s | 0.4–0.8 s | **4.1–5.1 s** | 0.8–1.2 s while reading; 1.9–2.6 s overall, in the table view's setup |

The target is met. Reading is mostly gzip, which runs off the main thread. Most of the 397 MB is the two million unique strings, which the Datagrok column holds anyway. The stand ran in Docker on the same machine (load average around 10), so the browser numbers vary from run to run.

Before the typed arrays (values in JS arrays, one boxed number or string per row), the same file in Node read in 2.6–3.3 s and held 459 MB, and the text column alone took 7–22 s to become a Datagrok column. Earlier measurements in this section came from a generator whose random numbers repeated after a short cycle (528 distinct long strings), so they understated the cost and are replaced by these.

---

## 9. Testing and Validation

### 9.1 Strategy

The model is the way Crux validates against RDKit: **JMPReader.jl is the reference implementation**. We don't hand-assert expected values for every file. Instead we diff our output against the reference's output.

### 9.2 Fixtures

| Source | Licence | Commit to repo? |
|---|---|---|
| JMPReader.jl `test/*.jmp` (19 files covering byte integers, compact pools of UInt8/16/32, compression, currencies, dates, durations, geographic, long column names, row states, `minusfour`, `bugMWE4`) | MIT, part of the upstream repo | **Yes.** Vendor into `files/` with attribution in `CREDITS.md`. |
| Files made here from the vendored ones (`files/tests/`: a column with a layout byte no JMP writes, a truncated file, CSV under a `.jmp` name) | Ours | **Yes** |
| Synthetic saved analyses (`files/tests/saved-bar-chart.jmp`, `saved-fit-model.jmp`): `compact_UInt8.jmp` with its `Source` script replaced by one of the same length, so nothing in the file moves (`tools/make-script-fixtures.ts`). They test what the package does with an analysis, not where JMP keeps one (§4.5). | Ours | **Yes**, until files saved by JMP with analyses replace them |
| A table saved by JMP 15 with a Graph Builder (`files/tests/jmp15-graph-builder.jmp`, 12 rows). Shared by a user and anonymized in place: the `Source` path replaced by one of the same length, the measured values moved by up to ±0.15. JMP 15 writes the script and column names bare (`HC Met255 %Oxidation vs. Week_5C(…)`, `:HC Met255 %Oxidation`). | Shared, anonymized | **Yes**, with the owner's permission |
| Research data published openly (`files/public/`, 8 files, 1.7–77 KB; CC0 from Dryad and Harvard Dataverse, one CC BY 4.0 from figshare): JMP 13 in Japanese with Oneway analyses; JMP 12 and 11 with formulas and notes; JMP 10 with 21 hidden columns, value labels and Value Colors; JMP 11 with value labels and a Graph Builder; JMP 15 with the p-value format and text value labels; JMP 18 with its new date format. Chosen from the corpus runs (§9.4) for what the other fixtures lack. | CC0, one CC BY 4.0; credited in `CREDITS.md` §5 | **Yes** |
| Files made in JMP for us (edge cases: empty table, 0 rows, 1 col, all-missing, very wide, Unicode names, every display format) | — | **No.** We don't use JMP or ask customers to make files in it for us (§10). Edge cases come from files people already have, or from editing vendored files as `files/tests/` does. |
| JMP's own Samples/Data library (~600 files per version) | SAS copyright | **No.** Not used: it comes with a JMP install (§10). |
| Customer files | Customer confidential | **No.** Run locally only, with permission. |
| Files published openly by JMP users (papers, GitHub, Kaggle) | Varies | Only with a clear licence |

### 9.3 Golden-output tests (CI)

1. `tools/jmp2json.jl` loads each fixture with JMPReader.jl (pinned to the ported commit) and writes `files/golden/<name>.json`: column name, element type and values (date-times as ISO strings to the millisecond, dates as `yyyy-mm-dd`, durations as seconds, row states as `<marker> #RRGGBB`, missing as `null`). JSON rather than CSV keeps type, null and float precision unambiguous.
2. The golden files are **committed**, so CI doesn't need Julia. Regenerating them is a manual dev step after bumping the pinned JMPReader.jl version (setup in the script's header).
3. `src/tests/golden.ts` compares a parsed table with its golden file: names and order, row count, type compatibility per §5, exact equality for integers, strings and floats (the same IEEE bits), dates and date-times as written, null positions. It has no platform imports, so the package tests and the corpus script share it.
4. `jmp-reader-tests.ts` holds four categories, all marked `{node: true}` so `grok test` runs them in its headless Node pass: **header** (every fixture, plus non-JMP bytes, the old format, a truncated header), **reference** (every fixture against its golden file), **tables** (`runtests.jl` ported as named assertions) and **robustness** (the crafted files, plus a seeded fuzz of 420 mutations).
5. `jmp-reader-tests.ts` also checks the format details: fixed decimals (on `example1` with one descriptor changed, following JMPReader.jl's note), currency codes, and the hidden-column index.
6. `jmp-script-tests.ts`, also Node: **JSL** (calls, column references in their three forms, names with spaces, statements and sends, raw strings, syntax errors, a nesting limit, a fuzz of 2,000 truncated or mutated scripts), **saved scripts** (the `Source` script JMP wrote, files without scripts, the synthetic Graph Builder, the `New Script` form) and **analyses** (each platform in §11 Phase 3 in the form JMP saves it, with the Big Class columns; what is noted; what isn't converted and why).
7. `jmp-dataframe-tests.ts` runs in the browser: column types and nulls, integer sentinels, tags, geographic semTypes, row states as marker and colour columns with their codings, a saved analysis added to a table view as a bar chart, a Graph Builder trellis with its inner scatter plot, currency units and a hidden column hidden in the grid, the handler called as `JmpReader:importJmp`, and the preview built through `JmpReader:previewJmp` (disposed with no live u2 scopes left; listing a saved analysis that isn't converted). It reaches the preview through the platform, so the test bundle carries no u2 (the Node runtime has none).

8. **JMP: performance** (browser) times the 100 MB file of §8. It runs only with `grok test --benchmark`.

Result on the local stand (2026-10-28): Node pass 101 passed; browser pass 15 passed, the benchmark skipped (and passed when run with `--benchmark`).

### 9.4 Corpus run (local, not CI)

`npm run corpus -- <folder> [--golden <folder>]` (`tools/jmp-corpus.ts`, through `tsx`) reads every `.jmp` file under the folder with the package's parser and reports files read, files with undecoded columns, and failures grouped by reason; the saved scripts found, by JMP platform, with how many convert to viewers and why the rest don't; and row-state columns whose undecoded bytes are set, with the values seen (§4.3). With `--golden`, it also compares each file with JMPReader.jl's output for it (run `tools/jmp2json.jl` over the same folder first). Run it on any customer corpus we're given access to; where files can't leave the customer, use `jmp-survey.mjs` below. It also counts column properties by name, with any it can't read, and columns with value labels (§4.6). It hasn't been run on a customer corpus yet.

For files that can't leave their owner, `tools/jmp-survey.ts` is bundled into a single `jmp-survey.mjs` (command at the top of the file; no dependencies, Node 18 or later) and run where the files are: `node jmp-survey.mjs <folder> --summary-only`. It reports files read by JMP version; saved analyses by platform, with how many become viewers and why the rest don't; what each platform's conversion leaves out; and the JSL keywords each platform's scripts use, as a tree with counts. `--summary-only` prints no file paths, script names, column names or values, so the output can be sent back. A first run at a customer (2026-10-27) covered 4 files saved by JMP 15, with one Graph Builder.

A **public corpus** stands in until then: 14 `.jmp` files from GitHub repositories, mostly course and lab data saved by JMP 14–17 (`gm105/JMP`, `whippo/Whippo-NCE-2024`, `heyselbi/Marketing-Analytics-auditing`, `pineapplefart/S5_Hackathon_Supplementary`, `ArvindRajen/Salary_Prediction`, `BreakOnExecute/JMPReader`). They have no licence, so they stay outside the repo. Run on 2026-10-23: 13 of 14 read with every column decoded; the 14th is the older `(3,4,3)` format, refused as expected. 45 saved scripts were found and all read (Graph Builder, Distribution, Oneway, Bivariate, Contingency and Tabulate convert; Fit Model, DOE, Multivariate, Partition, Neural and screening don't). 87 column properties in 12 kinds and 2 value-labelled columns, all read. A larger **open corpus** came from Zenodo (2026-10-27): every `.jmp` file in its 32 records that have one, 224 files (505 MB) of published research data, all openly licensed (157 CC-BY 4.0, 57 CC0, 9 MIT, 1 GPL-3), saved by JMP 11–19, mostly 16. It stays outside the repo apart from the four CC0 fixtures in §9.2. 206 of 223 read with every column decoded; the other 17 are the old big-endian format (§4.1), from three records of 2007–08 data, refused with the old-format message. Of 333 saved scripts, 331 read once the parser took `;` statements inside arguments, matrix literals and `i++` (2 use names JSL can't parse bare). 160 are hand-written programs. Of the analyses, Tabulate, Oneway, Graph Builder, Distribution and Contingency convert (49); Fit Model (36), Fit Group (22), Fit Curve (14) and DOE (10) don't yet. In Graph Builder, Smoother is the most used element (8 of 10 scripts) and isn't shown yet. 3,406 column properties in 14 kinds (2,401 formulas) and 55 value-labelled columns, all read.

A **sweep of other open archives** (2026-10-28) added 72 files with open licences: Figshare (59, mostly CC BY 4.0), Harvard Dataverse (10, CC0) and OSF (4, CC BY). Dryad lists 158 more CC0 files not on Zenodo, but downloading them needs an API token; files without a licence were left out. All 72 read; 31 of the public files are now from JMP 17–19, and none of the new ones is the old format. Two display formats were new (`0x61`, `0x90`, §4.3) and a few JSL forms (`#`, `%`, `'` and dots in names, a column given a format as `:Age(10, 1)`, associative arrays `["en" => "…"]`). With them, all 305 files read with every column decoded, and all 1,054 saved scripts and every column property parse. Graph Builder across the public files: 93 scripts, 86 convert; Points (43), Line of Fit (32), Bar (25), Box Plot (19), Line (15), Smoother (10); Overlay (49), Group X (22), Color (20); nearly all set a size and a title, about 20 set axis limits. No file had row-state flags set.

**Acceptance for Phase 1:** ≥ 99% of the files in customer corpora open (counted by `jmp-survey.mjs` where the files are), with zero files throwing on supported versions.


### 9.5 Behaviour (BDD)

What a user does and sees is in `bdd/features/`, compiled by `libraries/bdd` into Playwright specs in `bdd/generated/` (committed; `grok-bdd compile --check` fails on a stale one). Per that library's rules, a feature claims the UI and leaves file contents to the package tests:

| Feature | Scenarios |
|---|---|
| `open.feature` | a file from disk opens as a 4-row table when dropped, through the Open local file icon, and through File > Open > File... |
| `preview.feature` | the Browse preview shows version and size, Open opens the table; the "Columns not imported" section lists the unreadable column, sits above the grid folded and open |
| `problem-files.feature` | a dropped file with an unreadable column opens without it and warns; a file that isn't JMP shows why and opens no table |
| `analyses.feature` | a dropped file with a saved Graph Builder opens with a bar chart split by its column and titled after the script; one with a Fit Model opens and names it in a warning; the preview lists both, with what each opens as or why not, shows the script behind the Script link, and Open adds the viewer; a row-state file opens with marker and colour columns |

The preview uses the package's own published files (`App Data/JmpReader/...`), so nothing is written to the server. `bdd.config.json` sets the layout and console-error audits to fail and the leak audit to report. Run from the package folder:

```bash
node ../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=http://localhost:8080 node ../../libraries/bdd/bin/grok-bdd.js run --reporter=list
```

Result on the local stand (2026-10-23): 12 passed, 0 audit findings, run with `--workers=1`. With the default workers, the first scenario on each worker sometimes times out on the stand's cold start (server unavailable at login, or the Files node slow to expand) before any package step runs.
---

## 10. Licensing and Attribution

- **JMPReader.jl is MIT** (© 2024 Jaakko Ruohio). A TypeScript port is a derivative work. Reproduce the MIT notice in `CREDITS.md` §1, as Minitab does for JSZip. Name the upstream commit, and add a header comment in each ported source file.
- The vendored test fixtures carry the same MIT notice.
- There are **no runtime dependencies** (`DecompressionStream` and `TextDecoder` are built into the browser). If `pako` is added as a fallback, it is MIT and gets its own entry.
- `tools/jmp2json.jl` is a dev-only tool. Julia and JMPReader.jl are not bundled.
- **No JMP at Datagrok** (decided 2026-10-27). No one at Datagrok installs or uses JMP for this package, and we don't ask customers to make files in JMP or compare our output with JMP's for us. The [JMP Software License Agreement](https://jmp.com/en/store/jmp-software-license-agreement) binds its customers and users not to "reverse engineer, reverse assemble or decompile the Software" (§4) or use it to "create or assist a third party in creating an offering that competes with JMP" (§1.4.3(b)). What we know of the format comes only from JMPReader.jl, jmp-convert and `.jmp` files people already have: public files, and customer files read where they are, with permission.
- **Name.** JMP® is a trademark of JMP Statistical Discovery LLC, a SAS subsidiary. The plugin uses it only to say which files it reads: it shows as "JMP Reader", not "JMP", uses no JMP logo or styling, and its README carries the trademark notice and says it isn't affiliated with or endorsed by JMP or SAS. The package is `@datagrok/jmp-reader` (`JmpReader` on the server, so the generated `JMPReader:` calls resolve); the folder stays `packages/Jmp`.
- *Legal check (see §12):* confirm that reading `.jmp` files this way is outside that agreement before merging to master or publishing.

---

## 11. Roadmap

### Phase 1 — Data bridge (MVP)

- [x] Scaffold `packages/Jmp` from Minitab (package.json, tsconfig, README, CREDITS)
- [x] `ByteReader` + metadata parser (§4.1–4.2)
- [x] Column decoder: numeric, byte ints, dates/date-times/durations, currency, geo, the three character layouts, gzip, row states (§4.3)
- [x] DataFrame conversion with tags (§5)
- [x] `fileHandler` + `fileViewer`, u2 preview and error views (§6)
- [x] Vendored fixtures, golden JSON, unit + golden tests, fuzz test, BDD features (§9)
- [ ] Corpus run on customer files, with `jmp-survey.mjs` where they are (first run, 2026-10-27: 4 files, all read)
- [x] A 100 MB file against the §8 target: 3.0–6.1 s in the browser, on a synthetic file (a real one would still be worth timing)
- [x] Help page: `.jmp` row in `help/access/files/supported-formats.md` (help submodule; not committed yet)
- [x] Format details: the currency code (confirmed), `jmp.hidden` (confirmed on the public files) and fixed-decimal formats (provisional, §4.3–4.4)
- [x] Row-state columns as marker and colour columns (§5); flags still to decode

**Exit criteria:** all golden tests pass (done), the corpus run passes at ≥ 99%, and the 100 MB file meets the target in §8.

### Phase 2 — Fidelity

- Row-state flags → Datagrok **selection** (selected), **filter** (excluded and hidden, as a saved filter state), and **colour coding** (row colours/markers, where a viewer supports them). Needs files with flags set, and the place JMP keeps the table's own row states (§4.3).
- [x] Value labels → a string column of the labels, values kept in a tag (§4.6, §5)
- [x] Column properties: units → `units`, notes → `description`, value order → category order, spec limits → `jmp.specLimits`, all kept as JSL text (§4.6)
- [x] Formula text → a `jmp.formula` tag
- [ ] Formulas as Datagrok calculated columns, where a translation exists (Datagrok formula syntax: no ternary, use `If()`)
- [ ] Value colours → categorical colour coding, once JMP's colour numbers are known (§4.4)
- Table variables and notes → DataFrame tags
- Older format `(3,4,3)` support, if the corpus shows meaningful volume

Every new decoding needs a self-made fixture, and ideally an upstream PR to JMPReader.jl so the reference stays ahead of us.

### Phase 3 — Saved analyses → live viewers

JMP tables carry **table scripts** (JSL) that re-launch analyses. We read a narrow, declarative subset (platform name plus role assignments) and map it to Datagrok viewers; nothing is ever evaluated. This is the same idea as re-creating Minitab graphs (see `MINITAB_PLUGIN.md` §7.4).

**Built (2026-10-22), against JSL text; the file side is provisional (§4.5):**

| JMP platform | Datagrok | What doesn't come across (named in the preview) |
|---|---|---|
| Distribution | Histogram per continuous column; bar chart of counts per Nominal/Ordinal or text column | summary statistics, quantiles, fitted distributions |
| Bivariate | Scatter plot; Fit Line or Fit Polynomial(1) → regression line | other fits (spline, polynomial > 1, density ellipse) |
| Oneway | Box plot by group | means, comparisons, tests |
| Contingency | Bar chart, X split and Y stacked, counts | mosaic, tests |
| Fit Y by X | By column types, as JMP does: bivariate, oneway or contingency, per Y/X pair | logistic (categorical Y by continuous X): not converted |
| Graph Builder: Points, Line of Fit, Smoother, Line, Bar, Box Plot, Histogram, Pie | Scatter (colour by Color/Overlay, size by Size, Line of Fit as regression line, one per category when the Overlay is categorical and also the colour), line chart, bar chart (Summary Statistic as the aggregation), box plot, histogram, pie chart. With Group X, Group Y or Wrap: a trellis plot of that viewer (Group X and Wrap across, Group Y down; Wrap tiled) | Page; a second X or Y; further panels and elements; the smoother; a fit per group of a numeric Overlay, or of an Overlay that is not the colour |
| Tabulate | Pivot table: row-table grouping → rows, column-table grouping → pivot, analysis columns × statistics → aggregates | statistics without a Datagrok aggregate (shown as mean) |
| Control Chart Builder, Fit Model, DOE platforms, other Graph Builder elements | Not converted: listed with the reason | |

Every launch also notes a `By()`, `Where()`, `Weight` or `Freq` it doesn't apply. Viewers are titled after the script. A column a script names that isn't in the imported table makes that launch not convert, with the column named. The `Source` script (the import that made the table) isn't an analysis and isn't listed.

The preview shows each script's text behind a **Script** link.

**Next:** confirm §4.5 on customer files, through `jmp-survey.mjs`; a control chart once a Datagrok control-chart viewer exists.

### Phase 4 — Portal ingestion

The bigger payoff: **backfill historical JMP files into shared, searchable Datagrok spaces.**

- A bulk import job takes a file-share connection or a folder and calls `readJmp` on each `.jmp` file
- Each file → a Datagrok project/table with metadata indexed for search (source path, JMP version, save time, column names, table notes, saved analysis names)
- Runs on the server (Node or a Datagrok server function) so it isn't limited by browser memory; same TypeScript parser
- Report per file: imported, imported with undecoded columns, or failed (with reason)
- Feeds the Project Portal: imported tables and re-created viewers become portal artefacts, with source provenance kept

---

## 12. Open Questions

1. **Package boundary.** Keep `Jmp` standalone, or fold it into a combined "statistics formats" package with Minitab? *Recommendation:* standalone for v1, then reconsider once both exist.
2. **Format version coverage.** Which JMP versions write `(7,4,0x00010001)`, and how common are `(3,4,3)` files in real customer archives? We need a customer corpus to answer this.
3. **Test corpus.** Which customers will share representative (non-confidential or anonymised) `.jmp` files? Who makes the request?
4. **Legal.** Confirm that porting an MIT-licensed independent reader, and reading the rest of the format from `.jmp` files without JMP, is outside the JMP Software License Agreement (§10), before merging to master or publishing. The code is already on a public branch. Datagrok doesn't install JMP (§10).
5. **Server-side runtime for Phase 4.** Run the TS parser in Node inside a Datagrok server function, or rely on browser workers plus chunked upload?
6. **Other JMP file types.** Is there demand for `.jmpprj` (projects that bundle tables + reports) or `.jrn` journals? A JMP project probably holds the decision record that Phase 4 is after.
7. **Write support.** Confirm the non-goal. Is there any customer flow that needs Datagrok → `.jmp`, rather than CSV/Excel?
8. **Preview on double-click.** A `.jmp` double-clicked in the Browse tree previews (Open is one more click), as other formats with a package viewer do, while a file from disk opens straight as a table. Keep, or import on double-click too? The shell decides by whether a viewer exists, so importing on double-click means dropping the preview.
9. **Telling disk from server.** The viewer takes a file without a connection to be one from disk (§3). That matches the shell's code today; a platform flag saying "open, don't preview" would be sturdier.
10. **Analyses through the file handler.** `importJmp` returns tables only, so with the core fix (a viewer that declines falls through to the handler), a file from disk would open without its analyses. Either the handler adds them once its view opens, or the shell lets a handler return a view.

---

## References

- JMPReader.jl: https://github.com/jaakkor2/JMPReader.jl (MIT) · [announcement on Julia Discourse](https://discourse.julialang.org/t/ann-jmpreader-jl-reader-for-jmp-files-in-julia/120710)
- JMP's `jmp` Python package (runs only inside JMP, so it isn't usable as a reader): https://www.jmp.com/support/help/en/19.1/jmp/the-jmp-package-for-python.shtml
- Minitab plugin (structural template): [`../../Minitab/docs/MINITAB_PLUGIN.md`](../../Minitab/docs/MINITAB_PLUGIN.md), [`../../Minitab/src/package.ts`](../../Minitab/src/package.ts), [`../../Minitab/CREDITS.md`](../../Minitab/CREDITS.md)
