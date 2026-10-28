# JMP Reader

JMP is widely used in pharma R&D, CMC and manufacturing for exploratory statistics, DOE and quality work.
A lot of that work lives in `.jmp` files on scientists' drives and team shares: hard to find, hard to share,
and invisible to anyone without a JMP license. Today, the only way to get them into Datagrok is to open each
one in JMP and export it to CSV, which loses types and metadata and doesn't scale to hundreds of files.

We want Datagrok to open `.jmp` files directly, in the browser, with no JMP installation and no server
component.

## Data

Opening a `.jmp` file should give you the table as JMP shows it:
- **Columns**: numbers, integers, text, dates, times, durations, currencies, geographic coordinates
- **Missing values**: as nulls, not as sentinels
- **Row states**: markers and colors as columns
- **Column properties**: notes, units, formulas, value order, value labels
- **Table metadata**: JMP version, save time, hidden columns

A column we can't read should come back empty with a warning, never cost the rest of the table.
A file we can't read should say why.

## Saved analyses

A JMP table often carries the analyses done on it, as saved scripts (Distribution, Bivariate, Oneway,
Graph Builder, Tabulate, ...). We want them to open as the closest Datagrok viewers, next to the table.
What doesn't come across (a fit, a test, a font) should be named, not silently dropped.

The preview, shown when you click a file in the Browse tree, lists the columns and saved analyses before
you open anything.

## Bulk import

The bigger payoff is backfilling: point Datagrok at a file share, import every `.jmp` file into a searchable
space, keep where each one came from, and report what didn't import and why. Same parser, run on the server.

## Constraints

- **Read only.** No writing `.jmp` files; CSV and Excel cover the way back.
- **No JMP.** We don't install or use JMP, and don't ask customers to make files in it for us.
  What we know of the format comes from open-source readers (JMPReader.jl, jmp-convert) and from
  `.jmp` files people already have.
- **Real files decide.** Each format detail gets a test file. Where files can't leave a customer,
  `jmp-survey.mjs` runs where they are and reports counts only.

See [JMP_PLUGIN.md](JMP_PLUGIN.md) for the format, the design and the status.
