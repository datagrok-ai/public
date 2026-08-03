# EpsViewer

`EpsViewer` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai)
platform that previews Encapsulated PostScript (`.eps`) files in the file browser.

The viewer handles two common flavors of EPS:

* **ASCII EPS** — rendered by parsing a subset of PostScript onto an HTML canvas via the
  vendored [UDOC](https://github.com/photopea/UDOC) interpreter.
* **DOS EPS binary** (Adobe Illustrator and similar) — the embedded TIFF preview is extracted from the binary header and
  decoded with [utif2](https://www.npmjs.com/package/utif2).

Complex, real-world Illustrator files typically fall into the second category, so their embedded preview is used instead
of attempting to interpret the full PostScript program.
