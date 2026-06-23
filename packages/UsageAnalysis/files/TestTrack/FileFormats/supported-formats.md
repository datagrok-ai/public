Verify that Datagrok previews and opens every file in Browse > Files > My Files > all_formats.

**Prerequisite**: folder `My Files / all_formats` contains one file per format —
.bmp, .csv, .edf, .fasta, .feather, .geojson, .gz, .h5, .html, .ipynb, .ivp, .json,
.kml, .kmz, .kxl, .mat, .md, .nc, .parquet, .pdf, .rda, .rds, .sas7bdat, .sdf,
.sqlite, .tar, .topojson, .xlsx, .xml, .zip

1. Open Browse > Files > My Files > all_formats.

2. For each file in the folder — **single-click** (preview):
   the context panel on the right must show a non-blank preview.
   No error dialog, no blank/loading panel stuck.

3. For each file — **double-click** (open):
   the file must open in the main workspace — as a table, a viewer, or a document.
   No error dialog, no "Cannot open file" message.
   Close the view before moving to the next file.

4. Open the browser console (F12). No red errors related to file parsing
   should appear for any of the 30 formats.

---
{
  "order": 1
}
