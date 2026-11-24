# File Editors changelog

## v.next

* PDF, DOCX, RTF: The preview view is named after the file (its tab had no title)
* RTF: Fixed the preview failing to load its renderer (the bundle was dropped in the rspack migration)
* DOCX, RTF: The preview scrolls as one document instead of squeezing each page into its own box
* DOC, ODT: Added previews (ODT with formatting, DOC as plain text)
* DOCX, RTF, ODT: A document can no longer run script in the preview (DOCX symbol markup, javascript: links), and DOCX styles no longer leak into the app

## 1.3.2 (2026-07-12)

* GROK-18695: Dependency security updates — upgraded pdfjs-dist 2.16 → 4.10 (fixes arbitrary JavaScript execution upon opening a malicious PDF)

## 1.3.1 (2025-07-30)

Added support of TEX files

## 1.0.0 (2023-11-13)

Released package

### Features

* Added support for PDF files (open and preview)
* Added support for DOCX files (open and preview)
