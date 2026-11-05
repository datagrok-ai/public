@browse @realizes:views.browse
Feature: The file formats of the demo files preview and open from Browse
  A file of Browse > Files > Demo (System:DemoFiles) shows in a Browse preview on a single click and
  opens on a double click, whatever its format, with nothing logged. A table format previews as a
  table view that is not yet in the workspace, and the double click puts its table there with every
  row; a document previews in a view of its own, named after the file, that shows the file's content,
  and the double click keeps it. Translated from TestTrack FileFormats/supported-formats.md (and its
  spec), over the formats of that case that the demo files hold: CSV, JSON, XML, HTML, Markdown,
  PDF, SDF, Excel and ZIP. A format a package reads is gated on that package; a workbook opens a
  table per sheet, the demo ZIP a table per workbook it holds (PowerPack reads them), and one of
  them is claimed. The PDF preview was a view with no name, so its tab had no title (FileEditors
  previewPdf, fixed 2026-10-13, the DOCX preview with it).

  Not translated: the case's folder "My Files / all formats" is a hand-made fixture of 30 files that
  no stand provides; the formats the demo files do not hold (EDF, Feather, HDF5, NetCDF, Parquet,
  SAS, SQLite, MAT, IVP, ...) are the package tests' of the packages that read them, over their own
  test files. KML and KMZ preview on a map whose markers and tiles come from outside web services.
  What a reader makes of a file's content (the columns, the types, the values) is the reader's
  package test, not this feature's: the claims are the preview, the opening and the row count.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded

  Scenario: A CSV file previews as a table and opens into the workspace
    When user clicks on Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    And no table should be left in the workspace
    When user double-clicks on Files---Demo---cars.csv tree node inside browse tree
    Then table "cars" should have 30 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: A <format> file a package reads previews as a table and opens into the workspace
    Given the "<package>" package is installed
    And <parent> tree node inside browse tree is expanded
    And <folder> tree node inside browse tree is expanded
    When user clicks on <file> tree node inside browse tree
    Then the "<view>" view should be current
    And no table should be left in the workspace
    When user double-clicks on <file> tree node inside browse tree
    Then table "<table>" should have <rows> rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | format | package   | parent                   | folder                    | file                                         | view         | table        | rows |
      | SDF    | Chem      | Files---Demo---chem      | Files---Demo---chem---sdf | Files---Demo---chem---sdf---vitamin-d-2d.sdf | vitamin-d-2d | vitamin-d-2d | 1    |
      | Excel  | PowerPack | Files---Demo---northwind | Files---Demo---northwind  | Files---Demo---northwind---Northwind.xlsx    | Suppliers    | Orders       | 830  |
      | ZIP    | PowerPack | Files---Demo             | Files---Demo              | Files---Demo---chika-mexicana.zip            | Detailed Tickets from 2015-07-31 to 2015-08-07.xlsx - Worksheet | Detailed Tickets from 2015-05-08 to 2015-05-15.xlsx - Worksheet | 1859 |

  Scenario Outline: A <format> file previews as a document and stays open
    Given <parent> tree node inside browse tree is expanded
    And <folder> tree node inside browse tree is expanded
    When user clicks on <file> tree node inside browse tree
    Then the "<view>" view should be current
    And the current view should show the text "<text>"
    When user double-clicks on <file> tree node inside browse tree
    Then the "<view>" view should be current
    And the current view should show the text "<text>"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | format | parent                 | folder                             | file                                                          | view                     | text                         |
      | JSON   | Files---Demo---kaggle  | Files---Demo---kaggle              | Files---Demo---kaggle---insurance-metadata.json               | insurance-metadata.json  | @context                     |
      | XML    | Files---Demo---bio     | Files---Demo---bio---ngl-formats   | Files---Demo---bio---ngl-formats---3pqr_validation.xml        | 3pqr_validation.xml      | wwPDB-validation-information |
      | HTML   | Files---Demo---forms   | Files---Demo---forms               | Files---Demo---forms---tennis_players.html                    | tennis_players.html      | display: flex                |

  Scenario: A Markdown file previews rendered and stays open
    Given the "PowerPack" package is installed
    And Files---Demo---chem tree node inside browse tree is expanded
    And Files---Demo---chem---protac tree node inside browse tree is expanded
    When user clicks on Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree
    Then the "protac_degraders" view should be current
    And the current view should show the text "2,792 PROTAC degraders"
    When user double-clicks on Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree
    Then the "protac_degraders" view should be current
    And the current view should show the text "2,792 PROTAC degraders"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A PDF file previews its pages and stays open
    Given the "FileEditors" package is installed
    And Files---Demo---texts tree node inside browse tree is expanded
    And Files---Demo---texts---pdf tree node inside browse tree is expanded
    When user clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree
    Then the "research-paper.pdf" view should be current
    And pdf preview page should be visible
    When user double-clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree
    Then the "research-paper.pdf" view should be current
    And pdf preview page should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
