# Plate files — what the features could not say yet

## Plates

### Needs a signal (Plates)

- **A well hovered or clicked on the plate.** The Summary heat map is a grid with custom-drawn wells
  (`plate-widget.ts`); the plate widget reports no hit area per well, so neither a well's tooltip nor
  `onWellClick` (the well details panel) can be aimed at. A status provider on `PlateWidget`
  (`DG.Widget.addStatusProvider`, areas `well A01` … in px of its root) would let a feature say:

  ```gherkin
  Scenario: A click on a well shows its values
    When user clicks on the "well A01" area of plate widget
    Then the current view should show the text "Channel A"
    And no errors should have been logged
  ```

### Questions (Plates)

- **Open local file on a plate .txt shows the file viewer, not the importer.** The view that comes up
  is named after the file (`previewPlate`, the `fileViewer`), not "Plate" as `importPlate` (the
  `fileHandler` for the same extension) names its view. Both show the same plate; which one is meant
  for Open local file?
- **Excel plates** (`previewPlateXlsx`, `importPlateXlsx`) are not claimed: there is no Excel plate on
  the local stand (`System:DemoFiles/hts/xlsx_plates` is missing), and a binary fixture in the project is
  a decision for the lead.
