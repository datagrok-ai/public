# Ketcher

Adds [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html), the open-source 2D
molecule editor by EPAM, as an optional molecular sketcher in the Datagrok platform. Once the
package is installed, "Ketcher" appears in every sketcher picker across the platform — grid cell
editor, dialog inputs, substructure filter, and so on.

## Usage

Pick **Ketcher** from the sketcher dropdown. The selection is
per-user and persisted via platform preferences.

![Ketcher](./images/ketcher.gif)

## Known limitations

* A recursive SMARTS (`$(...)`) is not held: Ketcher (Indigo's SMARTS reader) drops it, so `[#6;$([#6]=[#8])]-[#7]`
  set on Ketcher is drawn, and searched, as `[#6]-[#7]` (on SPGI, 100 rows instead of 88).
