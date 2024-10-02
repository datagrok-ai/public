### Shape map viewer 
Shows a map that is applicable for the specified dataset. Typically, it would represent a geographical area (countries,
states, counties, etc), but also it can show an arbitrary shapes (such as a store floor plan, brain regions, or EEG
electrodes).

1. Go to Browse > Files > Demo > geo. Open "ua_population" dataset. 
2. Click the **Add viewer** icon, find **Shape map** viewer and press on it. 
*Expected result: **Shape map** viewer opens without errors in console. Map of Ukraine regions is drawn on viewer.
3. Close previously opened Shape map viewer. On the **Viewers tab**, click **Shape map** icon. 
*Expected result: **Shape map** viewer opens without errors in console. Map of Ukraine regions is drawn on viewer.
4. Interact with all elements directly on the viewer.
* Click on region on map. Expected Result: Row (s) selected.
* Select not matching rows (from viewer context menu or from "hamburger" menu). Expected Result: One row ("Kiev City 1") selected.
* Change display color (from viewer context menu or from "hamburger" menu). Expected Result: Color changes according to the selected palette.
5. On the **Shape map** viewer, click the **Gear** icon. The **Property Pane** opens.
6. Modify various properties in **Property Pane**:
* Expected Result: Changes to the properties should be reflected in the Shape map viewer without errors, and the viewer should be updated accordingly.

---
{
  "order": 26
}