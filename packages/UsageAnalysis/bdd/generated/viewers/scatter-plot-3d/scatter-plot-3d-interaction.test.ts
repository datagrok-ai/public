/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot-3d/scatter-plot-3d-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot3d]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, currentRowValue, filterPasses, filterTo, mouseOverRowIs, noneSelected, onlyOfSelected, resetFilter, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickAreaHolding, closeContextMenu, dragAcrossArea, hoverArea, noErrors, oneTooltip, pickFromContextMenu, propertyShouldBe, readingDiffers, readingIs, setProperty, takeSnapshot, tooltipValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("3D scatter plot interaction", () => {
  const session = feature(test, "features/viewers/scatter-plot-3d/scatter-plot-3d-interaction.feature", import.meta.url);
  test("3D scatter plot interaction", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
    await session.step(22, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(23, "And \"USUBJID\" of the current row should be \"X0273T21000300003\"", () => currentRowValue(page, "USUBJID", "X0273T21000300003"));
    await run.scenario("Hovering a point shows its row in a tooltip", async () => {
      await session.step(26, "When user hovers over the \"point\" area of 3d scatter plot viewer", () => hoverArea(page, "point", el("3d scatter plot viewer")));
      await session.step(27, "Then the mouse-over row of the table should be 983", () => mouseOverRowIs(page, 983));
      await session.step(28, "And exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(29, "And the tooltip should show \"USUBJID\" as \"X0273T51060300003\"", () => tooltipValue(page, "USUBJID", "X0273T51060300003"));
      await session.step(30, "And the tooltip should show \"AGE\" as \"50\"", () => tooltipValue(page, "AGE", "50"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl-click toggles a point's selection and Ctrl+Shift-click removes it", async () => {
      await session.step(34, "Given user clears the row selection", () => clearSelection(page));
      await session.step(35, "When user clicks on the \"point\" area of 3d scatter plot viewer holding Control", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Control"));
      await session.step(36, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(37, "And only rows where \"USUBJID\" is \"X0273T51060300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T51060300003"));
      await session.step(38, "And \"USUBJID\" of the current row should be \"X0273T21000300003\"", () => currentRowValue(page, "USUBJID", "X0273T21000300003"));
      await session.step(39, "When user clicks on the \"point\" area of 3d scatter plot viewer holding Control", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Control"));
      await session.step(40, "Then no rows should be selected", () => noneSelected(page));
      await session.step(41, "When user clicks on the \"point\" area of 3d scatter plot viewer holding Shift", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Shift"));
      await session.step(42, "Then only rows where \"USUBJID\" is \"X0273T51060300003\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T51060300003"));
      await session.step(43, "When user clicks on the \"point\" area of 3d scatter plot viewer holding Control+Shift", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Control+Shift"));
      await session.step(44, "Then no rows should be selected", () => noneSelected(page));
      await session.step(45, "When user clicks on the \"point\" area of 3d scatter plot viewer holding Control+Shift", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Control+Shift"));
      await session.step(46, "Then no rows should be selected", () => noneSelected(page));
      await session.step(47, "And \"USUBJID\" of the current row should be \"X0273T21000300003\"", () => currentRowValue(page, "USUBJID", "X0273T21000300003"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Filter > Show Filtered Out Points in the context menu brings the filtered-out points back", async () => {
      await session.step(51, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(52, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(53, "And the \"rows shown\" reading of 3d scatter plot viewer should be 499", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 499));
      await session.step(54, "When user picks \"Filter > Show Filtered Out Points\" from the context menu of 3d scatter plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Points", el("3d scatter plot viewer")));
      await session.step(55, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(56, "Then \"Show Filtered Out Points\" property of 3d scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Filtered Out Points", el("3d scatter plot viewer"), "true"));
      await session.step(57, "And the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(58, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(59, "When user picks \"Filter > Show Filtered Out Points\" from the context menu of 3d scatter plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Points", el("3d scatter plot viewer")));
      await session.step(60, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(61, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 499", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 499));
      await session.step(62, "When user resets the filter", () => resetFilter(page));
      await session.step(63, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click on empty space brings the camera home", async () => {
      await session.step(67, "When user drags across the \"view\" area of 3d scatter plot viewer", () => dragAcrossArea(page, "view", el("3d scatter plot viewer")));
      await session.step(68, "Then the \"camera x\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "camera x", el("3d scatter plot viewer")));
      await session.step(69, "When user double-clicks on empty plot space of 3d scatter plot viewer", () => doubleClickEmptySpace(page, el("3d scatter plot viewer")));
      await session.step(70, "Then the \"camera x\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "camera x", el("3d scatter plot viewer"), 0));
      await session.step(71, "And the \"camera y\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "camera y", el("3d scatter plot viewer"), 0));
      await session.step(72, "And the \"camera distance\" reading of 3d scatter plot viewer should be 4", () => readingIs(page, "camera distance", el("3d scatter plot viewer"), 4));
      await session.step(73, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected draws only the selected rows", async () => {
      await session.step(77, "Given user clears the row selection", () => clearSelection(page));
      await session.step(78, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(79, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(80, "When user sets \"Row Source\" property of 3d scatter plot viewer to \"Selected\"", () => setProperty(page, "Row Source", el("3d scatter plot viewer"), "Selected"));
      await session.step(81, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 499", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 499));
      await session.step(82, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
      await session.step(83, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(84, "And the \"rows shown\" reading of 3d scatter plot viewer should be 373", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 373));
      await session.step(85, "When user clears the row selection", () => clearSelection(page));
      await session.step(86, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 0));
      await session.step(87, "When user sets \"Row Source\" property of 3d scatter plot viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("3d scatter plot viewer"), "Filtered"));
      await session.step(88, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A string column on an axis draws every row", async () => {
      await session.step(92, "When user takes a snapshot of 3d scatter plot viewer", () => takeSnapshot(page, el("3d scatter plot viewer")));
      await session.step(93, "And user picks \"RACE\" in the \"X\" column selector of 3d scatter plot viewer", () => pickInColumnSelector(page, "RACE", "X", el("3d scatter plot viewer")));
      await session.step(94, "Then X column input in 3d scatter plot viewer should contain text \"RACE\"", () => shouldContainText(page, el("X column input in 3d scatter plot viewer"), "RACE"));
      await session.step(95, "And the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(96, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(97, "When user picks \"AGE\" in the \"X\" column selector of 3d scatter plot viewer", () => pickInColumnSelector(page, "AGE", "X", el("3d scatter plot viewer")));
      await session.step(98, "Then X column input in 3d scatter plot viewer should contain text \"AGE\"", () => shouldContainText(page, el("X column input in 3d scatter plot viewer"), "AGE"));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A boolean label column on 5850 rows draws without freezing (GROK-19048)", async () => {
      await session.step(102, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(103, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
      await session.step(104, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 5097", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 5097));
      await session.step(105, "When user sets \"Label\" property of 3d scatter plot viewer to \"CONTROL\"", () => setProperty(page, "Label", el("3d scatter plot viewer"), "CONTROL"));
      await session.step(106, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A SMILES column on an axis does not freeze the plot (public#398)", async () => {
      await session.step(110, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
      await session.step(111, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(112, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
      await session.step(113, "Then \"canonical_smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles", "Molecule"));
      await session.step(114, "When user takes a snapshot of 3d scatter plot viewer", () => takeSnapshot(page, el("3d scatter plot viewer")));
      await session.step(115, "And user picks \"canonical_smiles\" in the \"X\" column selector of 3d scatter plot viewer", () => pickInColumnSelector(page, "canonical_smiles", "X", el("3d scatter plot viewer")));
      await session.step(116, "Then X column input in 3d scatter plot viewer should contain text \"canonical_smiles\"", () => shouldContainText(page, el("X column input in 3d scatter plot viewer"), "canonical_smiles"));
      await session.step(117, "And the \"rows shown\" reading of 3d scatter plot viewer should be 1000", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 1000));
      await session.step(118, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(119, "When user takes a snapshot of 3d scatter plot viewer", () => takeSnapshot(page, el("3d scatter plot viewer")));
      await session.step(120, "And user drags across the \"view\" area of 3d scatter plot viewer", () => dragAcrossArea(page, "view", el("3d scatter plot viewer")));
      await session.step(121, "Then the \"camera x\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "camera x", el("3d scatter plot viewer")));
      await session.step(122, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
