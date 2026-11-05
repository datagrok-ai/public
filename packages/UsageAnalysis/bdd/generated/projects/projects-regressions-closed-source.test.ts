/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-closed-source.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#531]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, noProjectOnServer, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addToRow, clearSavedParameters} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {clickArea, closeContextMenu, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a pivot saved after its source table is closed", () => {
  const session = feature(test, "features/projects/projects-regressions-closed-source.feature", import.meta.url);
  test("Data sync is not offered for a pivot whose source table is closed", {tag: ["@realizes:views.projects", "@known-failure:49", "@realizes:core#531"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await knownFailure(async () => {
      await session.step(22, "Given no project named \"BDDRegClosedSrc{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegClosedSrc{time}")));
      await session.step(23, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
      await session.step(24, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(25, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(26, "And Files---Demo---northwind tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---northwind tree node inside browse tree")));
      await session.step(27, "When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---northwind---orders.csv tree node inside browse tree")));
      await session.step(28, "Then the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
      await session.step(29, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
      await session.step(30, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(31, "When user clicks on \"pivot table\" icon in toolbox", () => clickOn(page, el("\"pivot table\" icon in toolbox")));
      await session.step(32, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
      await session.step(33, "When user adds \"ShipCountry\" to the \"group by\" row of pivot table viewer", () => addToRow(page, "ShipCountry", "group by"));
      await session.step(34, "And user clicks on the \"remove group by chip ShipName\" area of pivot table viewer", () => clickArea(page, "remove group by chip ShipName", el("pivot table viewer")));
      await session.step(35, "And user adds \"ShipVia\" to the \"pivot\" row of pivot table viewer", () => addToRow(page, "ShipVia", "pivot"));
      await session.step(36, "And user clicks on the \"remove pivot chip CustomerID\" area of pivot table viewer", () => clickArea(page, "remove pivot chip CustomerID", el("pivot table viewer")));
      await session.step(37, "And user picks \"Aggregation > count\" from the context menu of the \"aggregate chip avg(OrderID)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > count", "aggregate chip avg(OrderID)", el("pivot table viewer")));
      await session.step(38, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(39, "Then the \"group by\" reading of pivot table viewer should be \"ShipCountry\"", () => readingReads(page, "group by", el("pivot table viewer"), "ShipCountry"));
      await session.step(40, "And the \"pivot\" reading of pivot table viewer should be \"ShipVia\"", () => readingReads(page, "pivot", el("pivot table viewer"), "ShipVia"));
      await session.step(41, "And the \"aggregate\" reading of pivot table viewer should be \"count(OrderID)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "count(OrderID)"));
      await session.step(42, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(43, "Then the \"orders aggregation\" view should be current", () => viewIsCurrent(page, "orders aggregation"));
      await session.step(44, "And the table should have 40 rows", () => rowCount(page, 40));
      await session.step(45, "When user closes \"orders\" view tab", () => close(page, el("\"orders\" view tab")));
      await session.step(46, "Then the \"orders aggregation\" view should be current", () => viewIsCurrent(page, "orders aggregation"));
      await session.step(47, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(48, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(49, "And Data sync switch in \"orders aggregation\" project table in \"Save project\" dialog should be disabled", () => shouldBe(page, el("Data sync switch in \"orders aggregation\" project table in \"Save project\" dialog"), "disabled"));
      await session.step(50, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(51, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    }, 49);
  });
});
