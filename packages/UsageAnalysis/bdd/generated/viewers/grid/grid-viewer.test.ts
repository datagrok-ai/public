/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaAtLeastTall, boundTable, closeContextMenu, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, propertyShouldNotBe, rightClickArea, setProperty, showsRows, viewerAdded, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A second grid as a viewer, and column tooltips", () => {
  const session = feature(test, "features/viewers/grid/grid-viewer.feature", import.meta.url);
  test("A second grid as a viewer, and column tooltips", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(14, "And the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
    await run.scenario("A second grid shows the same table and takes its own row height", async () => {
      await session.step(17, "When user adds a grid viewer", () => addViewer(page, "grid"));
      await session.step(18, "Then grid viewer should be added to the open tableview", () => viewerAdded(page, "grid"));
      await session.step(19, "And the open tableview should have 2 grid viewers", () => viewerCount(page, 2, "grid"));
      await session.step(20, "And second grid viewer should show 1000 rows", () => showsRows(page, el("second grid viewer"), 1000));
      await session.step(21, "When user sets \"Row Height\" property of second grid viewer to \"40\"", () => setProperty(page, "Row Height", el("second grid viewer"), "40"));
      await session.step(22, "Then the \"cell 1 of USUBJID\" area of second grid viewer should be at least 34 pixels tall", () => areaAtLeastTall(page, "cell 1 of USUBJID", el("second grid viewer"), 34));
      await session.step(23, "And \"Row Height\" property of grid should not be \"40\"", () => propertyShouldNotBe(page, "Row Height", el("grid"), "40"));
      await session.step(24, "When user sets \"Show Column Labels\" property of second grid viewer to \"false\"", () => setProperty(page, "Show Column Labels", el("second grid viewer"), "false"));
      await session.step(25, "Then second grid viewer should not have a \"header AGE\" area", () => hasNoArea(page, el("second grid viewer"), "header AGE"));
      await session.step(26, "And grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
      await session.step(27, "When user sets \"Show Column Labels\" property of second grid viewer to \"true\"", () => setProperty(page, "Show Column Labels", el("second grid viewer"), "true"));
      await session.step(28, "Then second grid viewer should have a \"header AGE\" area", () => hasArea(page, el("second grid viewer"), "header AGE"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second grid rebinds to another table and closes without taking the view's grid", async () => {
      await session.step(32, "When user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(33, "Then the table should have 100 rows", () => rowCount(page, 100));
      await session.step(34, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(35, "And user sets \"Table\" property of second grid viewer to \"spgi-100\"", () => setProperty(page, "Table", el("second grid viewer"), "spgi-100"));
      await session.step(36, "Then second grid viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("second grid viewer"), "spgi-100"));
      await session.step(37, "And second grid viewer should show 100 rows", () => showsRows(page, el("second grid viewer"), 100));
      await session.step(38, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(39, "When user clicks on close icon of second grid viewer", () => clickOn(page, el("close icon of second grid viewer")));
      await session.step(40, "Then the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
      await session.step(41, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(42, "And grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
      await session.step(43, "When user closes all views", () => closeAllViews(page));
      await session.step(44, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column tooltip menu marks the setting that is on", async () => {
      await session.step(47, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(48, "When user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(49, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(50, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(51, "Then \"Default\" menu item in context menu should be visible", () => shouldBe(page, el("\"Default\" menu item in context menu"), "visible"));
      await session.step(52, "And \"Form\" menu item in context menu should be visible", () => shouldBe(page, el("\"Form\" menu item in context menu"), "visible"));
      await session.step(53, "And \"Columns\" menu item in context menu should be visible", () => shouldBe(page, el("\"Columns\" menu item in context menu"), "visible"));
      await session.step(54, "And \"None\" menu item in context menu should be visible", () => shouldBe(page, el("\"None\" menu item in context menu"), "visible"));
      await session.step(55, "And \"Default\" menu item in context menu should be selected", () => shouldBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(56, "And \"None\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(57, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(58, "And user picks \"Tooltip > Current Column > None\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Tooltip > Current Column > None", "header AGE", el("grid")));
      await session.step(59, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(60, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(61, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(62, "Then \"None\" menu item in context menu should be selected", () => shouldBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(63, "And \"Default\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(64, "When user clicks on \"Default\" menu item in context menu", () => clickOn(page, el("\"Default\" menu item in context menu")));
      await session.step(65, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(66, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(67, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(68, "Then \"Default\" menu item in context menu should be selected", () => shouldBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(69, "And \"None\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(70, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(71, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
