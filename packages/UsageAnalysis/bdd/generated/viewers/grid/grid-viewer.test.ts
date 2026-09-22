/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-viewer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
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
import {closeAllViews, openDataset, openDatasetRowsAs, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaAtLeastTall, boundTable, closeContextMenu, hasArea, hasNoArea, hoverArea, noErrors, pickFromAreaContextMenu, pointerAway, propertyShouldBe, propertyShouldNotBe, readingAsRemembered, readingDoesNotRead, readingNotAsRemembered, rememberReading, rightClickArea, setProperty, showsRows, tooltipColumns, tooltipValue, viewerAdded, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A second grid as a viewer, and column tooltips", () => {
  const session = feature(test, "features/viewers/grid/grid-viewer.feature", import.meta.url);
  test("A second grid as a viewer, and column tooltips", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(37, "And the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
    await run.scenario("A second grid shows the same table and takes its own row height", async () => {
      await session.step(40, "When user adds a grid viewer", () => addViewer(page, "grid"));
      await session.step(41, "Then grid viewer should be added to the open tableview", () => viewerAdded(page, "grid"));
      await session.step(42, "And the open tableview should have 2 grid viewers", () => viewerCount(page, 2, "grid"));
      await session.step(43, "And second grid viewer should show 1000 rows", () => showsRows(page, el("second grid viewer"), 1000));
      await session.step(44, "When user sets \"Row Height\" property of second grid viewer to \"40\"", () => setProperty(page, "Row Height", el("second grid viewer"), "40"));
      await session.step(45, "Then the \"cell 1 of USUBJID\" area of second grid viewer should be at least 34 pixels tall", () => areaAtLeastTall(page, "cell 1 of USUBJID", el("second grid viewer"), 34));
      await session.step(46, "And \"Row Height\" property of grid should not be \"40\"", () => propertyShouldNotBe(page, "Row Height", el("grid"), "40"));
      await session.step(47, "When user sets \"Show Column Labels\" property of second grid viewer to \"false\"", () => setProperty(page, "Show Column Labels", el("second grid viewer"), "false"));
      await session.step(48, "Then second grid viewer should not have a \"header AGE\" area", () => hasNoArea(page, el("second grid viewer"), "header AGE"));
      await session.step(49, "And grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
      await session.step(50, "When user sets \"Show Column Labels\" property of second grid viewer to \"true\"", () => setProperty(page, "Show Column Labels", el("second grid viewer"), "true"));
      await session.step(51, "Then second grid viewer should have a \"header AGE\" area", () => hasArea(page, el("second grid viewer"), "header AGE"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second grid rebinds to another table and closes without taking the view's grid", async () => {
      await session.step(55, "When user opens demog-1000 dataset keeping the first 100 rows as \"demog-100\"", () => openDatasetRowsAs(page, ds("demog-1000"), 100, "demog-100"));
      await session.step(56, "Then the table should have 100 rows", () => rowCount(page, 100));
      await session.step(57, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(58, "And user sets \"Table\" property of second grid viewer to \"demog-100\"", () => setProperty(page, "Table", el("second grid viewer"), "demog-100"));
      await session.step(59, "Then second grid viewer should be bound to table \"demog-100\"", () => boundTable(page, el("second grid viewer"), "demog-100"));
      await session.step(60, "And second grid viewer should show 100 rows", () => showsRows(page, el("second grid viewer"), 100));
      await session.step(61, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(62, "When user clicks on close icon of second grid viewer", () => clickOn(page, el("close icon of second grid viewer")));
      await session.step(63, "Then the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
      await session.step(64, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(65, "And grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
      await session.step(66, "When user closes all views", () => closeAllViews(page));
      await session.step(67, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column tooltip menu marks the setting that is on, and Columns lists the chosen columns (GROK-20890)", async () => {
      await session.step(71, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(72, "When user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(73, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(74, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(75, "Then \"Default\" menu item in context menu should be visible", () => shouldBe(page, el("\"Default\" menu item in context menu"), "visible"));
      await session.step(76, "And \"Form\" menu item in context menu should be visible", () => shouldBe(page, el("\"Form\" menu item in context menu"), "visible"));
      await session.step(77, "And \"Columns\" menu item in context menu should be visible", () => shouldBe(page, el("\"Columns\" menu item in context menu"), "visible"));
      await session.step(78, "And \"None\" menu item in context menu should be visible", () => shouldBe(page, el("\"None\" menu item in context menu"), "visible"));
      await session.step(79, "And \"Default\" menu item in context menu should be selected", () => shouldBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(80, "And \"None\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(81, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(82, "And user picks \"Tooltip > Current Column > Columns\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Tooltip > Current Column > Columns", "header AGE", el("grid")));
      await session.step(83, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(84, "When user clicks on All link in \"Select columns...\" dialog", () => clickOn(page, el("All link in \"Select columns...\" dialog")));
      await session.step(85, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(86, "Then \"Select columns...\" dialog should be hidden", () => shouldBe(page, el("\"Select columns...\" dialog"), "hidden"));
      await session.step(87, "When user hovers over the \"cell 3 of AGE\" area of grid", () => hoverArea(page, "cell 3 of AGE", el("grid")));
      await session.step(88, "Then the tooltip should show columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => tooltipColumns(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(89, "And the tooltip should show \"AGE\" as \"58\"", () => tooltipValue(page, "AGE", "58"));
      await session.step(90, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(91, "And user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(92, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(93, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(94, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(95, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(96, "Then \"Columns\" menu item in context menu should be selected", () => shouldBe(page, el("\"Columns\" menu item in context menu"), "selected"));
      await session.step(97, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(98, "And user picks \"Tooltip > Current Column > None\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Tooltip > Current Column > None", "header AGE", el("grid")));
      await session.step(99, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(100, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(101, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(102, "Then \"None\" menu item in context menu should be selected", () => shouldBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(103, "And \"Default\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(104, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(105, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(106, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(107, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(108, "And user clicks on \"Default\" menu item in context menu", () => clickOn(page, el("\"Default\" menu item in context menu")));
      await session.step(109, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(110, "And user hovers over \"Tooltip\" menu item in context menu", () => hoverOver(page, el("\"Tooltip\" menu item in context menu")));
      await session.step(111, "And user hovers over \"Current Column\" menu item in context menu", () => hoverOver(page, el("\"Current Column\" menu item in context menu")));
      await session.step(112, "Then \"Default\" menu item in context menu should be selected", () => shouldBe(page, el("\"Default\" menu item in context menu"), "selected"));
      await session.step(113, "And \"None\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"None\" menu item in context menu"), "selected"));
      await session.step(114, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(115, "Then no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    await run.scenario("Pick Up and Apply carry the grid's look to a second grid", async () => {
      await session.step(118, "When user adds a grid viewer", () => addViewer(page, "grid"));
      await session.step(119, "Then the open tableview should have 2 grid viewers", () => viewerCount(page, 2, "grid"));
      await session.step(120, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 3 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 3 of AGE", el("grid")));
      await session.step(121, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(122, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(123, "And \"Color Coding\" property of second grid viewer should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("second grid viewer"), "Auto"));
      await session.step(124, "When user remembers the \"color of cell 2 of HEIGHT\" reading of grid", () => rememberReading(page, "color of cell 2 of HEIGHT", el("grid")));
      await session.step(125, "Then the \"color of cell 2 of HEIGHT\" reading of grid should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of cell 2 of HEIGHT", el("grid"), "#ffffff"));
      await session.step(126, "And the \"color of cell 2 of HEIGHT\" reading of second grid viewer should not be as remembered", () => readingNotAsRemembered(page, "color of cell 2 of HEIGHT", el("second grid viewer")));
      await session.step(127, "When user picks \"Pick Up / Apply > Pick Up\" from the context menu of the \"cell 3 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Pick Up", "cell 3 of AGE", el("grid")));
      await session.step(128, "And user picks \"Pick Up / Apply > Apply\" from the context menu of the \"cell 3 of AGE\" area of second grid viewer", () => pickFromAreaContextMenu(page, "Pick Up / Apply > Apply", "cell 3 of AGE", el("second grid viewer")));
      await session.step(129, "Then \"Color Coding\" property of second grid viewer should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("second grid viewer"), "All"));
      await session.step(130, "And the \"color of cell 2 of HEIGHT\" reading of second grid viewer should be as remembered", () => readingAsRemembered(page, "color of cell 2 of HEIGHT", el("second grid viewer")));
      await session.step(131, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 3 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 3 of AGE", el("grid")));
      await session.step(132, "And user clicks on close icon of second grid viewer", () => clickOn(page, el("close icon of second grid viewer")));
      await session.step(133, "Then the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
