/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/persistence-sources.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.cp.add-new-column-persists, powerpack.int.add-new-column-datasync-roundtrip, GROK-17109]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, pressKeyIn, shouldBe, shouldBeSwitchedOn, shouldBecomeVisibleWithin, switchOn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, everyValuePlus, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, noProjectOnServer, projectsOnServer, refreshBrowse, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {copyInUsersFiles, noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {doubleClickArea, noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Calculated columns over a file in My files follow a rename and an edit, and survive a project round trip", () => {
  const session = feature(test, "features/add-new-column/persistence-sources.feature", import.meta.url);
  test("Calculated columns over a file in My files follow a rename and an edit, and survive a project round trip", {tag: ["@journey", "@realizes:powerpack.cp.add-new-column-persists", "@realizes:powerpack.int.add-new-column-datasync-roundtrip", "@realizes:GROK-17109"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And no project named \"bdd-anc-home-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-anc-home-{run}")));
    await session.step(21, "And a copy of the \"System:DemoFiles/demog.csv\" file is in the user's files as \"bdd-anc-{run}.csv\"", () => copyInUsersFiles(page, "System:DemoFiles/demog.csv", session.text("bdd-anc-{run}.csv")));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(24, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(25, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(26, "When user double-clicks on Files---My-files---bdd-anc-{run}.csv tree node inside browse tree", () => doubleClickOn(page, el(session.text("Files---My-files---bdd-anc-{run}.csv tree node inside browse tree"))));
    await session.step(27, "Then the \"bdd-anc-{run}\" view should be current", () => viewIsCurrent(page, session.text("bdd-anc-{run}")));
    await session.step(28, "And the table should have 5850 rows", () => rowCount(page, 5850));
    await run.scenario("Two chained columns, then a rename and an edit of their source", async () => {
      await session.step(31, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(32, "And user types \"Weight2\" into column name input", () => typeInto(page, "Weight2", el("column name input")));
      await session.step(33, "And user types \"${WEIGHT} + 100\" into formula editor", () => typeInto(page, "${WEIGHT} + 100", el("formula editor")));
      await session.step(34, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(35, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(36, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(37, "And user types \"Weight3\" into column name input", () => typeInto(page, "Weight3", el("column name input")));
      await session.step(38, "And user types \"${Weight2} + 100\" into formula editor", () => typeInto(page, "${Weight2} + 100", el("formula editor")));
      await session.step(39, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(40, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(41, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 100", () => everyValuePlus(page, "Weight2", "WEIGHT", 100));
      await session.step(42, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValuePlus(page, "Weight3", "Weight2", 100));
      await session.step(43, "When user picks \"Column Properties...\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header WEIGHT", el("grid")));
      await session.step(44, "And user types \"BaseWeight\" into \"New name:\" input in \"WEIGHT\" dialog", () => typeInto(page, "BaseWeight", el("\"New name:\" input in \"WEIGHT\" dialog")));
      await session.step(45, "And user clicks on OK button in \"WEIGHT\" dialog", () => clickOn(page, el("OK button in \"WEIGHT\" dialog")));
      await session.step(46, "Then the table should have a column \"BaseWeight\"", () => hasColumn(page, "BaseWeight"));
      await session.step(47, "And \"Weight2\" column should have tag \"formula\" equal to \"${BaseWeight} + 100\"", () => columnTag(page, "Weight2", "formula", "${BaseWeight} + 100"));
      await session.step(48, "When user double-clicks on the \"cell 1 of BaseWeight\" area of grid", () => doubleClickArea(page, "cell 1 of BaseWeight", el("grid")));
      await session.step(49, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(50, "And user types \"500\" into cell editor", () => typeInto(page, "500", el("cell editor")));
      await session.step(51, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(52, "Then the value of \"BaseWeight\" column in row 1 should be \"500\"", () => valueInRow(page, "BaseWeight", 1, "500"));
      await session.step(53, "And the value of \"Weight2\" column in row 1 should be \"600\"", () => valueInRow(page, "Weight2", 1, "600"));
      await session.step(54, "And every value of \"Weight2\" column should equal \"BaseWeight\" column plus 100", () => everyValuePlus(page, "Weight2", "BaseWeight", 100));
      await session.step(55, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValuePlus(page, "Weight3", "Weight2", 100));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Saved with Data sync, closed and reopened from Dashboards, the columns are there", async () => {
      await session.step(60, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
      await session.step(61, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(62, "When user enters \"bdd-anc-home-{run}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("bdd-anc-home-{run}"), el("Name text input in \"Save project\" dialog")));
      await session.step(63, "And user switches on Data sync input in \"Save project\" dialog", () => switchOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(64, "Then Data sync input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(65, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(66, "Then \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Save project\" dialog"), "hidden"));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(68, "And 1 project named \"bdd-anc-home-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-anc-home-{run}")));
      await session.step(69, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(70, "And user closes all views", () => closeAllViews(page));
      await session.step(71, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(72, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(73, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(74, "When user types \"bdd-anc-home-{run}\" into gallery search", () => typeInto(page, session.text("bdd-anc-home-{run}"), el("gallery search")));
      await session.step(75, "Then \"bdd-anc-home-{run}\" gallery card should become visible within 60 seconds", () => shouldBecomeVisibleWithin(page, el(session.text("\"bdd-anc-home-{run}\" gallery card")), 60));
      await session.step(76, "When user double-clicks on \"bdd-anc-home-{run}\" gallery card", () => doubleClickOn(page, el(session.text("\"bdd-anc-home-{run}\" gallery card"))));
      await session.step(77, "Then the \"bdd-anc-{run}\" view should be current", () => viewIsCurrent(page, session.text("bdd-anc-{run}")));
      await session.step(78, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(79, "And the table should have a column \"Weight2\"", () => hasColumn(page, "Weight2"));
      await session.step(80, "And the table should have a column \"Weight3\"", () => hasColumn(page, "Weight3"));
      await session.step(81, "And the table should have a column \"BaseWeight\"", () => hasColumn(page, "BaseWeight"));
      await session.step(82, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
