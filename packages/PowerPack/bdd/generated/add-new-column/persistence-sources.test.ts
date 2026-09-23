/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/persistence-sources.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.cp.add-new-column-persists, powerpack.int.add-new-column-datasync-roundtrip, GROK-17109]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {everyValueEquals, fileInHome, noTablesOpen} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, pressKeyIn, shouldBe, shouldBeSwitchedOn, shouldBecomeVisibleWithin, switchOn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, noProjectOnServer, projectsOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Calculated columns over a file in My files follow a rename and an edit, and survive a project round trip", () => {
  const session = feature(test, "features/add-new-column/persistence-sources.feature", import.meta.url);
  test("Calculated columns over a file in My files follow a rename and an edit, and survive a project round trip", {tag: ["@journey", "@realizes:powerpack.cp.add-new-column-persists", "@realizes:powerpack.int.add-new-column-datasync-roundtrip", "@realizes:GROK-17109"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And no project named \"bdd-anc-home-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-anc-home-{run}")));
    await session.step(18, "And a copy of the \"System:DemoFiles/demog.csv\" file is in the home folder as \"bdd-anc-{run}.csv\"", () => fileInHome(page, "System:DemoFiles/demog.csv", session.text("bdd-anc-{run}.csv")));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(22, "When user double-clicks on Files---My-files---bdd-anc-{run}.csv tree node inside browse tree", () => doubleClickOn(page, el(session.text("Files---My-files---bdd-anc-{run}.csv tree node inside browse tree"))));
    await session.step(23, "Then the \"bdd-anc-{run}\" view should be current", () => viewIsCurrent(page, session.text("bdd-anc-{run}")));
    await session.step(24, "And the table should have 5850 rows", () => rowCount(page, 5850));
    await run.scenario("Two chained columns, then a rename and an edit of their source", async () => {
      await session.step(27, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(28, "And user types \"Weight2\" into column name input", () => typeInto(page, "Weight2", el("column name input")));
      await session.step(29, "And user types \"${WEIGHT} + 100\" into formula editor", () => typeInto(page, "${WEIGHT} + 100", el("formula editor")));
      await session.step(30, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(31, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(32, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(33, "And user types \"Weight3\" into column name input", () => typeInto(page, "Weight3", el("column name input")));
      await session.step(34, "And user types \"${Weight2} + 100\" into formula editor", () => typeInto(page, "${Weight2} + 100", el("formula editor")));
      await session.step(35, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(36, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(37, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 100", () => everyValueEquals(page, "Weight2", "WEIGHT", 100));
      await session.step(38, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValueEquals(page, "Weight3", "Weight2", 100));
      await session.step(39, "When user picks \"Column Properties...\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header WEIGHT", el("grid")));
      await session.step(40, "And user types \"BaseWeight\" into \"New name:\" input in \"WEIGHT\" dialog", () => typeInto(page, "BaseWeight", el("\"New name:\" input in \"WEIGHT\" dialog")));
      await session.step(41, "And user clicks on OK button in \"WEIGHT\" dialog", () => clickOn(page, el("OK button in \"WEIGHT\" dialog")));
      await session.step(42, "Then the table should have a column \"BaseWeight\"", () => hasColumn(page, "BaseWeight"));
      await session.step(43, "And \"Weight2\" column should have tag \"formula\" equal to \"${BaseWeight} + 100\"", () => columnTag(page, "Weight2", "formula", "${BaseWeight} + 100"));
      await session.step(44, "When user double-clicks on the \"cell 1 of BaseWeight\" area of grid", () => doubleClickArea(page, "cell 1 of BaseWeight", el("grid")));
      await session.step(45, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(46, "And user types \"500\" into cell editor", () => typeInto(page, "500", el("cell editor")));
      await session.step(47, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(48, "Then the value of \"BaseWeight\" column in row 1 should be \"500\"", () => valueInRow(page, "BaseWeight", 1, "500"));
      await session.step(49, "And the value of \"Weight2\" column in row 1 should be \"600\"", () => valueInRow(page, "Weight2", 1, "600"));
      await session.step(50, "And every value of \"Weight2\" column should equal \"BaseWeight\" column plus 100", () => everyValueEquals(page, "Weight2", "BaseWeight", 100));
      await session.step(51, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValueEquals(page, "Weight3", "Weight2", 100));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Saved with Data sync, closed and reopened, the columns keep their formulas", async () => {
      await session.step(56, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
      await session.step(57, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(58, "When user enters \"bdd-anc-home-{run}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("bdd-anc-home-{run}"), el("Name text input in \"Save project\" dialog")));
      await session.step(59, "And user switches on Data sync input in \"Save project\" dialog", () => switchOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(60, "Then Data sync input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(61, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(62, "Then \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Save project\" dialog"), "hidden"));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(64, "And 1 project named \"bdd-anc-home-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-anc-home-{run}")));
      await session.step(65, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(66, "And user closes all views", () => closeAllViews(page));
      await session.step(67, "Then no table should be open", () => noTablesOpen(page));
      await session.step(68, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(69, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(70, "When user types \"bdd-anc-home-{run}\" into gallery search", () => typeInto(page, session.text("bdd-anc-home-{run}"), el("gallery search")));
      await session.step(71, "Then \"bdd-anc-home-{run}\" project card should become visible within 60 seconds", () => shouldBecomeVisibleWithin(page, el(session.text("\"bdd-anc-home-{run}\" project card")), 60));
      await session.step(72, "When user double-clicks on \"bdd-anc-home-{run}\" project card", () => doubleClickOn(page, el(session.text("\"bdd-anc-home-{run}\" project card"))));
      await session.step(73, "Then the \"bdd-anc-{run}\" view should be current", () => viewIsCurrent(page, session.text("bdd-anc-{run}")));
      await session.step(74, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(75, "And the table should have a column \"Weight2\"", () => hasColumn(page, "Weight2"));
      await session.step(76, "And the table should have a column \"Weight3\"", () => hasColumn(page, "Weight3"));
      await session.step(77, "And the table should have a column \"BaseWeight\"", () => hasColumn(page, "BaseWeight"));
      await session.step(78, "And \"Weight2\" column should have tag \"formula\" equal to \"${BaseWeight} + 100\"", () => columnTag(page, "Weight2", "formula", "${BaseWeight} + 100"));
      await session.step(79, "And \"Weight3\" column should have tag \"formula\" equal to \"${Weight2} + 100\"", () => columnTag(page, "Weight3", "formula", "${Weight2} + 100"));
      await session.step(81, "And the value of \"BaseWeight\" column in row 1 should be \"73.19999694824219\"", () => valueInRow(page, "BaseWeight", 1, "73.19999694824219"));
      await session.step(82, "And every value of \"Weight2\" column should equal \"BaseWeight\" column plus 100", () => everyValueEquals(page, "Weight2", "BaseWeight", 100));
      await session.step(83, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValueEquals(page, "Weight3", "Weight2", 100));
      await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On the reopened table a second rename and edit are followed too", async () => {
      await session.step(88, "When user picks \"Column Properties...\" from the context menu of the \"header BaseWeight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header BaseWeight", el("grid")));
      await session.step(89, "And user types \"BaseWeight2\" into \"New name:\" input in \"BaseWeight\" dialog", () => typeInto(page, "BaseWeight2", el("\"New name:\" input in \"BaseWeight\" dialog")));
      await session.step(90, "And user clicks on OK button in \"BaseWeight\" dialog", () => clickOn(page, el("OK button in \"BaseWeight\" dialog")));
      await session.step(91, "Then the table should have a column \"BaseWeight2\"", () => hasColumn(page, "BaseWeight2"));
      await session.step(92, "And \"Weight2\" column should have tag \"formula\" equal to \"${BaseWeight2} + 100\"", () => columnTag(page, "Weight2", "formula", "${BaseWeight2} + 100"));
      await session.step(93, "And \"Weight3\" column should have tag \"formula\" equal to \"${Weight2} + 100\"", () => columnTag(page, "Weight3", "formula", "${Weight2} + 100"));
      await session.step(94, "When user double-clicks on the \"cell 2 of BaseWeight2\" area of grid", () => doubleClickArea(page, "cell 2 of BaseWeight2", el("grid")));
      await session.step(95, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(96, "And user types \"400\" into cell editor", () => typeInto(page, "400", el("cell editor")));
      await session.step(97, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(98, "Then the value of \"BaseWeight2\" column in row 2 should be \"400\"", () => valueInRow(page, "BaseWeight2", 2, "400"));
      await session.step(99, "And the value of \"Weight3\" column in row 2 should be \"600\"", () => valueInRow(page, "Weight3", 2, "600"));
      await session.step(100, "And every value of \"Weight2\" column should equal \"BaseWeight2\" column plus 100", () => everyValueEquals(page, "Weight2", "BaseWeight2", 100));
      await session.step(101, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 100", () => everyValueEquals(page, "Weight3", "Weight2", 100));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
