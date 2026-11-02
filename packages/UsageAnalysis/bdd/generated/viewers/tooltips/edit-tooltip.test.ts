/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/edit-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips, public#1208, public#373]
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
import {clearField, clickOn, shouldBe, shouldHaveValue, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {closeCurrentView, dialogCloses, openDataset, switchTableView, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, hoverArea, noErrors, pickFromContextMenu, pointerAway, propertiesShouldBe, setProperties, setProperty, tooltipColumnsInOrder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, columnListExactly, columnListLacks, columnListStartsWith, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing the table's tooltip from a viewer", () => {
  const session = feature(test, "features/viewers/tooltips/edit-tooltip.feature", import.meta.url);
  test("Editing the table's tooltip from a viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips", "@realizes:public#1208", "@realizes:public#373"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await session.step(34, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(35, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(36, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(38, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Tooltip > Edit... opens the Edit Tooltip dialog with every column ticked", async () => {
      await session.step(42, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(43, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(44, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Table\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Table"));
      await session.step(45, "And \"Show column names\" input in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Show column names\" input in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(46, "And checkbox in \"Edit Tooltip\" dialog should be checked", () => shouldBe(page, el("checkbox in \"Edit Tooltip\" dialog"), "checked"));
      await session.step(47, "And the column list of \"Edit Tooltip\" dialog should be exactly \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(48, "And the \"AGE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(49, "And the \"SEVERITY\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "SEVERITY", el("\"Edit Tooltip\" dialog")));
      await session.step(50, "And \"Reset group tooltip\" label in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Reset group tooltip\" label in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(51, "And \"Design custom tooltip...\" label in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Design custom tooltip...\" label in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(52, "And OK button in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("OK button in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(53, "And CANCEL button in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(54, "And \"History\" icon in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"History\" icon in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search narrows the list whatever the case of the text", async () => {
      await session.step(58, "When user types \"race\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "race", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(59, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"RACE\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "RACE"));
      await session.step(60, "And the column list of \"Edit Tooltip\" dialog should start with \"RACE\"", () => columnListStartsWith(page, el("\"Edit Tooltip\" dialog"), "RACE"));
      await session.step(61, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(62, "And user types \"SEVer\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "SEVer", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(63, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "SEVERITY"));
      await session.step(64, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(65, "And user types \"ht\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "ht", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(66, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"HEIGHT, WEIGHT\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "HEIGHT, WEIGHT"));
      await session.step(67, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(68, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The columns picked with OK become the tooltip of the plots and the grid, in a cloned view too (public#1208)", async () => {
      await session.step(73, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(74, "Then the \"AGE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(75, "And the \"RACE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
      await session.step(76, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(77, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
      await session.step(78, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
      await session.step(79, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(80, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(81, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","HEIGHT"],["Y","WEIGHT"]]), [["X","HEIGHT"],["Y","WEIGHT"]]);
      await session.step(84, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(85, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
      await session.step(86, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(87, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(88, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(89, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
      await session.step(90, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(91, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(92, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(93, "Then the tooltip should show columns in the order \"AGE, SEX, WEIGHT\"", () => tooltipColumnsInOrder(page, "AGE, SEX, WEIGHT"));
      await session.step(94, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(95, "And user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(96, "Then the \"demog-1000 copy\" view should be current", () => viewIsCurrent(page, "demog-1000 copy"));
      await session.step(97, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(98, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
      await session.step(99, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(100, "And user closes the current view", () => closeCurrentView(page));
      await session.step(101, "And user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(102, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(103, "Then the tooltip should show columns in the order \"HEIGHT, WEIGHT, AGE, SEX\"", () => tooltipColumnsInOrder(page, "HEIGHT, WEIGHT, AGE, SEX"));
      await session.step(104, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(105, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit Tooltip on a viewer bound to another table lists that table's columns (public#373)", async () => {
      await session.step(109, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(110, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(111, "When user sets \"table\" property of scatter plot viewer to \"spgi-100\"", () => setProperty(page, "table", el("scatter plot viewer"), "spgi-100"));
      await session.step(112, "Then scatter plot viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("scatter plot viewer"), "spgi-100"));
      await session.step(113, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(114, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(115, "And the column list of \"Edit Tooltip\" dialog should start with \"Id\"", () => columnListStartsWith(page, el("\"Edit Tooltip\" dialog"), "Id"));
      await session.step(116, "And the column list of \"Edit Tooltip\" dialog should not list \"USUBJID\"", () => columnListLacks(page, el("\"Edit Tooltip\" dialog"), "USUBJID"));
      await session.step(117, "And the column list of \"Edit Tooltip\" dialog should not list \"AGE\"", () => columnListLacks(page, el("\"Edit Tooltip\" dialog"), "AGE"));
      await session.step(118, "When user clicks on CANCEL button in \"Edit Tooltip\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Tooltip\" dialog")));
      await session.step(119, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
