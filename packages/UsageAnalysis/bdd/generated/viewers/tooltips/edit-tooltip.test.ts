/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/edit-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, shouldBe, shouldHaveValue, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverArea, noErrors, pickFromContextMenu, pointerAway, propertiesShouldBe, setProperties, tooltipColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, columnListExactly, columnListStartsWith, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing the table's tooltip from a viewer", () => {
  const session = feature(test, "features/viewers/tooltips/edit-tooltip.feature", import.meta.url);
  test("Editing the table's tooltip from a viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]), [["Show Tooltip","inherit from table"],["Show Column Names","Always"],["Show Visible Columns In Tooltip","true"]]);
    await session.step(32, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(33, "And user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(34, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await session.step(36, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["showLabels","Always"]]), [["showLabels","Always"]]);
    await run.scenario("Tooltip > Edit... opens the Edit Tooltip dialog with every column ticked", async () => {
      await session.step(40, "When user picks \"Tooltip > Edit...\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("scatter plot viewer")));
      await session.step(41, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
      await session.step(42, "And \"Use tooltip\" input in \"Edit Tooltip\" dialog should have value \"Table\"", () => shouldHaveValue(page, el("\"Use tooltip\" input in \"Edit Tooltip\" dialog"), "Table"));
      await session.step(43, "And \"Show column names\" input in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Show column names\" input in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(44, "And checkbox in \"Edit Tooltip\" dialog should be checked", () => shouldBe(page, el("checkbox in \"Edit Tooltip\" dialog"), "checked"));
      await session.step(45, "And the column list of \"Edit Tooltip\" dialog should be exactly \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(46, "And the \"AGE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(47, "And the \"SEVERITY\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "SEVERITY", el("\"Edit Tooltip\" dialog")));
      await session.step(48, "And \"Reset group tooltip\" label in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Reset group tooltip\" label in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(49, "And \"Design custom tooltip...\" label in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Design custom tooltip...\" label in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(50, "And OK button in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("OK button in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(51, "And CANCEL button in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(52, "And \"History\" icon in \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"History\" icon in \"Edit Tooltip\" dialog"), "visible"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search narrows the list whatever the case of the text", async () => {
      await session.step(56, "When user types \"race\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "race", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(57, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"RACE\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "RACE"));
      await session.step(58, "And the column list of \"Edit Tooltip\" dialog should start with \"RACE\"", () => columnListStartsWith(page, el("\"Edit Tooltip\" dialog"), "RACE"));
      await session.step(59, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(60, "And user types \"SEVer\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "SEVer", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(61, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "SEVERITY"));
      await session.step(62, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(63, "And user types \"ht\" into \"Search\" input in \"Edit Tooltip\" dialog", () => typeInto(page, "ht", el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(64, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"HEIGHT, WEIGHT\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "HEIGHT, WEIGHT"));
      await session.step(65, "When user clears \"Search\" input in \"Edit Tooltip\" dialog", () => clearField(page, el("\"Search\" input in \"Edit Tooltip\" dialog")));
      await session.step(66, "Then the column list of \"Edit Tooltip\" dialog should be exactly \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => columnListExactly(page, el("\"Edit Tooltip\" dialog"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The columns picked with OK become the tooltip of the plots and the grid", async () => {
      await session.step(70, "When user unchecks checkbox in \"Edit Tooltip\" dialog", () => uncheck(page, el("checkbox in \"Edit Tooltip\" dialog")));
      await session.step(71, "Then the \"AGE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(72, "And the \"RACE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
      await session.step(73, "When user toggles the \"AGE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
      await session.step(74, "And user toggles the \"SEX\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "SEX", el("\"Edit Tooltip\" dialog")));
      await session.step(75, "And user toggles the \"WEIGHT\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Edit Tooltip\" dialog")));
      await session.step(76, "And user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
      await session.step(77, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
      await session.step(78, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","HEIGHT"],["Y","WEIGHT"]]), [["X","HEIGHT"],["Y","WEIGHT"]]);
      await session.step(81, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(82, "Then the tooltip should show columns \"AGE, SEX, WEIGHT, HEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT, HEIGHT"));
      await session.step(83, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(84, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(85, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(86, "Then the tooltip should show columns \"AGE, SEX, WEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(87, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(88, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(89, "When user hovers over the \"cell 11 of AGE\" area of grid", () => hoverArea(page, "cell 11 of AGE", el("grid")));
      await session.step(90, "Then the tooltip should show columns \"AGE, SEX, WEIGHT\"", () => tooltipColumns(page, "AGE, SEX, WEIGHT"));
      await session.step(91, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(92, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
