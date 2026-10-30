/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/change-type.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {clickOn, enterInto, pressKey, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {allOfFiltered, filterPasses, filterPassesAll, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Change Type converts a column from its header menu or from Column Properties, and Ctrl+Z puts it back", () => {
  const session = feature(test, "features/table-view/change-type.feature", import.meta.url);
  test("Change Type converts a column from its header menu or from Column Properties, and Ctrl+Z puts it back", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens a table \"ct\" with:", () => openTableOf(page, "ct", [["name","code"],["a","10"],["b","20"],["c","n/a"],["d",""]]), [["name","code"],["a","10"],["b","20"],["c","n/a"],["d",""]]);
    await run.scenario("The header menu opens the dialog with the column's current type, and OK converts it to the new type", async () => {
      await session.step(26, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(27, "When user picks \"Change Type...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header code", el("grid")));
      await session.step(28, "Then \"Change column type\" dialog should be visible", () => shouldBe(page, el("\"Change column type\" dialog"), "visible"));
      await session.step(29, "And \"Type\" input in \"Change column type\" dialog should have value \"string\"", () => shouldHaveValue(page, el("\"Type\" input in \"Change column type\" dialog"), "string"));
      await session.step(30, "And \"Type\" input in \"Change column type\" dialog should be disabled", () => shouldBe(page, el("\"Type\" input in \"Change column type\" dialog"), "disabled"));
      await session.step(31, "When user selects \"int\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "int", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
      await session.step(33, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(34, "And \"code\" column should have type \"int\"", () => columnType(page, "code", "int"));
      await session.step(35, "And the value of \"code\" column in row 2 should be \"20\"", () => valueInRow(page, "code", 2, "20"));
      await session.step(36, "And the value of \"name\" column in row 2 should be \"b\"", () => valueInRow(page, "name", 2, "b"));
      await session.step(37, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(38, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(39, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(40, "And the value of \"code\" column in row 3 should be \"n/a\"", () => valueInRow(page, "code", 3, "n/a"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A value that cannot be converted takes the Not convertible value", async () => {
      await session.step(44, "When user picks \"Change Type...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header code", el("grid")));
      await session.step(45, "And user selects \"int\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "int", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(46, "And user enters \"-1\" into \"Not convertible\" input in \"Change column type\" dialog", () => enterInto(page, "-1", el("\"Not convertible\" input in \"Change column type\" dialog")));
      await session.step(47, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
      await session.step(48, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(49, "And \"code\" column should have type \"int\"", () => columnType(page, "code", "int"));
      await session.step(50, "And the value of \"code\" column in row 3 should be \"-1\"", () => valueInRow(page, "code", 3, "-1"));
      await session.step(51, "And the value of \"code\" column in row 1 should be \"10\"", () => valueInRow(page, "code", 1, "10"));
      await session.step(52, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(53, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(54, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(55, "And the value of \"code\" column in row 3 should be \"n/a\"", () => valueInRow(page, "code", 3, "n/a"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter icon beside Not convertible keeps only the rows that cannot be converted", async () => {
      await session.step(59, "When user picks \"Change Type...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header code", el("grid")));
      await session.step(60, "And user selects \"int\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "int", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(61, "And user clicks on second \"Filter\" icon in \"Change column type\" dialog", () => clickOn(page, el("second \"Filter\" icon in \"Change column type\" dialog")));
      await session.step(62, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(63, "And all rows where \"code\" is \"n/a\" should pass the filter", () => allOfFiltered(page, "code", "n/a"));
      await session.step(64, "When user clicks on CANCEL button in \"Change column type\" dialog", () => clickOn(page, el("CANCEL button in \"Change column type\" dialog")));
      await session.step(65, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(66, "When user resets the filter", () => resetFilter(page));
      await session.step(67, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Change type... in Column Properties opens the same dialog, and redo converts the column again", async () => {
      await session.step(71, "When user picks \"Column Properties...\" from the context menu of the \"header code\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header code", el("grid")));
      await session.step(72, "Then \"code\" dialog should be visible", () => shouldBe(page, el("\"code\" dialog"), "visible"));
      await session.step(73, "When user clicks on \"Change type...\" label in \"code\" dialog", () => clickOn(page, el("\"Change type...\" label in \"code\" dialog")));
      await session.step(74, "Then \"Change column type\" dialog should be visible", () => shouldBe(page, el("\"Change column type\" dialog"), "visible"));
      await session.step(75, "When user selects \"double\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "double", el("\"New type\" input in \"Change column type\" dialog")));
      await session.step(76, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
      await session.step(77, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
      await session.step(78, "And \"code\" column should have type \"double\"", () => columnType(page, "code", "double"));
      await session.step(79, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(80, "Then \"code\" dialog should be hidden", () => shouldBe(page, el("\"code\" dialog"), "hidden"));
      await session.step(81, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(82, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(83, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(84, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(85, "Then \"code\" column should have type \"double\"", () => columnType(page, "code", "double"));
      await session.step(86, "And the value of \"code\" column in row 1 should be \"10\"", () => valueInRow(page, "code", 1, "10"));
      await session.step(87, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(88, "Then \"code\" column should have type \"string\"", () => columnType(page, "code", "string"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
