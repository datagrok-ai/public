/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/column-properties.feature
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
import {clickOn, enterInto, hoverOver, pressKey, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, columnsExactly, displayedInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, pickFromAreaContextMenu, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Column Properties dialog edits a column's tags, and the grid shows what they mean", () => {
  const session = feature(test, "features/table-view/column-properties.feature", import.meta.url);
  test("The Column Properties dialog edits a column's tags, and the grid shows what they mean", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens a table \"props\" with:", () => openTableOf(page, "props", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]);
    await run.scenario("F2 opens the dialog of the current column with its name and the format the platform gave it", async () => {
      await session.step(25, "When user clicks on the \"cell 1 of weight\" area of grid", () => clickArea(page, "cell 1 of weight", el("grid")));
      await session.step(26, "And user presses F2", () => pressKey(page, "F2"));
      await session.step(27, "Then \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(28, "And \"New name:\" input in \"weight\" dialog should have value \"weight\"", () => shouldHaveValue(page, el("\"New name:\" input in \"weight\" dialog"), "weight"));
      await session.step(29, "And \"Key 1\" text input in \"weight\" dialog should have value \"format\"", () => shouldHaveValue(page, el("\"Key 1\" text input in \"weight\" dialog"), "format"));
      await session.step(30, "And \"Value 1\" text input in \"weight\" dialog should have value \"0.00\"", () => shouldHaveValue(page, el("\"Value 1\" text input in \"weight\" dialog"), "0.00"));
      await session.step(31, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(32, "Then \"weight\" dialog should be hidden", () => shouldBe(page, el("\"weight\" dialog"), "hidden"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A format tag changes how the grid shows the values, and Ctrl+Z and redo follow it", async () => {
      await session.step(36, "Then the \"weight\" cell of row 1 should be displayed as \"1.50\"", () => displayedInRow(page, "weight", 1, "1.50"));
      await session.step(37, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(38, "Then \"Key 1\" text input in \"weight\" dialog should have value \"format\"", () => shouldHaveValue(page, el("\"Key 1\" text input in \"weight\" dialog"), "format"));
      await session.step(39, "When user enters \"0.000\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "0.000", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(40, "And user clicks on OK button in \"weight\" dialog", () => clickOn(page, el("OK button in \"weight\" dialog")));
      await session.step(41, "Then the \"weight\" dialog should close", () => dialogCloses(page, "weight"));
      await session.step(42, "And the \"weight\" cell of row 1 should be displayed as \"1.500\"", () => displayedInRow(page, "weight", 1, "1.500"));
      await session.step(43, "And the \"weight\" cell of row 2 should be displayed as \"2.250\"", () => displayedInRow(page, "weight", 2, "2.250"));
      await session.step(44, "And the \"name\" cell of row 1 should be displayed as \"Ann\"", () => displayedInRow(page, "name", 1, "Ann"));
      await session.step(45, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(46, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(47, "Then the \"weight\" cell of row 1 should be displayed as \"1.50\"", () => displayedInRow(page, "weight", 1, "1.50"));
      await session.step(48, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(49, "Then the \"weight\" cell of row 1 should be displayed as \"1.500\"", () => displayedInRow(page, "weight", 1, "1.500"));
      await session.step(50, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(51, "Then the \"weight\" cell of row 1 should be displayed as \"1.50\"", () => displayedInRow(page, "weight", 1, "1.50"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Description and units added as two tag rows show in the header tooltip", async () => {
      await session.step(55, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(56, "And user clicks on first \"Add row\" button in \"weight\" dialog", () => clickOn(page, el("first \"Add row\" button in \"weight\" dialog")));
      await session.step(57, "And user enters \"description\" into \"Key 1\" text input in \"weight\" dialog", () => enterInto(page, "description", el("\"Key 1\" text input in \"weight\" dialog")));
      await session.step(58, "And user enters \"Body weight\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "Body weight", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(59, "And user clicks on first \"Add row\" button in \"weight\" dialog", () => clickOn(page, el("first \"Add row\" button in \"weight\" dialog")));
      await session.step(60, "Then \"Key 2\" text input in \"weight\" dialog should have value \"description\"", () => shouldHaveValue(page, el("\"Key 2\" text input in \"weight\" dialog"), "description"));
      await session.step(61, "When user enters \"units\" into \"Key 1\" text input in \"weight\" dialog", () => enterInto(page, "units", el("\"Key 1\" text input in \"weight\" dialog")));
      await session.step(62, "And user enters \"kg\" into \"Value 1\" text input in \"weight\" dialog", () => enterInto(page, "kg", el("\"Value 1\" text input in \"weight\" dialog")));
      await session.step(63, "And user clicks on OK button in \"weight\" dialog", () => clickOn(page, el("OK button in \"weight\" dialog")));
      await session.step(64, "Then the \"weight\" dialog should close", () => dialogCloses(page, "weight"));
      await session.step(65, "And \"weight\" column should have units \"kg\"", () => columnUnits(page, "weight", "kg"));
      await session.step(66, "When user hovers over the \"header weight\" area of grid", () => hoverArea(page, "header weight", el("grid")));
      await session.step(67, "Then tooltip should contain text \"Body weight\"", () => shouldContainText(page, el("tooltip"), "Body weight"));
      await session.step(68, "And tooltip should contain text \"units: kg\"", () => shouldContainText(page, el("tooltip"), "units: kg"));
      await session.step(69, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(70, "And user hovers over the \"header city\" area of grid", () => hoverArea(page, "header city", el("grid")));
      await session.step(71, "Then tooltip should not contain text \"Body weight\"", () => shouldNotContainText(page, el("tooltip"), "Body weight"));
      await session.step(72, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(73, "And user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(74, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(75, "Then \"weight\" column should have units \"\"", () => columnUnits(page, "weight", ""));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog shows a tag set earlier, and removing its row removes the tag", async () => {
      await session.step(79, "When user picks \"Column Properties...\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header city", el("grid")));
      await session.step(80, "And user clicks on first \"Add row\" button in \"city\" dialog", () => clickOn(page, el("first \"Add row\" button in \"city\" dialog")));
      await session.step(81, "And user enters \"units\" into \"Key 1\" text input in \"city\" dialog", () => enterInto(page, "units", el("\"Key 1\" text input in \"city\" dialog")));
      await session.step(82, "And user enters \"town\" into \"Value 1\" text input in \"city\" dialog", () => enterInto(page, "town", el("\"Value 1\" text input in \"city\" dialog")));
      await session.step(83, "And user clicks on OK button in \"city\" dialog", () => clickOn(page, el("OK button in \"city\" dialog")));
      await session.step(84, "Then \"city\" column should have units \"town\"", () => columnUnits(page, "city", "town"));
      await session.step(85, "When user picks \"Column Properties...\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header city", el("grid")));
      await session.step(86, "Then \"Key 1\" text input in \"city\" dialog should have value \"units\"", () => shouldHaveValue(page, el("\"Key 1\" text input in \"city\" dialog"), "units"));
      await session.step(87, "And \"Value 1\" text input in \"city\" dialog should have value \"town\"", () => shouldHaveValue(page, el("\"Value 1\" text input in \"city\" dialog"), "town"));
      await session.step(88, "When user clicks on first \"Remove row\" button in \"city\" dialog", () => clickOn(page, el("first \"Remove row\" button in \"city\" dialog")));
      await session.step(89, "And user clicks on OK button in \"city\" dialog", () => clickOn(page, el("OK button in \"city\" dialog")));
      await session.step(90, "Then the \"city\" dialog should close", () => dialogCloses(page, "city"));
      await session.step(91, "And \"city\" column should have units \"\"", () => columnUnits(page, "city", ""));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A name another column already has is refused", async () => {
      await session.step(95, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(96, "And user enters \"City\" into \"New name:\" input in \"weight\" dialog", () => enterInto(page, "City", el("\"New name:\" input in \"weight\" dialog")));
      await session.step(97, "Then \"New name:\" input in \"weight\" dialog should be invalid", () => shouldBe(page, el("\"New name:\" input in \"weight\" dialog"), "invalid"));
      await session.step(98, "When user hovers over \"New name:\" input in \"weight\" dialog", () => hoverOver(page, el("\"New name:\" input in \"weight\" dialog")));
      await session.step(99, "Then tooltip should contain text \"Column 'City' already exists.\"", () => shouldContainText(page, el("tooltip"), "Column 'City' already exists."));
      await session.step(100, "When user enters \"mass\" into \"New name:\" input in \"weight\" dialog", () => enterInto(page, "mass", el("\"New name:\" input in \"weight\" dialog")));
      await session.step(101, "Then \"New name:\" input in \"weight\" dialog should be valid", () => shouldBe(page, el("\"New name:\" input in \"weight\" dialog"), "valid"));
      await session.step(102, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(103, "Then \"weight\" dialog should be hidden", () => shouldBe(page, el("\"weight\" dialog"), "hidden"));
      await session.step(104, "And the table should have the columns \"name, city, weight\"", () => columnsExactly(page, "name, city, weight"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit category colors... is offered for a text column and not for a number column", async () => {
      await session.step(108, "When user picks \"Column Properties...\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header weight", el("grid")));
      await session.step(109, "Then \"Edit category colors...\" label in \"weight\" dialog should be disabled", () => shouldBe(page, el("\"Edit category colors...\" label in \"weight\" dialog"), "disabled"));
      await session.step(110, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(111, "And user picks \"Column Properties...\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header city", el("grid")));
      await session.step(112, "Then \"Edit category colors...\" label in \"city\" dialog should be enabled", () => shouldBe(page, el("\"Edit category colors...\" label in \"city\" dialog"), "enabled"));
      await session.step(113, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(114, "Then \"city\" dialog should be hidden", () => shouldBe(page, el("\"city\" dialog"), "hidden"));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
