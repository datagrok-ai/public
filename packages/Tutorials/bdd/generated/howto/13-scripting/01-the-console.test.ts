/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/01-the-console.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, pasteInto, pressKey, replaceCode, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, consoleShows, contextPanelClosed, noteConsole, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/01-the-console.feature", import.meta.url);
  test("The console", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "When user clicks on \"Console\" status bar toggle", () => clickOn(page, el("\"Console\" status bar toggle")), undefined, "The console: commands and their results");
    await session.step(22, "Then console should be visible", () => shouldBe(page, el("console"), "visible"));
    await session.step(23, "When user notes the console output", () => noteConsole(page));
    await session.step(25, "And user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")), undefined, "Add a column with the dialog, as usual");
    await session.step(26, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await session.step(27, "When user types \"BMI\" into column name input", () => typeInto(page, "BMI", el("column name input")));
    await session.step(28, "And user replaces the code of code editor in \"Add New Column\" dialog with \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"));
    await session.step(29, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
    await session.step(30, "Then the table should have a column \"BMI\"", () => hasColumn(page, "BMI"));
    await session.step(32, "And the console should show \"AddNewColumn(\\\"demog\\\", \\\"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\\\", \\\"BMI\\\")\"", () => consoleShows(page, "AddNewColumn(\"demog\", \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\", \"BMI\")"), undefined, "The console logged it as a command");
    await session.step(34, "When user clicks on console command line", () => clickOn(page, el("console command line")), undefined, "Type it here with another formula and name, then Enter");
    await session.step(35, "And user pastes \"AddNewColumn(\\\"demog\\\", \\\"${WEIGHT} * 2.2\\\", \\\"Weight in lb\\\")\" into console command line", () => pasteInto(page, "AddNewColumn(\"demog\", \"${WEIGHT} * 2.2\", \"Weight in lb\")", el("console command line")));
    await session.step(36, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(37, "Then the table should have a column \"Weight in lb\"", () => hasColumn(page, "Weight in lb"));
    await session.step(39, "And the \"Weight in lb\" cell of row 1 should be displayed as \"161.04\"", () => displayedInRow(page, "Weight in lb", 1, "161.04"));
    await session.step(41, "When user types \"SelectRandomRows\" into console command line", () => typeInto(page, "SelectRandomRows", el("console command line")), undefined, "A command's name alone opens a form for its inputs");
    await session.step(42, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(43, "Then \"selectRandomRows\" dialog should be visible", () => shouldBe(page, el("\"selectRandomRows\" dialog"), "visible"));
    await session.step(44, "When user enters \"100\" into Count input in \"selectRandomRows\" dialog", () => enterInto(page, "100", el("Count input in \"selectRandomRows\" dialog")));
    await session.step(45, "And user clicks on OK button in \"selectRandomRows\" dialog", () => clickOn(page, el("OK button in \"selectRandomRows\" dialog")));
    await session.step(46, "Then 100 rows should be selected", () => selectedRowCount(page, 100));
    await session.step(48, "And the console should show \"selectRandomRows(\\\"demog\\\", 100, 42, true)\"", () => consoleShows(page, "selectRandomRows(\"demog\", 100, 42, true)"), undefined, "The console shows the command with the values you chose");
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
