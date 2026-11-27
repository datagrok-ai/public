/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/01-add-a-calculated-column.feature
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
import {clickOn, expand, holdsCode, hoverText, pressKey, pressKeyIn, replaceCode, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/01-add-a-calculated-column.feature", import.meta.url);
  test("Add a calculated column", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")), undefined, "Add New Column: a column computed by a formula");
    await session.step(21, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await session.step(22, "When user types \"BMI\" into column name input", () => typeInto(page, "BMI", el("column name input")));
    await session.step(24, "And user types \"pow\" into \"Search function\" input in \"Add New Column\" dialog", () => typeInto(page, "pow", el("\"Search function\" input in \"Add New Column\" dialog")), undefined, "The table's columns and the functions are listed on the right");
    await session.step(25, "And user hovers over the text \"Pow\" in \"Add New Column\" dialog", () => hoverText(page, "Pow", el("\"Add New Column\" dialog")));
    await session.step(27, "And user replaces the code of code editor in \"Add New Column\" dialog with \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"), undefined, "Columns go in as ${NAME}, next to functions and numbers");
    await session.step(28, "Then code editor in \"Add New Column\" dialog should hold the code \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\"", () => holdsCode(page, el("code editor in \"Add New Column\" dialog"), "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"));
    await session.step(29, "When the video says \"The preview shows the first values before you add the column\"", () => videoSays(page, "The preview shows the first values before you add the column"));
    await session.step(30, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
    await session.step(31, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
    await session.step(32, "And the table should have a column \"BMI\"", () => hasColumn(page, "BMI"));
    await session.step(34, "And the \"BMI\" cell of row 1 should be displayed as \"28.42\"", () => displayedInRow(page, "BMI", 1, "28.42"));
    await session.step(36, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The formula stays with the column: edit it in the context panel");
    await session.step(37, "And user clicks on the \"header BMI\" area of grid", () => clickArea(page, "header BMI", el("grid")));
    await session.step(38, "Then the context panel should show \"BMI\"", () => contextPanelShows(page, "BMI"));
    await session.step(39, "When user expands Formula pane in context panel", () => expand(page, el("Formula pane in context panel")));
    await session.step(41, "And user replaces the code of code editor in Formula pane with \"Round(${WEIGHT} / Pow(${HEIGHT} / 100, 2))\"", () => replaceCode(page, el("code editor in Formula pane"), "Round(${WEIGHT} / Pow(${HEIGHT} / 100, 2))"), undefined, "Round to whole numbers, then Apply");
    await session.step(42, "And user clicks on \"Apply\" button in Formula pane", () => clickOn(page, el("\"Apply\" button in Formula pane")));
    await session.step(43, "Then the \"BMI\" cell of row 1 should be displayed as \"28.00\"", () => displayedInRow(page, "BMI", 1, "28.00"));
    await session.step(45, "When user double-clicks on the \"cell 1 of HEIGHT\" area of grid", () => doubleClickArea(page, "cell 1 of HEIGHT", el("grid")), undefined, "Change a value it reads, and the column follows");
    await session.step(46, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
    await session.step(47, "And user types \"170\" into cell editor", () => typeInto(page, "170", el("cell editor")));
    await session.step(48, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(49, "Then the \"HEIGHT\" cell of row 1 should be displayed as \"170.000\"", () => displayedInRow(page, "HEIGHT", 1, "170.000"));
    await session.step(51, "When user hovers over the \"left edge of cell 1 of BMI\" area of grid", () => hoverArea(page, "left edge of cell 1 of BMI", el("grid")), undefined, "Row 1 now reads 25: 73.2 kg at 1.70 m");
    await session.step(53, "Then the \"BMI\" cell of row 1 should be displayed as \"25.00\"", () => displayedInRow(page, "BMI", 1, "25.00"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
