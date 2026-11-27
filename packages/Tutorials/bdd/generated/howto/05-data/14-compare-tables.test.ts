/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/14-compare-tables.feature
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
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, deleteSelected, rowCount, selectWhereIs, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, clickViewTab, contextPanelClosed, openDatasetRowsAs, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/14-compare-tables.feature", import.meta.url);
  test("Compare two versions of a table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens cars dataset keeping the first 28 rows as \"cars v1\"", () => openDatasetRowsAs(page, ds("cars"), 28, "cars v1"));
    await session.step(17, "And user opens cars dataset keeping the first 30 rows as \"cars v2\"", () => openDatasetRowsAs(page, ds("cars"), 30, "cars v2"));
    await session.step(18, "And user sets \"price\" column in row 2 to \"18950\"", () => setCell(page, "price", 2, "18950"));
    await session.step(19, "And user sets \"price\" column in row 11 to \"7650\"", () => setCell(page, "price", 11, "7650"));
    await session.step(20, "And user selects rows where \"model\" is \"mercury\"", () => selectWhereIs(page, "model", "mercury"));
    await session.step(21, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(22, "And user clears the row selection", () => clearSelection(page));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(27, "When user clicks on the tab of the \"cars v1\" view", () => clickViewTab(page, "cars v1"), undefined, "Two versions of one table: cars v1, the older one");
    await session.step(29, "And user clicks on the tab of the \"cars v2\" view", () => clickViewTab(page, "cars v2"), undefined, "and cars v2, edited since");
    await session.step(31, "And user picks \"Data > Compare Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Tables..."), undefined, "Compare Tables: what changed from one to the other");
    await session.step(32, "Then Compare Tables dialog should be visible", () => shouldBe(page, el("Compare Tables dialog"), "visible"));
    await session.step(33, "And \"Tables\" input in Compare Tables dialog should have value \"cars v1\"", () => shouldHaveValue(page, el("\"Tables\" input in Compare Tables dialog"), "cars v1"));
    await session.step(35, "When user hovers over Key Columns input in Compare Tables dialog", () => hoverOver(page, el("Key Columns input in Compare Tables dialog")), undefined, "Key: the column that tells which row is which");
    await session.step(36, "Then Key Columns input in Compare Tables dialog should contain text \"model\"", () => shouldContainText(page, el("Key Columns input in Compare Tables dialog"), "model"));
    await session.step(37, "And Compare Tables dialog should not contain text \"is not unique\"", () => shouldNotContainText(page, el("Compare Tables dialog"), "is not unique"));
    await session.step(39, "When user hovers over \"Value Columns\" text in Compare Tables dialog", () => hoverOver(page, el("\"Value Columns\" text in Compare Tables dialog")), undefined, "Values: the columns to compare, all shared ones to start with");
    await session.step(40, "And user clicks on OK button in Compare Tables dialog", () => clickOn(page, el("OK button in Compare Tables dialog")));
    await session.step(41, "Then Compare Tables dialog should be hidden", () => shouldBe(page, el("Compare Tables dialog"), "hidden"));
    await session.step(42, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(43, "And the value of \"model\" column in row 3 should be \"mercury\"", () => valueInRow(page, "model", 3, "mercury"));
    await session.step(44, "And the value of \"model\" column in row 4 should be \"volvo2\"", () => valueInRow(page, "model", 4, "volvo2"));
    await session.step(47, "When the video says \"diff: only the rows that differ\" at the bottom", () => videoSaysAt(page, "diff: only the rows that differ", "bottom"));
    await session.step(49, "And user clicks on the \"cell 4 of model\" area of grid", () => clickArea(page, "cell 4 of model", el("grid")), undefined, "Green: rows new in cars v2");
    await session.step(50, "Then row 4 should be current", () => currentRowIs(page, 4));
    await session.step(52, "When user hovers over the \"cell 3 of model\" area of grid", () => hoverArea(page, "cell 3 of model", el("grid")), undefined, "Red: a row cars v2 no longer has");
    await session.step(54, "And user hovers over the \"cell 1 of price\" area of grid", () => hoverArea(page, "cell 1 of price", el("grid")), undefined, "Orange: a changed value; point at it for the old one");
    await session.step(55, "Then no errors should have been logged", () => noErrors(page));
  });
});
