/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/03-run-a-script-on-your-table.feature
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
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {scriptsSearchPutBack} from '../../../bindings/howto-scripting.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverText, selectIn, shouldBe, shouldContainText, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, contextPanelClosed, dialogCloses, openDataset, scriptOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/03-run-a-script-on-your-table.feature", import.meta.url);
  test("Run a script on your table", {tag: ["@guide", "@tutorial", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And a script \"Add BMI column\" is on the server:", () => scriptOnServer(page, "Add BMI column", "//language: javascript\n//description: Adds a BMI column from height in cm and weight in kg\n//input: dataframe table\n//input: column height {type: numerical}\n//input: column weight {type: numerical}\ntable.columns.addNewFloat('BMI').init((i) => weight.get(i) / (height.get(i) / 100) ** 2);"));
    await session.step(27, "And the Scripts gallery search is emptied at feature end", () => scriptsSearchPutBack(page));
    await session.step(28, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(29, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")), undefined, "A saved script will add a BMI column to demog");
    await session.step(34, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
    await session.step(35, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
    await session.step(36, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
    await session.step(37, "When user types \"BMI\" into gallery search", () => typeInto(page, "BMI", el("gallery search")));
    await session.step(38, "Then \"Add BMI column\" link in gallery should be visible", () => shouldBe(page, el("\"Add BMI column\" link in gallery"), "visible"));
    await session.step(40, "When user picks \"Edit...\" from the context menu of \"Add BMI column\" link in gallery", () => pickFromContextMenu(page, "Edit...", el("\"Add BMI column\" link in gallery")), undefined, "Edit... opens its code");
    await session.step(41, "Then the \"Add BMI column\" view should be current", () => viewIsCurrent(page, "Add BMI column"));
    await session.step(42, "And code editor should contain the text \"//input: dataframe table\"", () => shouldContainText(page, el("code editor"), "//input: dataframe table"));
    await session.step(43, "When the video says \"A dataframe input takes a table, a column input one of its columns\"", () => videoSays(page, "A dataframe input takes a table, a column input one of its columns"));
    await session.step(44, "And user hovers over the text \"//input: column height\" in code editor", () => hoverText(page, "//input: column height", el("code editor")));
    await session.step(46, "And user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")), undefined, "Run, then pick the table and the columns");
    await session.step(47, "Then \"Add BMI column\" dialog should be visible", () => shouldBe(page, el("\"Add BMI column\" dialog"), "visible"));
    await session.step(48, "When user selects \"demog\" in Table input in \"Add BMI column\" dialog", () => selectIn(page, "demog", el("Table input in \"Add BMI column\" dialog")));
    await session.step(49, "And user selects \"HEIGHT\" in Height input in \"Add BMI column\" dialog", () => selectIn(page, "HEIGHT", el("Height input in \"Add BMI column\" dialog")));
    await session.step(50, "And user selects \"WEIGHT\" in Weight input in \"Add BMI column\" dialog", () => selectIn(page, "WEIGHT", el("Weight input in \"Add BMI column\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Add BMI column\" dialog", () => clickOn(page, el("OK button in \"Add BMI column\" dialog")));
    await session.step(52, "Then the \"Add BMI column\" dialog should close", () => dialogCloses(page, "Add BMI column"));
    await session.step(53, "When user clicks on \"demog\" view tab", () => clickOn(page, el("\"demog\" view tab")));
    await session.step(54, "Then the table should have a column \"BMI\"", () => hasColumn(page, "BMI"));
    await session.step(57, "And the \"BMI\" cell of row 1 should be displayed as \"28.42\"", () => displayedInRow(page, "BMI", 1, "28.42"), undefined, "The script added a BMI column to the table");
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
  });
});
