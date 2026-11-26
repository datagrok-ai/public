/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/13-scripting/06-give-a-script-a-friendly-form.feature
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
import {clickOn, dragSliderTo, expand, hoverText, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelClosed, dialogCloses, openDataset, scriptOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Scripting", () => {
  const session = feature(test, "features/howto/13-scripting/06-give-a-script-a-friendly-form.feature", import.meta.url);
  test("Give a script a friendly form", {tag: ["@guide", "@tutorial", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And a script \"Flag outliers\" is on the server:", () => scriptOnServer(page, "Flag outliers", "//language: javascript\n//description: Selects or filters out the rows where a value is far from the column's average\n//input: dataframe table\n//input: column column {type: numerical}\n//input: double deviations = 3 {caption: Standard deviations; min: 1; max: 4; step: 0.5}\n//input: string action = \"Select\" {caption: Then; choices: [\"Select\", \"Filter out\"]}\nconst mean = column.stats.avg;\nconst sd = column.stats.stdev;\nconst far = (i) => !column.isNone(i) && Math.abs(column.get(i) - mean) > deviations * sd;\nif (action === 'Select')\n  table.selection.init(far);\nelse\n  table.filter.init((i) => !far(i));"));
    await session.step(35, "And the Scripts gallery search is emptied at feature end", () => scriptsSearchPutBack(page));
    await session.step(36, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(37, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(38, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(39, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(41, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")), undefined, "A saved script that finds values far from the average");
    await session.step(42, "And user expands \"Platform > Functions\" tree node inside browse tree", () => expand(page, el("\"Platform > Functions\" tree node inside browse tree")));
    await session.step(43, "And user clicks on \"Platform > Functions > Scripts\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Functions > Scripts\" tree node inside browse tree")));
    await session.step(44, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
    await session.step(45, "When user types \"outliers\" into gallery search", () => typeInto(page, "outliers", el("gallery search")));
    await session.step(46, "And user picks \"Edit...\" from the context menu of \"Flag outliers\" link in gallery", () => pickFromContextMenu(page, "Edit...", el("\"Flag outliers\" link in gallery")));
    await session.step(47, "Then the \"Flag outliers\" view should be current", () => viewIsCurrent(page, "Flag outliers"));
    await session.step(49, "When user hovers over the text \"//input: column column\" in code editor", () => hoverText(page, "//input: column column", el("code editor")), undefined, "{type: numerical} offers only number columns");
    await session.step(51, "And user hovers over the text \"//input: double deviations\" in code editor", () => hoverText(page, "//input: double deviations", el("code editor")), undefined, "min and max add a slider");
    await session.step(53, "And user hovers over the text \"//input: string action\" in code editor", () => hoverText(page, "//input: string action", el("code editor")), undefined, "choices make a drop-down list");
    await session.step(55, "And user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")), undefined, "Run: the form follows the header");
    await session.step(56, "Then \"Flag outliers\" dialog should be visible", () => shouldBe(page, el("\"Flag outliers\" dialog"), "visible"));
    await session.step(57, "When user selects \"demog\" in Table input in \"Flag outliers\" dialog", () => selectIn(page, "demog", el("Table input in \"Flag outliers\" dialog")));
    await session.step(58, "And user selects \"WEIGHT\" in Column input in \"Flag outliers\" dialog", () => selectIn(page, "WEIGHT", el("Column input in \"Flag outliers\" dialog")));
    await session.step(59, "And user drags the slider of \"Standard deviations\" input in \"Flag outliers\" dialog to 2.5", () => dragSliderTo(page, el("\"Standard deviations\" input in \"Flag outliers\" dialog"), 2.5));
    await session.step(60, "And user selects \"Select\" in Then input in \"Flag outliers\" dialog", () => selectIn(page, "Select", el("Then input in \"Flag outliers\" dialog")));
    await session.step(61, "And user clicks on OK button in \"Flag outliers\" dialog", () => clickOn(page, el("OK button in \"Flag outliers\" dialog")));
    await session.step(62, "Then the \"Flag outliers\" dialog should close", () => dialogCloses(page, "Flag outliers"));
    await session.step(63, "When user clicks on \"demog\" view tab", () => clickOn(page, el("\"demog\" view tab")));
    await session.step(64, "Then 114 rows should be selected", () => selectedRowCount(page, 114));
    await session.step(66, "When user hovers over the text \"Selected: 114\" in status bar", () => hoverText(page, "Selected: 114", el("status bar")), undefined, "114 weights are over 2.5 standard deviations from the average");
    await session.step(67, "Then no errors should have been logged", () => noErrors(page));
  });
});
