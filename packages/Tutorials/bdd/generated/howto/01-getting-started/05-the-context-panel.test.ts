/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/05-the-context-panel.feature
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
import {clickOn, doubleClickOn, expand} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShows, contextPanelShowsCurrentCell, openDataset, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/05-the-context-panel.feature", import.meta.url);
  test("The context panel", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "When user clicks on \"demog\" text in status bar", () => clickOn(page, el("\"demog\" text in status bar")), undefined, "The panel on the right shows what you click: here, the table's name");
    await session.step(24, "Then the context panel should show \"demog\"", () => contextPanelShows(page, "demog"));
    await session.step(26, "When user clicks on \"General\" accordion header in context panel", () => clickOn(page, el("\"General\" accordion header in context panel")), undefined, "Its panes open with a click");
    await session.step(28, "And user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")), undefined, "A column: click its header");
    await session.step(29, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(31, "When user clicks on \"Stats\" accordion header in context panel", () => clickOn(page, el("\"Stats\" accordion header in context panel")), undefined, "Stats: the column's numbers at a glance");
    await session.step(33, "And user clicks on \"Plots\" accordion header in context panel", () => clickOn(page, el("\"Plots\" accordion header in context panel")), undefined, "Plots: how the values are spread");
    await session.step(35, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Now a table of molecules, opened from Browse");
    await session.step(36, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(37, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(38, "And user expands Files---Demo---chem tree node inside browse tree", () => expand(page, el("Files---Demo---chem tree node inside browse tree")));
    await session.step(39, "And user double-clicks on Files---Demo---chem---smiles.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---chem---smiles.csv tree node inside browse tree")));
    await session.step(40, "Then the \"smiles\" view should be current", () => viewIsCurrent(page, "smiles"));
    await session.step(42, "When user clicks on the \"cell 3 of canonical_smiles\" area of grid", () => clickArea(page, "cell 3 of canonical_smiles", el("grid")), undefined, "A molecule: click its cell");
    await session.step(43, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(45, "When user clicks on \"Chemistry\" accordion header in context panel", () => clickOn(page, el("\"Chemistry\" accordion header in context panel")), undefined, "Chemistry > Properties: computed for this molecule");
    await session.step(46, "And user clicks on \"Properties\" accordion header in context panel", () => clickOn(page, el("\"Properties\" accordion header in context panel")));
    await session.step(47, "Then no errors should have been logged", () => noErrors(page));
  });
});
