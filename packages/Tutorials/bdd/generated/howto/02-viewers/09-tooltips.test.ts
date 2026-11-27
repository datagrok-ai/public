/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/09-tooltips.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, pickFromAreaContextMenu, tooltipColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/09-tooltips.feature", import.meta.url);
  test("Choose what a tooltip shows", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","WEIGHT"]]), [["X","HEIGHT"],["Y","WEIGHT"]]);
    await session.step(24, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(25, "When the video says \"Point at a marker to see its row\"", () => videoSays(page, "Point at a marker to see its row"));
    await session.step(27, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")), undefined, "The tooltip lists the row's values");
    await session.step(28, "Then the tooltip should show columns \"HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP\"", () => tooltipColumns(page, "HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP"));
    await session.step(30, "When user picks \"Tooltip > Edit...\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tooltip > Edit...", "view", el("scatter plot viewer")), undefined, "Right-click > Tooltip > Edit... picks the columns");
    await session.step(31, "Then \"Edit Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Tooltip\" dialog"), "visible"));
    await session.step(33, "When user toggles the \"USUBJID\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "USUBJID", el("\"Edit Tooltip\" dialog")), undefined, "Untick the columns you don't need");
    await session.step(34, "And user toggles the \"RACE\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(35, "And user toggles the \"DEMOG\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "DEMOG", el("\"Edit Tooltip\" dialog")));
    await session.step(36, "And user toggles the \"CONTROL\" column in the column list of \"Edit Tooltip\" dialog", () => toggleInColumnList(page, "CONTROL", el("\"Edit Tooltip\" dialog")));
    await session.step(37, "Then the \"AGE\" column should be checked in the column list of \"Edit Tooltip\" dialog", () => checkedInColumnList(page, "AGE", el("\"Edit Tooltip\" dialog")));
    await session.step(38, "And the \"RACE\" column should not be checked in the column list of \"Edit Tooltip\" dialog", () => uncheckedInColumnList(page, "RACE", el("\"Edit Tooltip\" dialog")));
    await session.step(39, "When user clicks on OK button in \"Edit Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Tooltip\" dialog")));
    await session.step(40, "Then the \"Edit Tooltip\" dialog should close", () => dialogCloses(page, "Edit Tooltip"));
    await session.step(42, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")), undefined, "Now the tooltip shows just those");
    await session.step(43, "Then the tooltip should show columns \"HEIGHT, WEIGHT, AGE, SEX, DIS_POP, STARTED, SEVERITY\"", () => tooltipColumns(page, "HEIGHT, WEIGHT, AGE, SEX, DIS_POP, STARTED, SEVERITY"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
