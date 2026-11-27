/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/05-compare-groups-with-a-control.feature
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
import {clickOn, hoverOver, selectIn, shouldContainText, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hoverArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/05-compare-groups-with-a-control.feature", import.meta.url);
  test("Compare each group with a control", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","HEIGHT"],["Category 1","RACE"]]), [["Value","HEIGHT"],["Category 1","RACE"]]);
    await session.step(23, "And box plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("box plot viewer"), 65));
    await session.step(24, "When the video says \"Test each group against one control group\"", () => videoSays(page, "Test each group against one control group"));
    await session.step(26, "And user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")), undefined, "The p-value at the top tests all groups at once");
    await session.step(27, "Then tooltip should contain text \"Alexander and Govern\"", () => shouldContainText(page, el("tooltip"), "Alexander and Govern"));
    await session.step(29, "When user clicks on show group stats icon in box plot viewer", () => clickOn(page, el("show group stats icon in box plot viewer")), undefined, "The icon beside it compares each group with one");
    await session.step(30, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
    await session.step(32, "And user selects \"Caucasian\" in control group choice input in box plot viewer", () => selectIn(page, "Caucasian", el("control group choice input in box plot viewer")), undefined, "Pick the control group");
    await session.step(33, "Then box plot viewer should have a \"p value of Asian\" area", () => hasArea(page, el("box plot viewer"), "p value of Asian"));
    await session.step(34, "When the video says \"Under each group: *** differs from the control, ns does not\" at the top", () => videoSaysAt(page, "Under each group: *** differs from the control, ns does not", "top"));
    await session.step(35, "And the video says \"Point at a mark to see its p-value\" at the top", () => videoSaysAt(page, "Point at a mark to see its p-value", "top"));
    await session.step(36, "And user hovers over the \"p value of Black\" area of box plot viewer", () => hoverArea(page, "p value of Black", el("box plot viewer")));
    await session.step(37, "Then tooltip should contain text \"Black\"", () => shouldContainText(page, el("tooltip"), "Black"));
    await session.step(39, "When user picks \"Add Control Comparisons Table\" from the context menu of the \"group comparison\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Add Control Comparisons Table", "group comparison", el("box plot viewer")), undefined, "Right-click the comparison to get the numbers as a table");
    await session.step(41, "Then table \"Control Comparisons: HEIGHT by RACE vs Caucasian\" should be open", () => tableOpen(page, "Control Comparisons: HEIGHT by RACE vs Caucasian"), undefined, "Asian and Other differ from Caucasian, Black does not");
    await session.step(42, "And table \"Control Comparisons: HEIGHT by RACE vs Caucasian\" should have 3 rows", () => tableRows(page, "Control Comparisons: HEIGHT by RACE vs Caucasian", 3));
    await session.step(43, "And the value of \"Group\" column in row 1 should be \"Asian\"", () => valueInRow(page, "Group", 1, "Asian"));
    await session.step(44, "And the value of \"Conclusion\" column in row 1 should be \"Significant\"", () => valueInRow(page, "Conclusion", 1, "Significant"));
    await session.step(45, "And the value of \"Group\" column in row 2 should be \"Black\"", () => valueInRow(page, "Group", 2, "Black"));
    await session.step(46, "And the value of \"Conclusion\" column in row 2 should be \"Not significant\"", () => valueInRow(page, "Conclusion", 2, "Not significant"));
    await session.step(47, "And the value of \"Group\" column in row 3 should be \"Other\"", () => valueInRow(page, "Group", 3, "Other"));
    await session.step(48, "And the value of \"Conclusion\" column in row 3 should be \"Significant\"", () => valueInRow(page, "Conclusion", 3, "Significant"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
