/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/03-label-points-with-structures.feature
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
import {selectOnAxes} from '../../../bindings/howto-annotations.js';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, noErrors, pickFromAreaContextMenu, propertyShouldBe, propertyShouldContain, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {hoverEmptySpace, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/03-label-points-with-structures.feature", import.meta.url);
  test("Label points with structures", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","Chemical Space X"],["Y","Chemical Space Y"]]), [["X","Chemical Space X"],["Y","Chemical Space Y"]]);
    await session.step(26, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(27, "When the video says \"Labels put a column's values next to the points\"", () => videoSays(page, "Labels put a column's values next to the points"));
    await session.step(29, "And user hovers over empty plot space of scatter plot viewer", () => hoverEmptySpace(page, el("scatter plot viewer")), undefined, "The gear icon opens the plot's settings");
    await session.step(30, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(31, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(33, "When user expands \"Labels\" category in context panel", () => expand(page, el("\"Labels\" category in context panel")), undefined, "Labels > Label Columns: pick what to show");
    await session.step(34, "And user clicks on \"...\" button in \"Label Columns\" property in context panel", () => clickOn(page, el("\"...\" button in \"Label Columns\" property in context panel")));
    await session.step(35, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(36, "When user toggles the \"Structure\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Structure", el("\"Select columns...\" dialog")));
    await session.step(37, "And user toggles the \"Id\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Id", el("\"Select columns...\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(40, "Then \"Label Columns\" property of scatter plot viewer should contain \"Structure\"", () => propertyShouldContain(page, "Label Columns", el("scatter plot viewer"), "Structure"), undefined, "A molecule column is drawn as its structure");
    await session.step(41, "And \"Label Columns\" property of scatter plot viewer should contain \"Id\"", () => propertyShouldContain(page, "Label Columns", el("scatter plot viewer"), "Id"));
    await session.step(42, "And the \"labels shown\" reading of scatter plot viewer should be at least 5", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 5));
    await session.step(44, "When user enters \"1.6\" into \"Label Content Size\" property in context panel", () => enterInto(page, "1.6", el("\"Label Content Size\" property in context panel")), undefined, "Label Content Size: bigger structures");
    await session.step(45, "Then \"Label Content Size\" property of scatter plot viewer should be \"1.6\"", () => propertyShouldBe(page, "Label Content Size", el("scatter plot viewer"), "1.6"));
    await session.step(47, "When user picks \"Labels > Show Labels For > Selected\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Labels > Show Labels For > Selected", "view", el("scatter plot viewer")), undefined, "The plot's menu has them too. Show Labels For: Selected");
    await session.step(48, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(49, "Then \"Show Labels For\" property of scatter plot viewer should be \"Selected\"", () => propertyShouldBe(page, "Show Labels For", el("scatter plot viewer"), "Selected"));
    await session.step(50, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
    await session.step(52, "When user drags a selection box on scatter plot viewer from 6.3, 6.5 to 7.5, 5.4 on its axes", () => selectOnAxes(page, el("scatter plot viewer"), 6.3, 6.5, 7.5, 5.4), undefined, "Shift+drag round some points: their labels appear");
    await session.step(53, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(54, "And the \"labels shown\" reading of scatter plot viewer should be 4", () => readingIs(page, "labels shown", el("scatter plot viewer"), 4));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
});
