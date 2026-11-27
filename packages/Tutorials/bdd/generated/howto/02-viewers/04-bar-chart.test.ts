/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/04-bar-chart.feature
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
import {pickAggregation, viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, hoverArea, noErrors, pickFromAreaContextMenu, pickFromOpenMenu, propertyShouldBe, readingIs, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/04-bar-chart.feature", import.meta.url);
  test("Configure a bar chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(24, "And bar chart viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("bar chart viewer"), 65));
    await session.step(25, "When the video says \"A bar chart compares one number across categories\"", () => videoSays(page, "A bar chart compares one number across categories"));
    await session.step(27, "And user picks \"DIS_POP\" in the \"split\" column selector of bar chart viewer", () => pickInColumnSelector(page, "DIS_POP", "split", el("bar chart viewer")), undefined, "The categories: the selector beside the bars");
    await session.step(28, "Then the \"bars\" reading of bar chart viewer should be 6", () => readingIs(page, "bars", el("bar chart viewer"), 6));
    await session.step(30, "When user picks the \"avg\" aggregation in the \"value\" column selector of bar chart viewer", () => pickAggregation(page, "avg", "value", el("bar chart viewer")), undefined, "The number: how to aggregate, and which column");
    await session.step(31, "And user picks \"WEIGHT\" in the \"value\" column selector of bar chart viewer", () => pickInColumnSelector(page, "WEIGHT", "value", el("bar chart viewer")));
    await session.step(32, "Then \"Value Aggr Type\" property of bar chart viewer should be \"avg\"", () => propertyShouldBe(page, "Value Aggr Type", el("bar chart viewer"), "avg"));
    await session.step(33, "And \"Value\" property of bar chart viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "Value", el("bar chart viewer"), "WEIGHT"));
    await session.step(35, "When user picks the \"count\" aggregation in the \"value\" column selector of bar chart viewer", () => pickAggregation(page, "count", "value", el("bar chart viewer")), undefined, "Stacking needs a count or a sum");
    await session.step(37, "And user picks \"SEX\" in the \"stack\" column selector of bar chart viewer", () => pickInColumnSelector(page, "SEX", "stack", el("bar chart viewer")), undefined, "Stack: split each bar by another category");
    await session.step(38, "Then the \"stack segments\" reading of bar chart viewer should be 12", () => readingIs(page, "stack segments", el("bar chart viewer"), 12));
    await session.step(39, "And the \"bar RA\" area of bar chart viewer should lie above the \"bar AS\" area", () => areaLies(page, "bar RA", el("bar chart viewer"), "above", "bar AS"));
    await session.step(41, "When user picks \"Order > by category\" from the context menu of the \"view\" area of bar chart viewer", () => pickFromAreaContextMenu(page, "Order > by category", "view", el("bar chart viewer")), undefined, "Right-click to sort the bars by name");
    await session.step(42, "And user picks \"Order > asc\" from the open menu", () => pickFromOpenMenu(page, "Order > asc"));
    await session.step(43, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(44, "Then the \"bar AS\" area of bar chart viewer should lie above the \"bar RA\" area", () => areaLies(page, "bar AS", el("bar chart viewer"), "above", "bar RA"));
    await session.step(46, "When user picks \"Orientation > vertical\" from the context menu of the \"view\" area of bar chart viewer", () => pickFromAreaContextMenu(page, "Orientation > vertical", "view", el("bar chart viewer")), undefined, "or to turn them upright");
    await session.step(47, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(48, "Then the \"bar AS\" area of bar chart viewer should lie to the left of the \"bar RA\" area", () => areaLiesBeside(page, "bar AS", el("bar chart viewer"), "left", "bar RA"));
    await session.step(50, "When user hovers over the \"x axis\" area of bar chart viewer", () => hoverArea(page, "x axis", el("bar chart viewer")), undefined, "The gear opens every setting, grouped in sections");
    await session.step(51, "And user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(52, "Then the context panel should show \"Bar chart\"", () => contextPanelShows(page, "Bar chart"));
    await session.step(54, "When user takes a snapshot of bar chart viewer", () => takeSnapshot(page, el("bar chart viewer")), undefined, "Y Axis: show each bar as 100%");
    await session.step(55, "And user expands \"Y Axis\" category in context panel", () => expand(page, el("\"Y Axis\" category in context panel")));
    await session.step(56, "And user checks \"Relative Values\" property in context panel", () => check(page, el("\"Relative Values\" property in context panel")));
    await session.step(57, "Then \"Relative Values\" property of bar chart viewer should be \"true\"", () => propertyShouldBe(page, "Relative Values", el("bar chart viewer"), "true"));
    await session.step(58, "And bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
