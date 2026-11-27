/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/03-trellis-plot.feature
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
import {check, clickOn, expand, hoverOver, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, pickInAreaSelector, pickInColumnSelector, pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/03-trellis-plot.feature", import.meta.url);
  test("Small multiples with a trellis plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(14, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(15, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(16, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(17, "When the video says \"A trellis plot draws one small chart per category\"", () => videoSays(page, "A trellis plot draws one small chart per category"));
    await session.step(18, "And user clicks on trellis plot icon on toolbox", () => clickOn(page, el("trellis plot icon on toolbox")));
    await session.step(19, "Then trellis plot viewer should be added to the open tableview", () => viewerAdded(page, "trellis plot"));
    await session.step(21, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(22, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(24, "When user picks \"SEX\" in the column selector at the \"x selector 1\" area of trellis plot viewer", () => pickInAreaSelector(page, "SEX", "x selector 1", el("trellis plot viewer")), undefined, "One column splits the cells across");
    await session.step(26, "And user picks \"RACE\" in the column selector at the \"y selector 1\" area of trellis plot viewer", () => pickInAreaSelector(page, "RACE", "y selector 1", el("trellis plot viewer")), undefined, "Another splits them down");
    await session.step(27, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(29, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")), undefined, "Choose the chart that goes in each cell");
    await session.step(30, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    await session.step(32, "When user picks \"AGE\" in the \"color\" column selector of trellis plot viewer", () => pickInColumnSelector(page, "AGE", "color", el("trellis plot viewer")), undefined, "Its columns apply to every cell");
    await session.step(34, "When user hovers over trellis plot viewer", () => hoverOver(page, el("trellis plot viewer")), undefined, "The rest is in the settings");
    await session.step(35, "And user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
    await session.step(36, "Then the context panel should show \"Trellis plot\"", () => contextPanelShows(page, "Trellis plot"));
    await session.step(38, "When user expands \"Axes\" category in context panel", () => expand(page, el("\"Axes\" category in context panel")), undefined, "Axes: one scale for every cell, so they compare");
    await session.step(39, "And user checks \"Global Scale\" property in context panel", () => check(page, el("\"Global Scale\" property in context panel")));
    await session.step(41, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: click a cell to filter the table to it");
    await session.step(42, "And user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
    await session.step(43, "And user clicks on the \"cell M | Black\" area of trellis plot viewer", () => clickArea(page, "cell M | Black", el("trellis plot viewer")));
    await session.step(44, "Then 53 rows should pass the filter", () => filterPasses(page, 53));
    await session.step(46, "When user clicks on \"Scatter plot\" tab in context panel", () => clickOn(page, el("\"Scatter plot\" tab in context panel")), undefined, "The second tab holds the settings of the chart inside");
    await session.step(47, "And user expands \"Statistics\" category in context panel", () => expand(page, el("\"Statistics\" category in context panel")));
    await session.step(48, "And user checks \"Show Regression Line\" property in context panel", () => check(page, el("\"Show Regression Line\" property in context panel")));
    await session.step(49, "Then \"showRegressionLine\" inner property of trellis plot viewer should be \"true\"", () => innerPropertyShouldBe(page, "showRegressionLine", el("trellis plot viewer"), "true"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
});
