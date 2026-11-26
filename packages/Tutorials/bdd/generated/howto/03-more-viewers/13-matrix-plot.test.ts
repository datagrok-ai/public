/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/13-matrix-plot.feature
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
import {clickOn, expand, hoverOver, selectIn, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/13-matrix-plot.feature", import.meta.url);
  test("See every pair of columns with a matrix plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(19, "When the video says \"A matrix plot: every pair of numeric columns in one grid\"", () => videoSays(page, "A matrix plot: every pair of numeric columns in one grid"));
    await session.step(20, "And user clicks on matrix plot icon on toolbox", () => clickOn(page, el("matrix plot icon on toolbox")));
    await session.step(21, "Then matrix plot viewer should be added to the open tableview", () => viewerAdded(page, "matrix plot"));
    await session.step(23, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(24, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(25, "When the video says \"Each cell plots one column against another; the diagonal shows histograms\" at the top", () => videoSaysAt(page, "Each cell plots one column against another; the diagonal shows histograms", "top"));
    await session.step(26, "And user hovers over the \"cell body Petal.Length x Sepal.Length\" area of matrix plot viewer", () => hoverArea(page, "cell body Petal.Length x Sepal.Length", el("matrix plot viewer")));
    await session.step(28, "And user hovers over matrix plot viewer", () => hoverOver(page, el("matrix plot viewer")), undefined, "The rest is in the settings");
    await session.step(29, "And user clicks on settings icon of matrix plot viewer", () => clickOn(page, el("settings icon of matrix plot viewer")));
    await session.step(30, "Then the context panel should show \"Matrix plot\"", () => contextPanelShows(page, "Matrix plot"));
    await session.step(32, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: scatter plots instead of density plots");
    await session.step(33, "And user selects \"Scatter plot\" in \"Cell Plot Type\" property in context panel", () => selectIn(page, "Scatter plot", el("\"Cell Plot Type\" property in context panel")));
    await session.step(34, "Then the \"cell viewer type\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Scatter plot"));
    await session.step(36, "When user clicks on \"...\" button in \"X\" property in context panel", () => clickOn(page, el("\"...\" button in \"X\" property in context panel")), undefined, "Data: choose the columns for each axis");
    await session.step(37, "And user toggles the \"col 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "col 1", el("\"Select columns...\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(39, "And user clicks on \"...\" button in \"Y\" property in context panel", () => clickOn(page, el("\"...\" button in \"Y\" property in context panel")));
    await session.step(40, "And user toggles the \"col 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "col 1", el("\"Select columns...\" dialog")));
    await session.step(41, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(42, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(44, "When user clicks on \"Scatter plot\" tab in context panel", () => clickOn(page, el("\"Scatter plot\" tab in context panel")), undefined, "The second tab sets up the plot inside every cell");
    await session.step(45, "And user expands \"Color\" category in context panel", () => expand(page, el("\"Color\" category in context panel")));
    await session.step(46, "And user selects \"Species\" in \"Color\" property in context panel", () => selectIn(page, "Species", el("\"Color\" property in context panel")));
    await session.step(47, "Then \"colorColumnName\" inner property of matrix plot viewer should be \"Species\"", () => innerPropertyShouldBe(page, "colorColumnName", el("matrix plot viewer"), "Species"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
