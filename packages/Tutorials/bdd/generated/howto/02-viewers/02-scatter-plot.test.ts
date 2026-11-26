/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/02-scatter-plot.feature
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
import {check, clickOn, expand, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaWider, closeContextMenu, dragAcrossArea, hasArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, propertyShouldContain, readingAtLeast, takeSnapshot, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace, hoverEmptySpace, pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/02-scatter-plot.feature", import.meta.url);
  test("Configure a scatter plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","wheel.base"],["Y","length"]]), [["X","wheel.base"],["Y","length"]]);
    await session.step(27, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(28, "When the video says \"A scatter plot draws one dot per row\"", () => videoSays(page, "A scatter plot draws one dot per row"));
    await session.step(30, "And user picks \"horsepower\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "horsepower", "x", el("scatter plot viewer")), undefined, "Pick X and Y right on the axes");
    await session.step(31, "And user picks \"price\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "price", "y", el("scatter plot viewer")));
    await session.step(32, "Then \"X\" property of scatter plot viewer should be \"horsepower\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "horsepower"));
    await session.step(33, "And \"Y\" property of scatter plot viewer should be \"price\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "price"));
    await session.step(35, "When user picks \"city.mpg\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "city.mpg", "color", el("scatter plot viewer")), undefined, "Color and Size: the selectors at the top right");
    await session.step(36, "Then scatter plot viewer should have a \"color scale\" area", () => hasArea(page, el("scatter plot viewer"), "color scale"));
    await session.step(37, "When user picks \"curb.weight\" in the \"size\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "curb.weight", "size", el("scatter plot viewer")));
    await session.step(38, "Then \"Size\" property of scatter plot viewer should be \"curb.weight\"", () => propertyShouldBe(page, "Size", el("scatter plot viewer"), "curb.weight"));
    await session.step(40, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")), undefined, "The wheel zooms, and the table keeps the rows in view");
    await session.step(41, "Then fewer than 30 rows should pass the filter", () => filterPassesFewer(page, 30));
    await session.step(43, "When user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")), undefined, "Drag to move around");
    await session.step(45, "And user double-clicks on empty plot space of scatter plot viewer", () => doubleClickEmptySpace(page, el("scatter plot viewer")), undefined, "Double-click an empty spot to see all the points again");
    await session.step(46, "Then 30 rows should pass the filter", () => filterPasses(page, 30));
    await session.step(48, "When user picks \"Tools > Show Regression Line\" from the context menu of the \"view\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Show Regression Line", "view", el("scatter plot viewer")), undefined, "Right-click for more, like a regression line");
    await session.step(49, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(50, "Then scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
    await session.step(52, "When user hovers over empty plot space of scatter plot viewer", () => hoverEmptySpace(page, el("scatter plot viewer")), undefined, "The gear opens every setting, grouped in sections");
    await session.step(53, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(54, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(56, "When user expands \"Marker\" category in context panel", () => expand(page, el("\"Marker\" category in context panel")), undefined, "Marker: the shape of the dots");
    await session.step(57, "And user selects \"square\" in \"Marker Type\" property in context panel", () => selectIn(page, "square", el("\"Marker Type\" property in context panel")));
    await session.step(58, "Then \"Marker Type\" property of scatter plot viewer should be \"square\"", () => propertyShouldBe(page, "Marker Type", el("scatter plot viewer"), "square"));
    await session.step(60, "When user expands \"Labels\" category in context panel", () => expand(page, el("\"Labels\" category in context panel")), undefined, "Labels: write a column's value next to each dot");
    await session.step(61, "And user clicks on \"...\" button in \"Label Columns\" property in context panel", () => clickOn(page, el("\"...\" button in \"Label Columns\" property in context panel")));
    await session.step(62, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(63, "When user toggles the \"model\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "model", el("\"Select columns...\" dialog")));
    await session.step(64, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(65, "Then \"Label Columns\" property of scatter plot viewer should contain \"model\"", () => propertyShouldContain(page, "Label Columns", el("scatter plot viewer"), "model"));
    await session.step(66, "And the \"labels shown\" reading of scatter plot viewer should be at least 10", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 10));
    await session.step(68, "When user expands \"X Axis\" category in context panel", () => expand(page, el("\"X Axis\" category in context panel")), undefined, "X Axis: a histogram along the axis");
    await session.step(69, "And user checks \"Show X Histogram\" property in context panel", () => check(page, el("\"Show X Histogram\" property in context panel")));
    await session.step(70, "Then scatter plot viewer should have an \"x histogram\" area", () => hasArea(page, el("scatter plot viewer"), "x histogram"));
    await session.step(72, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")), undefined, "Statistics: what to show beside the fit");
    await session.step(73, "And user expands \"Statistics\" category in context panel", () => expand(page, el("\"Statistics\" category in context panel")));
    await session.step(74, "And user checks \"Show Pearson Correlation\" property in context panel", () => check(page, el("\"Show Pearson Correlation\" property in context panel")));
    await session.step(75, "Then the \"regression stats\" area of scatter plot viewer should be wider than before", () => areaWider(page, "regression stats", el("scatter plot viewer")));
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
  });
});
