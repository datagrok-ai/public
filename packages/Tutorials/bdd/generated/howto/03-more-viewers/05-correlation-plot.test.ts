/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/05-correlation-plot.feature
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
import {doubleClickArea, hoverArea, noErrors, pointerAway, readingReads, viewerAdded, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/05-correlation-plot.feature", import.meta.url);
  test("Find related columns with a correlation plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(19, "When the video says \"A correlation plot compares every pair of numeric columns\"", () => videoSays(page, "A correlation plot compares every pair of numeric columns"));
    await session.step(20, "And user clicks on correlation plot icon on toolbox", () => clickOn(page, el("correlation plot icon on toolbox")));
    await session.step(21, "Then correlation plot viewer should be added to the open tableview", () => viewerAdded(page, "correlation plot"));
    await session.step(23, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(24, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(25, "When the video says \"Red: the two rise together. Blue: one falls as the other rises\" at the bottom", () => videoSaysAt(page, "Red: the two rise together. Blue: one falls as the other rises", "bottom"));
    await session.step(27, "And user hovers over the \"cell horsepower x price\" area of correlation plot viewer", () => hoverArea(page, "cell horsepower x price", el("correlation plot viewer")), undefined, "Hover a cell to see that pair as a scatter plot");
    await session.step(29, "When user hovers over correlation plot viewer", () => hoverOver(page, el("correlation plot viewer")), undefined, "The rest is in the settings");
    await session.step(30, "And user clicks on settings icon of correlation plot viewer", () => clickOn(page, el("settings icon of correlation plot viewer")));
    await session.step(31, "Then the context panel should show \"Correlation plot\"", () => contextPanelShows(page, "Correlation plot"));
    await session.step(33, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: Spearman compares ranks, so outliers count less");
    await session.step(34, "And user selects \"Spearman\" in \"Correlation Type\" property in context panel", () => selectIn(page, "Spearman", el("\"Correlation Type\" property in context panel")));
    await session.step(35, "Then the \"correlation type\" reading of correlation plot viewer should be \"Spearman\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Spearman"));
    await session.step(37, "When user clicks on \"...\" button in \"Y\" property in context panel", () => clickOn(page, el("\"...\" button in \"Y\" property in context panel")), undefined, "Data: keep just mileage and price on one axis");
    await session.step(38, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(39, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(40, "And user toggles the \"city.mpg\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "city.mpg", el("\"Select columns...\" dialog")));
    await session.step(41, "And user toggles the \"price\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "price", el("\"Select columns...\" dialog")));
    await session.step(42, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(43, "Then the \"y columns\" reading of correlation plot viewer should be \"city.mpg, price\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "city.mpg, price"));
    await session.step(45, "When user double-clicks on the \"cell horsepower x price\" area of correlation plot viewer", () => doubleClickArea(page, "cell horsepower x price", el("correlation plot viewer")), undefined, "Double-click a cell to open the pair as a scatter plot");
    await session.step(46, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(47, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
    await session.step(48, "Then no errors should have been logged", () => noErrors(page));
  });
});
