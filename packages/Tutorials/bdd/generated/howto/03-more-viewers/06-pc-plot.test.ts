/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/06-pc-plot.feature
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
import {check, clickOn, expand, hoverOver, shouldBe, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, readingAtLeast, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle, pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/06-pc-plot.feature", import.meta.url);
  test("Compare many columns with a PC plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(19, "When the video says \"A PC plot draws each row as a line across many columns\"", () => videoSays(page, "A PC plot draws each row as a line across many columns"));
    await session.step(20, "And user clicks on pc plot icon on toolbox", () => clickOn(page, el("pc plot icon on toolbox")));
    await session.step(21, "Then pc plot viewer should be added to the open tableview", () => viewerAdded(page, "pc plot"));
    await session.step(23, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(24, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(26, "When user picks \"Species\" in the \"color\" column selector of pc plot viewer", () => pickInColumnSelector(page, "Species", "color", el("pc plot viewer")), undefined, "Color the lines by a column");
    await session.step(28, "And user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")), undefined, "Drag an axis handle to keep only the rows in range");
    await session.step(29, "And user drags the max handle of the \"Petal.Length\" range slider of pc plot viewer by 300 pixels", () => dragRangeHandle(page, "max", "Petal.Length", el("pc plot viewer"), 300));
    await session.step(30, "Then fewer than 150 rows should pass the filter", () => filterPassesFewer(page, 150));
    await session.step(32, "When user clicks on settings icon of pc plot viewer", () => clickOn(page, el("settings icon of pc plot viewer")), undefined, "The rest is in the settings");
    await session.step(33, "Then the context panel should show \"PC Plot\"", () => contextPanelShows(page, "PC Plot"));
    await session.step(35, "When user checks \"Show Filtered Out Lines\" property in context panel", () => check(page, el("\"Show Filtered Out Lines\" property in context panel")), undefined, "Data: keep the filtered-out rows in view, in grey");
    await session.step(36, "Then the \"filtered out lines drawn\" reading of pc plot viewer should be at least 1", () => readingAtLeast(page, "filtered out lines drawn", el("pc plot viewer"), 1));
    await session.step(38, "When user expands \"Y Axis\" category in context panel", () => expand(page, el("\"Y Axis\" category in context panel")), undefined, "Y Axis: choose the columns to draw");
    await session.step(39, "And user clicks on \"...\" button in \"Column Names\" property in context panel", () => clickOn(page, el("\"...\" button in \"Column Names\" property in context panel")));
    await session.step(40, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(41, "When user toggles the \"col 1\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "col 1", el("\"Select columns...\" dialog")));
    await session.step(42, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(43, "Then the \"axis order\" reading of pc plot viewer should be \"Sepal.Length, Sepal.Width, Petal.Length, Petal.Width\"", () => readingReads(page, "axis order", el("pc plot viewer"), "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"));
    await session.step(45, "When user unchecks \"Normalize Each Column\" property in context panel", () => uncheck(page, el("\"Normalize Each Column\" property in context panel")), undefined, "All four are in cm, so one scale fits them all");
    await session.step(46, "Then the \"normalization\" reading of pc plot viewer should be \"global\"", () => readingReads(page, "normalization", el("pc plot viewer"), "global"));
    await session.step(48, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")), undefined, "Reset View clears the axis filters");
    await session.step(49, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
});
