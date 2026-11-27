/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/07-viewer-settings.feature
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
import {clearField, clickOn, enterInto, expand, selectIn, shouldHaveText, typeInto, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, legendSide, noErrors, pickColorSwatch, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace, hoverEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/07-viewer-settings.feature", import.meta.url);
  test("Change a viewer's settings", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["Color","RACE"]]), [["Color","RACE"]]);
    await session.step(23, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(24, "When the video says \"Every viewer has its settings in the context panel\"", () => videoSays(page, "Every viewer has its settings in the context panel"));
    await session.step(26, "And user hovers over empty plot space of scatter plot viewer", () => hoverEmptySpace(page, el("scatter plot viewer")), undefined, "The gear on the title bar opens them");
    await session.step(27, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(28, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(30, "When user expands \"Description\" category in context panel", () => expand(page, el("\"Description\" category in context panel")), undefined, "Sections group the settings. Description: title and notes");
    await session.step(31, "And user enters \"Height and weight\" into \"Title\" property in context panel", () => enterInto(page, "Height and weight", el("\"Title\" property in context panel")));
    await session.step(32, "Then title of scatter plot viewer should have text \"Height and weight\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Height and weight"));
    await session.step(33, "When user enters \"5,850 patients of the demog study\" into \"Description\" property in context panel", () => enterInto(page, "5,850 patients of the demog study", el("\"Description\" property in context panel")));
    await session.step(34, "Then description of scatter plot viewer should have text \"5,850 patients of the demog study\"", () => shouldHaveText(page, el("description of scatter plot viewer"), "5,850 patients of the demog study"));
    await session.step(36, "When user expands \"Legend\" category in context panel", () => expand(page, el("\"Legend\" category in context panel")), undefined, "Legend: where it goes");
    await session.step(37, "And user selects \"Bottom\" in \"Legend Position\" property in context panel", () => selectIn(page, "Bottom", el("\"Legend Position\" property in context panel")));
    await session.step(38, "Then the legend of scatter plot viewer should be on the bottom", () => legendSide(page, el("scatter plot viewer"), "bottom"));
    await session.step(40, "When user clicks on context panel search", () => clickOn(page, el("context panel search")), undefined, "Search finds a setting in any section");
    await session.step(41, "And user types \"grid lines\" into context panel search", () => typeInto(page, "grid lines", el("context panel search")));
    await session.step(42, "And user unchecks \"Show Vertical Grid Lines\" property in context panel", () => uncheck(page, el("\"Show Vertical Grid Lines\" property in context panel")));
    await session.step(43, "And user unchecks \"Show Horizontal Grid Lines\" property in context panel", () => uncheck(page, el("\"Show Horizontal Grid Lines\" property in context panel")));
    await session.step(44, "Then \"Show Vertical Grid Lines\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Vertical Grid Lines", el("scatter plot viewer"), "false"));
    await session.step(45, "And \"Show Horizontal Grid Lines\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Horizontal Grid Lines", el("scatter plot viewer"), "false"));
    await session.step(47, "When user clears context panel search", () => clearField(page, el("context panel search")), undefined, "Style: colors, fonts and sizes");
    await session.step(48, "And user types \"back color\" into context panel search", () => typeInto(page, "back color", el("context panel search")));
    await session.step(49, "And user clicks on editor of \"Back Color\" property in context panel", () => clickOn(page, el("editor of \"Back Color\" property in context panel")));
    await session.step(50, "And user picks the color \"#9EDAE5\" in the color picker", () => pickColorSwatch(page, "#9EDAE5"));
    await session.step(51, "And user clicks on empty plot space of scatter plot viewer", () => clickEmptySpace(page, el("scatter plot viewer")));
    await session.step(52, "Then the \"view\" area of scatter plot viewer should contain the color \"#9EDAE5\"", () => areaColor(page, "view", el("scatter plot viewer"), "#9EDAE5"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
});
