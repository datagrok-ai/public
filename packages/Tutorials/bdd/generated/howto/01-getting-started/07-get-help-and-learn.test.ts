/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/07-get-help-and-learn.feature
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
import {closeTutorial, stepDone, tutorialNotCompleted, tutorialProgress, tutorialsClosed} from '../../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKey, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelClosed, openDataset, simpleModeOff, standServesHelp, toolboxPaneHidden, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/07-get-help-and-learn.feature", import.meta.url);
  test("Get help and learn", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the stand serves the help pages", () => standServesHelp(page));
    await session.step(22, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(23, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(24, "And the \"Grid Customization\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Grid Customization"));
    await session.step(25, "And the Tutorials app is closed", () => tutorialsClosed(page));
    await session.step(26, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(27, "And user adds scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(28, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(29, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "When user clicks on \"Context Help\" status bar toggle", () => clickOn(page, el("\"Context Help\" status bar toggle")), undefined, "Context Help: the help page of what you click");
    await session.step(32, "Then help panel should be visible", () => shouldBe(page, el("help panel"), "visible"));
    await session.step(34, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")), undefined, "Click a column: help on tables and columns");
    await session.step(35, "Then help panel should contain text \"Dataframe\"", () => shouldContainText(page, el("help panel"), "Dataframe"));
    await session.step(37, "When user clicks on empty plot space of scatter plot viewer", () => clickEmptySpace(page, el("scatter plot viewer")), undefined, "Click a viewer: help on that viewer");
    await session.step(38, "Then help panel should contain text \"scatterplot\"", () => shouldContainText(page, el("help panel"), "scatterplot"));
    await session.step(40, "When user clicks on \"Context Help\" status bar toggle", () => clickOn(page, el("\"Context Help\" status bar toggle")), undefined, "The same toggle closes it");
    await session.step(41, "Then help panel should be hidden", () => shouldBe(page, el("help panel"), "hidden"));
    await session.step(43, "When user clicks on \"Help\" tab", () => clickOn(page, el("\"Help\" tab")), undefined, "The sidebar's Help: docs, community, videos, feedback");
    await session.step(44, "And user hovers over \"Help\" menu item inside context menu", () => hoverOver(page, el("\"Help\" menu item inside context menu")));
    await session.step(46, "And user hovers over \"Community\" menu item inside context menu", () => hoverOver(page, el("\"Community\" menu item inside context menu")), undefined, "Each one opens in a new browser tab");
    await session.step(47, "And user hovers over \"YouTube\" menu item inside context menu", () => hoverOver(page, el("\"YouTube\" menu item inside context menu")));
    await session.step(48, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(49, "When the video says \"To learn by doing, take a tutorial\"", () => videoSays(page, "To learn by doing, take a tutorial"));
    await session.step(51, "And user clicks on \"Menu\" tab", () => clickOn(page, el("\"Menu\" tab")), undefined, "The menu behind the logo: Help > Tutorials");
    await session.step(52, "And user picks \"Help > Tutorials\" from the open menu", () => pickFromOpenMenu(page, "Help > Tutorials"));
    await session.step(53, "Then Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(55, "When user hovers over \"Grid Customization\" tutorial card", () => hoverOver(page, el("\"Grid Customization\" tutorial card")), undefined, "Tutorials come in tracks, from data analysis to chemistry");
    await session.step(57, "And user clicks on \"Grid Customization\" tutorial card", () => clickOn(page, el("\"Grid Customization\" tutorial card")), undefined, "Click one to start it");
    await session.step(58, "Then the tutorial progress should be 1 of 30", () => tutorialProgress(page, 1, 30));
    await session.step(60, "When user hovers over tutorial bar", () => hoverOver(page, el("tutorial bar")), undefined, "The bar on top says what to do next");
    await session.step(62, "And user clicks on all steps icon", () => clickOn(page, el("all steps icon")), undefined, "This icon lists all the steps");
    await session.step(63, "Then Tutorials panel should be visible", () => shouldBe(page, el("Tutorials panel"), "visible"));
    await session.step(65, "When user clicks on the \"cell 1 of AGE\" area of grid", () => clickArea(page, "cell 1 of AGE", el("grid")), undefined, "Do it, and the step is ticked off");
    await session.step(66, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(67, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(68, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(69, "Then the tutorial step \"Go a few rows down in the grid\" should be done", () => stepDone(page, "Go a few rows down in the grid"));
    await session.step(71, "When user closes the tutorial", () => closeTutorial(page), undefined, "Close it any time from its bar");
    await session.step(72, "Then no errors should have been logged", () => noErrors(page));
  });
});
