/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/06-a-file-in-the-context-panel.feature
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
import {clickOn, downloadContains, downloadThrough, expand, fileDownloaded, hoverText, shouldBe, videoSays, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/06-a-file-in-the-context-panel.feature", import.meta.url);
  test("A file's details and actions", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user watches downloads", () => watchDownloads(page));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(20, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(21, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(22, "And user clicks on Files---Demo---cars.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(23, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(24, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")));
    await session.step(25, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(26, "And the context panel should show \"cars.csv\"", () => contextPanelShows(page, "cars.csv"));
    await session.step(27, "When the video says \"The context panel shows the file you clicked\"", () => videoSays(page, "The context panel shows the file you clicked"));
    await session.step(28, "And user expands Details pane in context panel", () => expand(page, el("Details pane in context panel")));
    await session.step(29, "And user hovers over the text \"Size\" in Details pane in context panel", () => hoverText(page, "Size", el("Details pane in context panel")));
    await session.step(30, "And the video says \"Details: where it is, its size, when it changed\"", () => videoSays(page, "Details: where it is, its size, when it changed"));
    await session.step(31, "And user expands Actions pane in context panel", () => expand(page, el("Actions pane in context panel")));
    await session.step(32, "And the video says \"Actions: what you can do with it\"", () => videoSays(page, "Actions: what you can do with it"));
    await session.step(33, "And user downloads a file through \"Download\" label in Actions pane in context panel", () => downloadThrough(page, el("\"Download\" label in Actions pane in context panel")));
    await session.step(34, "Then a file \"cars.csv\" should have been downloaded", () => fileDownloaded(page, "cars.csv"));
    await session.step(35, "And the downloaded file \"cars.csv\" should contain text \"alfaromeo\"", () => downloadContains(page, "cars.csv", "alfaromeo"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "When the video says \"Download saves it to your computer\"", () => videoSays(page, "Download saves it to your computer"));
  });
});
