/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/01-find-your-way-around.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/01-find-your-way-around.feature", import.meta.url);
  test("Find your way around Datagrok", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(14, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(15, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(16, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(17, "When the video says \"Your data opens in the middle. The panels sit around it\"", () => videoSays(page, "Your data opens in the middle. The panels sit around it"));
    await session.step(19, "And user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Browse: your files, databases, queries and projects");
    await session.step(20, "Then browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(22, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "The toolbox: viewers, filters and layouts for this table");
    await session.step(23, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(25, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The status bar turns panels on and off");
    await session.step(27, "And user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")), undefined, "Click a column to see its details on the right");
    await session.step(28, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(29, "When user clicks on the \"header RACE\" area of grid", () => clickArea(page, "header RACE", el("grid")));
    await session.step(30, "Then the context panel should show \"RACE\"", () => contextPanelShows(page, "RACE"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
});
