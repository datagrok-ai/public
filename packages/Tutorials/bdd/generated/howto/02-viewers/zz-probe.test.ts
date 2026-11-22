/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/zz-probe.feature
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
import {clickOn, hoverOver, pressKey, shouldBe, waitSeconds} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/zz-probe.feature", import.meta.url);
  test("Probe full screen", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And simple mode is off", () => simpleModeOff(page));
    await session.step(10, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(11, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(12, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(13, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(14, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(15, "When user hovers over title of scatter plot viewer", () => hoverOver(page, el("title of scatter plot viewer")));
    await session.step(16, "And user clicks on \"expand-arrows\" icon in \"Scatter plot\" dock panel", () => clickOn(page, el("\"expand-arrows\" icon in \"Scatter plot\" dock panel")));
    await session.step(17, "And user waits 2 seconds", () => waitSeconds(page, 2));
    await session.step(18, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(19, "And user waits 2 seconds", () => waitSeconds(page, 2));
    await session.step(20, "And user hovers over title of scatter plot viewer", () => hoverOver(page, el("title of scatter plot viewer")));
    await session.step(21, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
  });
});
