/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/01-run-a-model-in-diff-studio.feature
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
import {clickOn, doubleClickOn, dragSliderTo, enterInto, hoverOver, shouldBe, shouldHaveValue, shouldHaveValueBetween, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/01-run-a-model-in-diff-studio.feature", import.meta.url);
  test("Run a model in Diff Studio", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When the video says \"Diff Studio is one of the apps in Browse\"", () => videoSays(page, "Diff Studio is one of the apps in Browse"));
    await session.step(23, "And user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(24, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(25, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(26, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(28, "When user hovers over \"Lotka-Volterra\" model card", () => hoverOver(page, el("\"Lotka-Volterra\" model card")), undefined, "Templates to start from, and a library of ready models");
    await session.step(29, "And user double-clicks on \"Lotka-Volterra\" model card", () => doubleClickOn(page, el("\"Lotka-Volterra\" model card")));
    await session.step(30, "Then the \"Lotka-Volterra\" view should be current", () => viewIsCurrent(page, "Lotka-Volterra"));
    await session.step(31, "And \"Prey\" input should be visible", () => shouldBe(page, el("\"Prey\" input"), "visible"));
    await session.step(32, "When the video says \"Inputs on the left, the solution on the right\"", () => videoSays(page, "Inputs on the left, the solution on the right"));
    await session.step(34, "And user enters \"30\" into \"Prey\" input", () => enterInto(page, "30", el("\"Prey\" input")), undefined, "Change a starting value: the curves are recomputed at once");
    await session.step(35, "Then \"Prey\" input should have the value \"30\"", () => shouldHaveValue(page, el("\"Prey\" input"), "30"));
    await session.step(37, "When user hovers over \"alpha\" input", () => hoverOver(page, el("\"alpha\" input")), undefined, "Point at an input to see what it stands for");
    await session.step(39, "And user drags the slider of \"alpha\" input to 0.6", () => dragSliderTo(page, el("\"alpha\" input"), 0.6), undefined, "Drag its slider and watch the curves follow");
    await session.step(40, "Then \"alpha\" input should have a value between 0.5 and 0.7", () => shouldHaveValueBetween(page, el("\"alpha\" input"), 0.5, 0.7));
    await session.step(42, "When user enters \"150\" into \"Finish\" input", () => enterInto(page, "150", el("\"Finish\" input")), undefined, "Finish: how far in time to solve");
    await session.step(43, "Then \"Finish\" input should have the value \"150\"", () => shouldHaveValue(page, el("\"Finish\" input"), "150"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
