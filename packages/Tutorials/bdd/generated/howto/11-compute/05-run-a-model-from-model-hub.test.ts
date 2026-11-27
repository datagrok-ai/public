/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/05-run-a-model-from-model-hub.feature
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
import {clickOn, doubleClickOn, dragSliderTo, enterInto, hoverOver, shouldBe, shouldHaveValue, shouldHaveValueBetween, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/05-run-a-model-from-model-hub.feature", import.meta.url);
  test("Run a model from Model Hub", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(16, "And the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Model Hub is in Browse > Apps");
    await session.step(23, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(24, "And user double-clicks on Model-Hub gallery card", () => doubleClickOn(page, el("Model-Hub gallery card")));
    await session.step(25, "Then \"Ball Flight Simulation\" link should be visible", () => shouldBe(page, el("\"Ball Flight Simulation\" link"), "visible"));
    await session.step(27, "When user hovers over \"Ball Flight Simulation\" link", () => hoverOver(page, el("\"Ball Flight Simulation\" link")), undefined, "Each card is a model; point at one to read about it");
    await session.step(29, "And user double-clicks on \"Ball Flight Simulation\" link", () => doubleClickOn(page, el("\"Ball Flight Simulation\" link")), undefined, "Double-click to run it");
    await session.step(30, "Then \"Angle\" input should be visible", () => shouldBe(page, el("\"Angle\" input"), "visible"));
    await session.step(31, "When the video says \"Inputs on the left; the results update as you change them\"", () => videoSays(page, "Inputs on the left; the results update as you change them"));
    await session.step(32, "And user enters \"30\" into \"Angle\" input", () => enterInto(page, "30", el("\"Angle\" input")));
    await session.step(33, "Then \"Angle\" input should have the value \"30\"", () => shouldHaveValue(page, el("\"Angle\" input"), "30"));
    await session.step(34, "When user drags the slider of \"Velocity\" input to 58", () => dragSliderTo(page, el("\"Velocity\" input"), 58));
    await session.step(35, "Then \"Velocity\" input should have a value between 57 and 59", () => shouldHaveValueBetween(page, el("\"Velocity\" input"), 57, 59));
    await session.step(36, "When the video says \"Max distance and height, at the bottom, follow each change\" at the top", () => videoSaysAt(page, "Max distance and height, at the bottom, follow each change", "top"));
    await session.step(37, "And the video says \"Fitting and sensitivity analysis are on its ribbon too\" at the top", () => videoSaysAt(page, "Fitting and sensitivity analysis are on its ribbon too", "top"));
    await session.step(38, "And user hovers over \"Fit inputs\" icon", () => hoverOver(page, el("\"Fit inputs\" icon")));
    await session.step(39, "Then no errors should have been logged", () => noErrors(page));
  });
});
