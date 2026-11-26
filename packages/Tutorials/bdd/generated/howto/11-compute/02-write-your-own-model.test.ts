/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/02-write-your-own-model.feature
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
import {clickOn, doubleClickOn, dragBy, dragSliderTo, insertLineAfter, replaceLine, shouldBe, shouldHaveValueBetween, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/02-write-your-own-model.feature", import.meta.url);
  test("Write your own model", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Diff Studio is in Browse > Apps");
    await session.step(25, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(26, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(27, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(29, "When user clicks on Create button", () => clickOn(page, el("Create button")), undefined, "Create starts a new model from a template");
    await session.step(30, "Then code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
    await session.step(32, "When user drags equations editor divider by 300 and 0 pixels", () => dragBy(page, el("equations editor divider"), 300, 0), undefined, "Drag the divider to widen the editor");
    await session.step(33, "When the video says \"A model is plain text: equations, time, starting values, parameters\"", () => videoSays(page, "A model is plain text: equations, time, starting values, parameters"));
    await session.step(35, "And user replaces the line starting with \"#name\" in code editor with \"#name: Population growth\"", () => replaceLine(page, "#name", el("code editor"), "#name: Population growth"), undefined, "Give it a name");
    await session.step(37, "And user replaces the line starting with \"dy/dt\" in code editor with \"dN/dt = r * N * (1 - N / K)\"", () => replaceLine(page, "dy/dt", el("code editor"), "dN/dt = r * N * (1 - N / K)"), undefined, "One line per equation: logistic growth");
    await session.step(39, "And user replaces the line starting with \"y = 0\" in code editor with \"N = 10 {min: 1; max: 100} [Starting population]\"", () => replaceLine(page, "y = 0", el("code editor"), "N = 10 {min: 1; max: 100} [Starting population]"), undefined, "The starting value, with a range for its slider and a hint");
    await session.step(41, "And user puts \"#parameters:\" on a new line after the line starting with \"N = 10\" in code editor", () => insertLineAfter(page, "#parameters:", "N = 10", el("code editor")), undefined, "Parameters get inputs of their own");
    await session.step(42, "And user puts \"  r = 0.5 {min: 0.1; max: 2} [Growth rate]\" on a new line after the line starting with \"#parameters\" in code editor", () => insertLineAfter(page, "  r = 0.5 {min: 0.1; max: 2} [Growth rate]", "#parameters", el("code editor")));
    await session.step(43, "And user puts \"  K = 1000 {min: 100; max: 2000} [Carrying capacity]\" on a new line after the line starting with \"r = 0.5\" in code editor", () => insertLineAfter(page, "  K = 1000 {min: 100; max: 2000} [Carrying capacity]", "r = 0.5", el("code editor")));
    await session.step(45, "And user clicks on Refresh ribbon item", () => clickOn(page, el("Refresh ribbon item")), undefined, "Refresh solves it");
    await session.step(46, "Then the \"Population growth\" view should be current", () => viewIsCurrent(page, "Population growth"));
    await session.step(48, "When user clicks on Edit ribbon item", () => clickOn(page, el("Edit ribbon item")), undefined, "Turn Edit off to get back to the inputs");
    await session.step(49, "Then code editor should be absent", () => shouldBe(page, el("code editor"), "absent"));
    await session.step(50, "And \"r\" input should be visible", () => shouldBe(page, el("\"r\" input"), "visible"));
    await session.step(52, "When user drags the slider of \"r\" input to 1.5", () => dragSliderTo(page, el("\"r\" input"), 1.5), undefined, "Every parameter is now an input with a slider");
    await session.step(53, "Then \"r\" input should have a value between 1.4 and 1.6", () => shouldHaveValueBetween(page, el("\"r\" input"), 1.4, 1.6));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
