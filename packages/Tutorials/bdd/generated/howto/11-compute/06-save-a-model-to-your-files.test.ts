/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/06-save-a-model-to-your-files.feature
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
import {clickOn, doubleClickOn, enterInto, expand, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, dialogCloses, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/06-save-a-model-to-your-files.feature", import.meta.url);
  test("Save a model to your files", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(21, "And the file \"Predator-prey.ivp\" of the user's home is put back at feature end", () => homeFilePutBack(page, "Predator-prey.ivp"));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Diff Studio is in Browse > Apps");
    await session.step(27, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(28, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(29, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(30, "When user double-clicks on \"Lotka-Volterra\" model card", () => doubleClickOn(page, el("\"Lotka-Volterra\" model card")));
    await session.step(31, "Then \"Prey\" input should be visible", () => shouldBe(page, el("\"Prey\" input"), "visible"));
    await session.step(33, "When user clicks on Save button", () => clickOn(page, el("Save button")), undefined, "SAVE keeps a copy of the model in your files");
    await session.step(34, "Then \"Save to My files\" dialog should be visible", () => shouldBe(page, el("\"Save to My files\" dialog"), "visible"));
    await session.step(36, "When user enters \"Predator-prey\" into \"Name\" input in \"Save to My files\" dialog", () => enterInto(page, "Predator-prey", el("\"Name\" input in \"Save to My files\" dialog")), undefined, "Give it a name");
    await session.step(37, "And user clicks on Save button in \"Save to My files\" dialog", () => clickOn(page, el("Save button in \"Save to My files\" dialog")));
    await session.step(38, "Then the \"Save to My files\" dialog should close", () => dialogCloses(page, "Save to My files"));
    await session.step(40, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")), undefined, "It is in Browse > Files > My files");
    await session.step(41, "And user clicks on \"Files > My files\" tree node inside browse tree", () => clickOn(page, el("\"Files > My files\" tree node inside browse tree")));
    await session.step(42, "Then Predator-prey.ivp link in gallery should be visible", () => shouldBe(page, el("Predator-prey.ivp link in gallery"), "visible"));
    await session.step(44, "When user clicks on Predator-prey.ivp link in gallery", () => clickOn(page, el("Predator-prey.ivp link in gallery")), undefined, "Click it, and the model opens with its inputs");
    await session.step(45, "Then \"Prey\" input should be visible", () => shouldBe(page, el("\"Prey\" input"), "visible"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
