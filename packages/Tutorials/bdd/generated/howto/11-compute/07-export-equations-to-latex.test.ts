/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/07-export-equations-to-latex.feature
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
import {clickOn, clipboardContains, doubleClickOn, shouldBe, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, dialogCloses, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/07-export-equations-to-latex.feature", import.meta.url);
  test("Export a model's equations to LaTeX", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Diff Studio is in Browse > Apps");
    await session.step(24, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(25, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(26, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(27, "When user double-clicks on \"Lotka-Volterra\" model card", () => doubleClickOn(page, el("\"Lotka-Volterra\" model card")));
    await session.step(28, "Then \"Prey\" input should be visible", () => shouldBe(page, el("\"Prey\" input"), "visible"));
    await session.step(30, "When user clicks on \"Save or export model\" icon", () => clickOn(page, el("\"Save or export model\" icon")), undefined, "This icon saves or exports the model");
    await session.step(31, "And user picks \"As LaTeX...\" from the open menu", () => pickFromOpenMenu(page, "As LaTeX..."));
    await session.step(32, "Then \"Export as LaTeX\" dialog should be visible", () => shouldBe(page, el("\"Export as LaTeX\" dialog"), "visible"));
    await session.step(33, "When the video says \"The preview is the LaTeX source, ready to use\"", () => videoSays(page, "The preview is the LaTeX source, ready to use"));
    await session.step(35, "And user unchecks \"Standalone\" input in \"Export as LaTeX\" dialog", () => uncheck(page, el("\"Standalone\" input in \"Export as LaTeX\" dialog")), undefined, "Standalone off: only the part to paste into your document");
    await session.step(37, "And user clicks on \"Copy the generated source to the clipboard.\" icon in \"Export as LaTeX\" dialog", () => clickOn(page, el("\"Copy the generated source to the clipboard.\" icon in \"Export as LaTeX\" dialog")), undefined, "Copy it, or download it as a .tex file");
    await session.step(38, "Then the clipboard should contain the text \"alpha\"", () => clipboardContains(page, "alpha"));
    await session.step(39, "When user clicks on CANCEL button in \"Export as LaTeX\" dialog", () => clickOn(page, el("CANCEL button in \"Export as LaTeX\" dialog")));
    await session.step(40, "Then the \"Export as LaTeX\" dialog should close", () => dialogCloses(page, "Export as LaTeX"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
});
