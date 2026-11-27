/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/03-fit-a-model-to-data.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
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
import {clickAndStepOff, switchBeside} from '../../../bindings/howto-compute.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, expand, hoverOver, scrollIntoView, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, currentViewShows, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/03-fit-a-model-to-data.feature", import.meta.url);
  test("Fit a model to data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And simple mode is off", () => simpleModeOff(page));
    await session.step(23, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(24, "And the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(27, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(28, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(29, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(31, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Diff Studio is in Browse > Apps");
    await session.step(32, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(33, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(34, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(36, "When user double-clicks on \"PK\" model card", () => doubleClickOn(page, el("\"PK\" model card")), undefined, "A two-compartment PK model from the library");
    await session.step(37, "Then Fit ribbon item should be visible", () => shouldBe(page, el("Fit ribbon item"), "visible"));
    await session.step(38, "When user hovers over browse panel", () => hoverOver(page, el("browse panel")));
    await session.step(40, "And user clicks on browse panel close icon", () => clickOn(page, el("browse panel close icon")), undefined, "Close Browse for more room");
    await session.step(41, "Then browse panel should be hidden", () => shouldBe(page, el("browse panel"), "hidden"));
    await session.step(43, "When user clicks on Fit ribbon item and steps off the ribbon", () => clickAndStepOff(page, el("Fit ribbon item")), undefined, "Fit opens a view of its own");
    await session.step(44, "Then the current view should show the text \"Target\"", () => currentViewShows(page, "Target"));
    await session.step(45, "When the video says \"Target: what the model should match\"", () => videoSays(page, "Target: what the model should match"));
    await session.step(47, "And user clicks on \"Add file from Files\" icon", () => clickOn(page, el("\"Add file from Files\" icon")), undefined, "Here, a table of drug levels over time, from Files");
    await session.step(48, "Then \"Select a file\" dialog should be visible", () => shouldBe(page, el("\"Select a file\" dialog"), "visible"));
    await session.step(49, "When user expands \"Files > App Data\" tree node inside \"Select a file\" dialog", () => expand(page, el("\"Files > App Data\" tree node inside \"Select a file\" dialog")));
    await session.step(50, "And user expands \"Files > App Data > DiffStudio\" tree node inside \"Select a file\" dialog", () => expand(page, el("\"Files > App Data > DiffStudio\" tree node inside \"Select a file\" dialog")));
    await session.step(51, "And user expands \"Files > App Data > DiffStudio > library\" tree node inside \"Select a file\" dialog", () => expand(page, el("\"Files > App Data > DiffStudio > library\" tree node inside \"Select a file\" dialog")));
    await session.step(52, "And user scrolls \"Files > App Data > DiffStudio > library > pk.csv\" tree node inside \"Select a file\" dialog into view", () => scrollIntoView(page, el("\"Files > App Data > DiffStudio > library > pk.csv\" tree node inside \"Select a file\" dialog")));
    await session.step(53, "And user double-clicks on \"Files > App Data > DiffStudio > library > pk.csv\" tree node inside \"Select a file\" dialog", () => doubleClickOn(page, el("\"Files > App Data > DiffStudio > library > pk.csv\" tree node inside \"Select a file\" dialog")));
    await session.step(54, "Then \"PK\" input should have the value \"pk\"", () => shouldHaveValue(page, el("\"PK\" input"), "pk"));
    await session.step(55, "And \"argument\" input should have the value \"Time, h\"", () => shouldHaveValue(page, el("\"argument\" input"), "Time, h"));
    await session.step(57, "When user enters \"3\" into \"Count\" input", () => enterInto(page, "3", el("\"Count\" input")), undefined, "The data covers three doses: so must the model");
    await session.step(59, "And user switches on the switch beside \"Absorption rate\" input", () => switchBeside(page, el("\"Absorption rate\" input")), undefined, "Switch on the parameters to fit; each gets a range");
    await session.step(60, "And user switches on the switch beside \"Clearance\" input", () => switchBeside(page, el("\"Clearance\" input")));
    await session.step(61, "And user switches on the switch beside \"Central volume\" input", () => switchBeside(page, el("\"Central volume\" input")));
    await session.step(62, "And user switches on the switch beside \"Peripheral volume\" input", () => switchBeside(page, el("\"Peripheral volume\" input")));
    await session.step(63, "And user switches on the switch beside \"Intercompartmental clearance\" input", () => switchBeside(page, el("\"Intercompartmental clearance\" input")));
    await session.step(64, "Then \"Clearance (min)\" input should be visible", () => shouldBe(page, el("\"Clearance (min)\" input"), "visible"));
    await session.step(66, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")), undefined, "Run the fit");
    await session.step(67, "Then the table should have a column \"RMSE by iterations\"", () => hasColumn(page, "RMSE by iterations"));
    await session.step(68, "When the video says \"The best fit comes first: its values, and the model over the data\"", () => videoSays(page, "The best fit comes first: its values, and the model over the data"));
    await session.step(69, "Then no errors should have been logged", () => noErrors(page));
  });
});
