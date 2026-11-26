/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/11-compute/04-sensitivity-analysis.feature
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
import {clickAndStepOff} from '../../../bindings/howto-compute.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, hoverOver, selectIn, shouldBe, shouldHaveValue, switchOn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFilePutBack} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {hoverArea, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Compute", () => {
  const session = feature(test, "features/howto/11-compute/04-sensitivity-analysis.feature", import.meta.url);
  test("Find which inputs matter with sensitivity analysis", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(20, "And the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the file \"diff-studio-recent.d42\" of the user's home is put back at feature end", () => homeFilePutBack(page, "diff-studio-recent.d42"));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(27, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Diff Studio is in Browse > Apps");
    await session.step(28, "And user clicks on Apps tree node inside browse tree", () => clickOn(page, el("Apps tree node inside browse tree")));
    await session.step(29, "And user double-clicks on Diff-Studio gallery card", () => doubleClickOn(page, el("Diff-Studio gallery card")));
    await session.step(30, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(32, "When user double-clicks on \"PK\" model card", () => doubleClickOn(page, el("\"PK\" model card")), undefined, "A two-compartment PK model from the library");
    await session.step(33, "Then Sensitivity ribbon item should be visible", () => shouldBe(page, el("Sensitivity ribbon item"), "visible"));
    await session.step(34, "When user hovers over browse panel", () => hoverOver(page, el("browse panel")));
    await session.step(36, "And user clicks on browse panel close icon", () => clickOn(page, el("browse panel close icon")), undefined, "Close Browse for more room");
    await session.step(37, "Then browse panel should be hidden", () => shouldBe(page, el("browse panel"), "hidden"));
    await session.step(39, "When user clicks on Sensitivity ribbon item and steps off the ribbon", () => clickAndStepOff(page, el("Sensitivity ribbon item")), undefined, "Sensitivity opens a view of its own");
    await session.step(40, "Then \"Method\" input should be visible", () => shouldBe(page, el("\"Method\" input"), "visible"));
    await session.step(42, "When user selects \"Sobol\" in \"Method\" input", () => selectIn(page, "Sobol", el("\"Method\" input")), undefined, "Sobol: how much of the spread each input explains");
    await session.step(44, "And user enters \"1000\" into \"Samples\" input", () => enterInto(page, "1000", el("\"Samples\" input")), undefined, "Samples: more give steadier answers");
    await session.step(45, "Then \"Samples\" input should have the value \"1000\"", () => shouldHaveValue(page, el("\"Samples\" input"), "1000"));
    await session.step(47, "When user switches on \"Absorption rate\" input", () => switchOn(page, el("\"Absorption rate\" input")), undefined, "Switch on the inputs to vary, each within its range");
    await session.step(48, "And user switches on \"Clearance\" input", () => switchOn(page, el("\"Clearance\" input")));
    await session.step(49, "And user switches on \"Central volume\" input", () => switchOn(page, el("\"Central volume\" input")));
    await session.step(50, "And user switches on \"Peripheral volume\" input", () => switchOn(page, el("\"Peripheral volume\" input")));
    await session.step(51, "And user switches on \"Intercompartmental clearance\" input", () => switchOn(page, el("\"Intercompartmental clearance\" input")));
    await session.step(52, "Then \"Clearance min\" input should be visible", () => shouldBe(page, el("\"Clearance min\" input"), "visible"));
    await session.step(54, "When user hovers over \"PK\" input", () => hoverOver(page, el("\"PK\" input")), undefined, "The output: drug levels at the last time point, 12 h");
    await session.step(55, "Then \"PK\" input should have the value \"Last row\"", () => shouldHaveValue(page, el("\"PK\" input"), "Last row"));
    await session.step(57, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")), undefined, "Run: the model is solved for every sample");
    await session.step(58, "Then the open tableview should have 2 bar chart viewers", () => viewerCount(page, 2, "bar chart"));
    await session.step(59, "When the video says \"The longer the bar, the more it matters: here, clearance\"", () => videoSays(page, "The longer the bar, the more it matters: here, clearance"));
    await session.step(61, "And user hovers over the \"cell Clearance x Central conc., mg/L\" area of correlation plot viewer", () => hoverArea(page, "cell Clearance x Central conc., mg/L", el("correlation plot viewer")), undefined, "Correlations give the direction: more clearance, lower level");
    await session.step(62, "Then no errors should have been logged", () => noErrors(page));
  });
});
