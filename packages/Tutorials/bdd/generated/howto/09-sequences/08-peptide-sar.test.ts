/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/08-peptide-sar.feature
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
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {peptidesReady, sarAnalysisReady} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelClosed, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/08-peptide-sar.feature", import.meta.url);
  test("Analyze peptide SAR", {tag: ["@guide", "@tutorial", "@viewport:1920x1080"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the Peptides package is initialized", () => peptidesReady(page));
    await session.step(19, "And user opens peptides dataset keeping the first 200 rows", () => openDatasetRows(page, ds("peptides"), 200));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "When the video says \"200 aligned peptides with their IC50\"", () => videoSays(page, "200 aligned peptides with their IC50"));
    await session.step(24, "And user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."), undefined, "Bio > Analyze > SAR");
    await session.step(25, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(26, "And editor of Sequence input in \"Analyze Peptides\" dialog should have text \"AlignedSequence\"", () => shouldHaveText(page, el("editor of Sequence input in \"Analyze Peptides\" dialog"), "AlignedSequence"));
    await session.step(27, "And editor of Activity input in \"Analyze Peptides\" dialog should have text \"IC50\"", () => shouldHaveText(page, el("editor of Activity input in \"Analyze Peptides\" dialog"), "IC50"));
    await session.step(29, "When user selects \"-lg\" in Scaling input in \"Analyze Peptides\" dialog", () => selectIn(page, "-lg", el("Scaling input in \"Analyze Peptides\" dialog")), undefined, "Sequence and activity are picked; -lg turns IC50 into pIC50");
    await session.step(31, "And user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")), undefined, "Clusters: how alike peptides must be to group");
    await session.step(32, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(33, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(34, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(35, "Then the SAR analysis should be ready", () => sarAnalysisReady(page));
    await session.step(36, "When the video says \"Bottom left: every residue at every position. Top right: the clusters\"", () => videoSays(page, "Bottom left: every residue at every position. Top right: the clusters"));
    await session.step(38, "And user hovers over the \"cell A at 2\" area of Sequence Variability Map viewer", () => hoverArea(page, "cell A at 2", el("Sequence Variability Map viewer")), undefined, "Circle size: peptides in cliffs; color: how activity changes");
    await session.step(39, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(41, "When user clicks on the \"cell A at 2\" area of Sequence Variability Map viewer", () => clickArea(page, "cell A at 2", el("Sequence Variability Map viewer")), undefined, "Click a residue: its peptides are selected and listed");
    await session.step(42, "Then some rows should be selected", () => someSelected(page));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
