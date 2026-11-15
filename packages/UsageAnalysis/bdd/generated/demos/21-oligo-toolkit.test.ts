/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/21-oligo-toolkit.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openAppOnView, translatorStructureHidden} from '../../bindings/demo-oligo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldContainText, shouldHaveValue, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, dialogCloses, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("3 · Oligonucleotides · Oligo Toolkit", () => {
  const session = feature(test, "features/demos/21-oligo-toolkit.feature", import.meta.url);
  test("One siRNA from notation to pattern", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens the \"Oligo Toolkit\" app on its \"Oligo Translator\" view", () => openAppOnView(page, "Oligo Toolkit", "Oligo Translator"));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And the structure box of the \"Oligo Translator\" view is hidden", () => translatorStructureHidden(page, "Oligo Translator"));
    await session.step(22, "When user types \"csAfsgCfaUfuGfaAfgGfcUfaCfaAfusUfsu\" into text area", () => typeInto(page, "csAfsgCfaUfuGfaAfgGfcUfaCfaAfusUfsu", el("text area")));
    await session.step(23, "Then \"Single sequence format\" choice input should have the value \"Axolabs\"", () => shouldHaveValue(page, el("\"Single sequence format\" choice input"), "Axolabs"));
    await session.step(24, "And \"Nucleotides\" table row should contain the text \"CAGCAUUGAAGGCUACAAUUU\"", () => shouldContainText(page, el("\"Nucleotides\" table row"), "CAGCAUUGAAGGCUACAAUUU"));
    await session.step(25, "When the video says \"One strand, every vendor notation\"", () => videoSays(page, "One strand, every vendor notation"));
    await session.step(26, "And user clicks on PATTERN tab", () => clickOn(page, el("PATTERN tab")));
    await session.step(27, "And user selects \"2'-O-Methyl\" in \"Sequence basis\" choice input", () => selectIn(page, "2'-O-Methyl", el("\"Sequence basis\" choice input")));
    await session.step(28, "And user clicks on \"Edit strands\" button", () => clickOn(page, el("\"Edit strands\" button")));
    await session.step(29, "And user selects \"2'-Fluoro\" in \"Sense modification 2\" choice input in \"Edit strands\" dialog", () => selectIn(page, "2'-Fluoro", el("\"Sense modification 2\" choice input in \"Edit strands\" dialog")));
    await session.step(30, "And user selects \"2'-Fluoro\" in \"Sense modification 4\" choice input in \"Edit strands\" dialog", () => selectIn(page, "2'-Fluoro", el("\"Sense modification 4\" choice input in \"Edit strands\" dialog")));
    await session.step(31, "And user clicks on OK button in \"Edit strands\" dialog", () => clickOn(page, el("OK button in \"Edit strands\" dialog")));
    await session.step(32, "Then the \"Edit strands\" dialog should close", () => dialogCloses(page, "Edit strands"));
    await session.step(33, "When the video says \"2'-fluoro at positions 2 and 4\"", () => videoSays(page, "2'-fluoro at positions 2 and 4"));
    await session.step(34, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
