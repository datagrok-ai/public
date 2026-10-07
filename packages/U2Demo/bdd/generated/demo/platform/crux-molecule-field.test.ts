/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/crux-molecule-field.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, packageInstalled, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A u2 form's molecule field, edited in Crux", () => {
  const session = feature(test, "features/demo/platform/crux-molecule-field.feature", import.meta.url);
  test("Tab to the form's Structure field and Enter opens Crux, and the drawn molecule writes back through the bound property", {tag: ["@demo", "@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(16, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "Given user opens the \"Molecules\" demo page", () => openDemoPage(page, "Molecules"));
    enter(page, "U2 Demo");
    await session.step(22, "When user clicks on Name input in object form", () => clickOn(page, el("Name input in object form")));
    await session.step(23, "And user presses Tab", () => pressKey(page, "Tab"));
    await session.step(24, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(25, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(26, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"CC(=O)OC1=CC=CC=C1C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CC(=O)OC1=CC=CC=C1C(=O)O"));
    await session.step(27, "When user clicks on crux clear button", () => clickOn(page, el("crux clear button")));
    await session.step(28, "And user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(29, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(30, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccccc1"));
    await session.step(31, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(32, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(33, "And value of compound readout should contain text \"\\\"smiles\\\":\\\"c1ccccc1\\\"\"", () => shouldContainText(page, el("value of compound readout"), "\"smiles\":\"c1ccccc1\""));
    await session.step(34, "And value of compound readout should contain text \"\\\"name\\\":\\\"Aspirin\\\"\"", () => shouldContainText(page, el("value of compound readout"), "\"name\":\"Aspirin\""));
  });
});
