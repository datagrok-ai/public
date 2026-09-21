/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/msa-workbench.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dialog, u2.tabs]
--- */
import {test} from '@playwright/test';
import '../bindings/demo.js';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openWorkbench} from '../bindings/steps.js';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveItems, shouldHaveText, shouldHaveValue, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MSA workbench", () => {
  const session = feature(test, "features/msa-workbench.feature", import.meta.url);
  test("Running an alignment from the dialog", {tag: ["@demo", "@realizes:u2.dialog", "@realizes:u2.tabs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(11, "When user types \"Kinase panel\" into name input in alignment panel", () => typeInto(page, "Kinase panel", el("name input in alignment panel")));
    await session.step(12, "And user selects \"peptide\" in sequence column input in alignment panel", () => selectIn(page, "peptide", el("sequence column input in alignment panel")));
    await session.step(13, "And user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await session.step(14, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(15, "And sequence column input in MSA dialog should have value \"peptide\"", () => shouldHaveValue(page, el("sequence column input in MSA dialog"), "peptide"));
    await session.step(16, "When user selects \"muscle\" in method input in MSA dialog", () => selectIn(page, "muscle", el("method input in MSA dialog")));
    await session.step(17, "And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await session.step(18, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(19, "And results should be visible", () => shouldBe(page, el("results"), "visible"));
    await session.step(20, "And status line should contain text \"Aligned 5 sequences with muscle\"", () => shouldContainText(page, el("status line"), "Aligned 5 sequences with muscle"));
    await session.step(21, "And aligned sequences list should have 5 items", () => shouldHaveItems(page, el("aligned sequences list"), 5));
    await session.step(22, "And notification should contain text \"Alignment finished\"", () => shouldContainText(page, el("notification"), "Alignment finished"));
  });
  test("Cancelling the dialog keeps the results untouched", {tag: ["@demo", "@realizes:u2.dialog", "@realizes:u2.tabs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(25, "When user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await session.step(26, "And user clicks on cancel button in MSA dialog", () => clickOn(page, el("cancel button in MSA dialog")));
    await session.step(27, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(28, "And status line should have text \"No alignment yet\"", () => shouldHaveText(page, el("status line"), "No alignment yet"));
    await session.step(29, "And aligned sequences list should have 0 items", () => shouldHaveItems(page, el("aligned sequences list"), 0));
  });
  test("Switching result tabs", {tag: ["@demo", "@realizes:u2.dialog", "@realizes:u2.tabs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await session.step(32, "When user selects \"muscle\" in method input in alignment panel", () => selectIn(page, "muscle", el("method input in alignment panel")));
    await session.step(33, "And user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await session.step(34, "And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await session.step(35, "And user clicks on Log tab in results", () => clickOn(page, el("Log tab in results")));
    await session.step(36, "Then Log tab in results should be selected", () => shouldBe(page, el("Log tab in results"), "selected"));
    await session.step(37, "And log should contain text \"muscle\"", () => shouldContainText(page, el("log"), "muscle"));
    await session.step(38, "And Alignment tab in results should not be selected", () => shouldNotBe(page, el("Alignment tab in results"), "selected"));
  });
});
