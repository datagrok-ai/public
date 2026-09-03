/* ---
generated: features/forms-and-tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.form]
--- */
import {test} from '@playwright/test';
import '../bindings/demo.js';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openWorkbench} from '../bindings/steps.js';
import {clickOn, fillIn, followingShouldBe, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms, tables and outlines", () => {
  const session = feature(test);
  test("Filling the form from a table", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user fills in:", () => fillIn(page, [["name input in alignment panel","Batch 7"],["sequence column input in alignment panel","helm"],["method input in alignment panel","clustal"],["gap penalty input in alignment panel","3"],["keep gaps checkbox in alignment panel","no"]]));
    await test.step("Then name input in alignment panel should have value \"Batch 7\"", () => shouldHaveValue(page, el("name input in alignment panel"), "Batch 7"));
    await test.step("And method input in alignment panel should have value \"clustal\"", () => shouldHaveValue(page, el("method input in alignment panel"), "clustal"));
    await test.step("And gap penalty input in alignment panel should have value \"3\"", () => shouldHaveValue(page, el("gap penalty input in alignment panel"), "3"));
    await test.step("And keep gaps checkbox in alignment panel should be unchecked", () => shouldBe(page, el("keep gaps checkbox in alignment panel"), "unchecked"));
    await test.step("And label of name input in alignment panel should have text \"Name\"", () => shouldHaveText(page, el("label of name input in alignment panel"), "Name"));
    await test.step("And the following elements should be visible:", () => followingShouldBe(page, "visible", [["toolbar"],["alignment panel"],["results"]]));
  });
  test("Every method reports itself [method=kalign]", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user selects \"kalign\" in method input in alignment panel", () => selectIn(page, "kalign", el("method input in alignment panel")));
    await test.step("And user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await test.step("And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await test.step("Then status line should contain text \"with kalign\"", () => shouldContainText(page, el("status line"), "with kalign"));
  });
  test("Every method reports itself [method=muscle]", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user selects \"muscle\" in method input in alignment panel", () => selectIn(page, "muscle", el("method input in alignment panel")));
    await test.step("And user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await test.step("And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await test.step("Then status line should contain text \"with muscle\"", () => shouldContainText(page, el("status line"), "with muscle"));
  });
  test("Every method reports itself [method=clustal]", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the MSA workbench", () => openWorkbench(page));
    enter(page, "MSA workbench");
    await test.step("When user selects \"clustal\" in method input in alignment panel", () => selectIn(page, "clustal", el("method input in alignment panel")));
    await test.step("And user clicks on run msa button in alignment panel", () => clickOn(page, el("run msa button in alignment panel")));
    await test.step("And user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await test.step("Then status line should contain text \"with clustal\"", () => shouldContainText(page, el("status line"), "with clustal"));
  });
});
