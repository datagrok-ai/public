/* ---
generated: features/demo/forms/object-form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.object-form]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {check, clickOn, shouldContainText, shouldHaveText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The object form", () => {
  const session = feature(test);
  test("Saving and deleting a group through the generated form", {tag: ["@demo", "@realizes:u2.object-form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Object form\" demo page", () => openDemoPage(page, "Object form"));
    enter(page, "U2 Demo");
    await test.step("Then Name input in object form should not be empty", () => shouldNotBe(page, el("Name input in object form"), "empty"));
    await test.step("And value of result readout should have text \"(not saved)\"", () => shouldHaveText(page, el("value of result readout"), "(not saved)"));
    await test.step("When user checks Personal checkbox in object form", () => check(page, el("Personal checkbox in object form")));
    await test.step("And user clicks on SAVE button in object form", () => clickOn(page, el("SAVE button in object form")));
    await test.step("Then value of result readout should contain text \"Saved:\"", () => shouldContainText(page, el("value of result readout"), "Saved:"));
    await test.step("When user clicks on Delete button in object form", () => clickOn(page, el("Delete button in object form")));
    await test.step("Then value of result readout should have text \"Deleted\"", () => shouldHaveText(page, el("value of result readout"), "Deleted"));
  });
});
