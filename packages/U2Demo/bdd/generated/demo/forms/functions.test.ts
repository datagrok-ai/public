/* ---
generated: features/demo/forms/functions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.functions-browser, u2.func-form]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, enterInto, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Functions", () => {
  const session = feature(test);
  test("Picking a platform function and running it", {tag: ["@demo", "@realizes:u2.functions-browser", "@realizes:u2.func-form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Functions\" demo page", () => openDemoPage(page, "Functions"));
    enter(page, "U2 Demo");
    await test.step("When user types \"Abs\" into search of functions browser", () => typeInto(page, "Abs", el("search of functions browser")));
    await test.step("And user clicks on Abs item in list of functions browser", () => clickOn(page, el("Abs item in list of functions browser")));
    await test.step("Then Abs heading should be visible", () => shouldBe(page, el("Abs heading"), "visible"));
    await test.step("And function form should be visible", () => shouldBe(page, el("function form"), "visible"));
    await test.step("When user enters \"7.5\" into x input in function form", () => enterInto(page, "7.5", el("x input in function form")));
    await test.step("And user clicks on Run button", () => clickOn(page, el("Run button")));
    await test.step("Then value of result readout should have text \"7.5\"", () => shouldHaveText(page, el("value of result readout"), "7.5"));
  });
});
