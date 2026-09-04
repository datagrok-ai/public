/* ---
generated: features/demo/forms/run-history.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.function-input, u2.func-call-history-browser]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, enterInto, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Run history", () => {
  const session = feature(test);
  test("Running a function and opening its history", {tag: ["@demo", "@realizes:u2.function-input", "@realizes:u2.func-call-history-browser"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Run history\" demo page", () => openDemoPage(page, "Run history"));
    enter(page, "U2 Demo");
    await test.step("When user selects \"Abs\" in Function input", () => selectIn(page, "Abs", el("Function input")));
    await test.step("Then function form should be visible", () => shouldBe(page, el("function form"), "visible"));
    await test.step("And run button of function form should be visible", () => shouldBe(page, el("run button of function form"), "visible"));
    await test.step("When user enters \"3\" into x input in function form", () => enterInto(page, "3", el("x input in function form")));
    await test.step("And user clicks on run button of function form", () => clickOn(page, el("run button of function form")));
    await test.step("And user clicks on history icon of function form", () => clickOn(page, el("history icon of function form")));
    await test.step("Then history browser should be visible", () => shouldBe(page, el("history browser"), "visible"));
  });
});
