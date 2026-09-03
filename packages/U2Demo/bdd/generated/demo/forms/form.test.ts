/* ---
generated: features/demo/forms/form.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.form]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, fillIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms", () => {
  const session = feature(test);
  test("Validity aggregates across the form and Submit reports the values", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Form\" demo page", () => openDemoPage(page, "Form"));
    enter(page, "U2 Demo");
    await test.step("Then value of \"form validity\" readout should have text \"Required\"", () => shouldHaveText(page, el("value of \"form validity\" readout"), "Required"));
    await test.step("And error of \"Last name\" input should have text \"Required\"", () => shouldHaveText(page, el("error of \"Last name\" input"), "Required"));
    await test.step("When user types \"Lovelace\" into \"Last name\" input", () => typeInto(page, "Lovelace", el("\"Last name\" input")));
    await test.step("Then value of \"form validity\" readout should have text \"valid\"", () => shouldHaveText(page, el("value of \"form validity\" readout"), "valid"));
    await test.step("When user clicks on Submit button", () => clickOn(page, el("Submit button")));
    await test.step("Then notification should contain text \"\\\"first\\\":\\\"Ada\\\"\"", () => shouldContainText(page, el("notification"), "\"first\":\"Ada\""));
  });
  test("Filling and resetting", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Form\" demo page", () => openDemoPage(page, "Form"));
    enter(page, "U2 Demo");
    await test.step("When user fills in:", () => fillIn(page, [["\"First name\" input","Grace"],["\"Last name\" input","Hopper"],["Age input","85"],["Role input","Admin"],["Subscribe checkbox","no"]]));
    await test.step("Then Age input should have value \"85\"", () => shouldHaveValue(page, el("Age input"), "85"));
    await test.step("And Role input should have value \"Admin\"", () => shouldHaveValue(page, el("Role input"), "Admin"));
    await test.step("And Subscribe checkbox should be unchecked", () => shouldBe(page, el("Subscribe checkbox"), "unchecked"));
    await test.step("When user clicks on Reset button", () => clickOn(page, el("Reset button")));
    await test.step("Then \"First name\" input should have value \"Ada\"", () => shouldHaveValue(page, el("\"First name\" input"), "Ada"));
    await test.step("And \"Last name\" input should be empty", () => shouldBe(page, el("\"Last name\" input"), "empty"));
    await test.step("And Role input should have value \"Editor\"", () => shouldHaveValue(page, el("Role input"), "Editor"));
    await test.step("And Subscribe checkbox should be checked", () => shouldBe(page, el("Subscribe checkbox"), "checked"));
  });
});
