/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/forms/form.feature", import.meta.url);
  test("Validity aggregates across the form and Submit reports the values", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user opens the \"Form\" demo page", () => openDemoPage(page, "Form"));
    enter(page, "U2 Demo");
    await session.step(9, "Then value of \"form validity\" readout should have text \"Required\"", () => shouldHaveText(page, el("value of \"form validity\" readout"), "Required"));
    await session.step(10, "And error of \"Last name\" input should have text \"Required\"", () => shouldHaveText(page, el("error of \"Last name\" input"), "Required"));
    await session.step(11, "When user types \"Lovelace\" into \"Last name\" input", () => typeInto(page, "Lovelace", el("\"Last name\" input")));
    await session.step(12, "Then value of \"form validity\" readout should have text \"valid\"", () => shouldHaveText(page, el("value of \"form validity\" readout"), "valid"));
    await session.step(13, "When user clicks on Submit button", () => clickOn(page, el("Submit button")));
    await session.step(14, "Then notification should contain text \"\\\"first\\\":\\\"Ada\\\"\"", () => shouldContainText(page, el("notification"), "\"first\":\"Ada\""));
  });
  test("Filling and resetting", {tag: ["@demo", "@realizes:u2.form"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user opens the \"Form\" demo page", () => openDemoPage(page, "Form"));
    enter(page, "U2 Demo");
    await session.step(17, "When user fills in:", () => fillIn(page, [["\"First name\" input","Grace"],["\"Last name\" input","Hopper"],["Age input","85"],["Role input","Admin"],["Subscribe checkbox","no"]]));
    await session.step(23, "Then Age input should have value \"85\"", () => shouldHaveValue(page, el("Age input"), "85"));
    await session.step(24, "And Role input should have value \"Admin\"", () => shouldHaveValue(page, el("Role input"), "Admin"));
    await session.step(25, "And Subscribe checkbox should be unchecked", () => shouldBe(page, el("Subscribe checkbox"), "unchecked"));
    await session.step(26, "When user clicks on Reset button", () => clickOn(page, el("Reset button")));
    await session.step(27, "Then \"First name\" input should have value \"Ada\"", () => shouldHaveValue(page, el("\"First name\" input"), "Ada"));
    await session.step(28, "And \"Last name\" input should be empty", () => shouldBe(page, el("\"Last name\" input"), "empty"));
    await session.step(29, "And Role input should have value \"Editor\"", () => shouldHaveValue(page, el("Role input"), "Editor"));
    await session.step(30, "And Subscribe checkbox should be checked", () => shouldBe(page, el("Subscribe checkbox"), "checked"));
  });
});
