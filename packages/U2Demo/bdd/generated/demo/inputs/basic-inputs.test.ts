/* ---
generated: features/demo/inputs/basic-inputs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.inputs]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {check, clearField, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, toggle, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Basic inputs", () => {
  const session = feature(test);
  test("Text inputs and a bound search", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("Then Name input should have value \"Aspirin\"", () => shouldHaveValue(page, el("Name input"), "Aspirin"));
    await test.step("And label of Name input should have text \"Name\"", () => shouldHaveText(page, el("label of Name input"), "Name"));
    await test.step("And value of search readout should have text \"(empty)\"", () => shouldHaveText(page, el("value of search readout"), "(empty)"));
    await test.step("When user types \"Kinase\" into Search input", () => typeInto(page, "Kinase", el("Search input")));
    await test.step("Then value of search readout should have text \"Kinase\"", () => shouldHaveText(page, el("value of search readout"), "Kinase"));
    await test.step("And Preview input should contain text \"Kinase\"", () => shouldContainText(page, el("Preview input"), "Kinase"));
    await test.step("When user clears Search input", () => clearField(page, el("Search input")));
    await test.step("Then value of search readout should have text \"(empty)\"", () => shouldHaveText(page, el("value of search readout"), "(empty)"));
  });
  test("A text area", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("When user types \"first line\" into Notes text area", () => typeInto(page, "first line", el("Notes text area")));
    await test.step("Then Notes text area should have value \"first line\"", () => shouldHaveValue(page, el("Notes text area"), "first line"));
  });
  test("Checkboxes and switches", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("Then Enabled checkbox should be checked", () => shouldBe(page, el("Enabled checkbox"), "checked"));
    await test.step("And Notifications checkbox should be checked", () => shouldBe(page, el("Notifications checkbox"), "checked"));
    await test.step("When user unchecks Notifications checkbox", () => uncheck(page, el("Notifications checkbox")));
    await test.step("Then Notifications checkbox should be unchecked", () => shouldBe(page, el("Notifications checkbox"), "unchecked"));
    await test.step("And value of notifications readout should have text \"false\"", () => shouldHaveText(page, el("value of notifications readout"), "false"));
    await test.step("When user toggles Enabled checkbox", () => toggle(page, el("Enabled checkbox")));
    await test.step("Then Enabled checkbox should be unchecked", () => shouldBe(page, el("Enabled checkbox"), "unchecked"));
  });
  test("Numbers feed a computed readout", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("Then value of \"dose * replicates\" readout should have text \"750\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "750"));
    await test.step("When user enters \"4\" into Replicates input", () => enterInto(page, "4", el("Replicates input")));
    await test.step("Then value of \"dose * replicates\" readout should have text \"1000\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "1000"));
    await test.step("When user enters \"100\" into \"Dose, mg\" input", () => enterInto(page, "100", el("\"Dose, mg\" input")));
    await test.step("Then value of \"dose * replicates\" readout should have text \"400\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "400"));
  });
  test("Single and multiple choices", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("Then Stage input should have value \"Discovery\"", () => shouldHaveValue(page, el("Stage input"), "Discovery"));
    await test.step("When user selects \"Phase II\" in Stage input", () => selectIn(page, "Phase II", el("Stage input")));
    await test.step("Then Stage input should have value \"Phase II\"", () => shouldHaveValue(page, el("Stage input"), "Phase II"));
    await test.step("And Kinase checkbox in Targets input should be checked", () => shouldBe(page, el("Kinase checkbox in Targets input"), "checked"));
    await test.step("When user checks Protease checkbox in Targets input", () => check(page, el("Protease checkbox in Targets input")));
    await test.step("And user unchecks Kinase checkbox in Targets input", () => uncheck(page, el("Kinase checkbox in Targets input")));
    await test.step("Then Protease checkbox in Targets input should be checked", () => shouldBe(page, el("Protease checkbox in Targets input"), "checked"));
    await test.step("And Kinase checkbox in Targets input should be unchecked", () => shouldBe(page, el("Kinase checkbox in Targets input"), "unchecked"));
  });
  test("Validation shows on the input and aggregates into a readout", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("Then error of Code input should have text \"Value is required\"", () => shouldHaveText(page, el("error of Code input"), "Value is required"));
    await test.step("And value of \"code validity\" readout should have text \"Value is required\"", () => shouldHaveText(page, el("value of \"code validity\" readout"), "Value is required"));
    await test.step("When user types \"ABC-123\" into Code input", () => typeInto(page, "ABC-123", el("Code input")));
    await test.step("Then error of Code input should be hidden", () => shouldBe(page, el("error of Code input"), "hidden"));
    await test.step("And value of \"code validity\" readout should have text \"valid\"", () => shouldHaveText(page, el("value of \"code validity\" readout"), "valid"));
    await test.step("When user types \"ABCDEFGHIJKLMNOP\" into Code input", () => typeInto(page, "ABCDEFGHIJKLMNOP", el("Code input")));
    await test.step("Then error of Code input should have text \"At most 10 characters\"", () => shouldHaveText(page, el("error of Code input"), "At most 10 characters"));
  });
  test("Pickers that open a popup", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await test.step("When user selects \"acorn\" in Icon input", () => selectIn(page, "acorn", el("Icon input")));
    await test.step("Then value of icon readout should have text \"acorn\"", () => shouldHaveText(page, el("value of icon readout"), "acorn"));
    await test.step("When user types \"Ky\" into City input", () => typeInto(page, "Ky", el("City input")));
    await test.step("And user selects \"Kyiv\" in City input", () => selectIn(page, "Kyiv", el("City input")));
    await test.step("Then City input should have value \"Kyiv\"", () => shouldHaveValue(page, el("City input"), "Kyiv"));
    await test.step("And value of city readout should have text \"Kyiv\"", () => shouldHaveText(page, el("value of city readout"), "Kyiv"));
    await test.step("When user selects \"Abs\" in Scorer input", () => selectIn(page, "Abs", el("Scorer input")));
    await test.step("Then value of scorer readout should contain text \"Abs\"", () => shouldContainText(page, el("value of scorer readout"), "Abs"));
  });
});
