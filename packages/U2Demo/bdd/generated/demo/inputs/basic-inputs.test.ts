/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/inputs/basic-inputs.feature", import.meta.url);
  test("Text inputs and a bound search", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(11, "Then Name input should have value \"Aspirin\"", () => shouldHaveValue(page, el("Name input"), "Aspirin"));
    await session.step(12, "And label of Name input should have text \"Name\"", () => shouldHaveText(page, el("label of Name input"), "Name"));
    await session.step(13, "And value of search readout should have text \"(empty)\"", () => shouldHaveText(page, el("value of search readout"), "(empty)"));
    await session.step(14, "When user types \"Kinase\" into Search input", () => typeInto(page, "Kinase", el("Search input")));
    await session.step(15, "Then value of search readout should have text \"Kinase\"", () => shouldHaveText(page, el("value of search readout"), "Kinase"));
    await session.step(16, "And Preview input should contain text \"Kinase\"", () => shouldContainText(page, el("Preview input"), "Kinase"));
    await session.step(17, "When user clears Search input", () => clearField(page, el("Search input")));
    await session.step(18, "Then value of search readout should have text \"(empty)\"", () => shouldHaveText(page, el("value of search readout"), "(empty)"));
  });
  test("A text area", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(21, "When user types \"first line\" into Notes text area", () => typeInto(page, "first line", el("Notes text area")));
    await session.step(22, "Then Notes text area should have value \"first line\"", () => shouldHaveValue(page, el("Notes text area"), "first line"));
  });
  test("Checkboxes and switches", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(25, "Then Enabled checkbox should be checked", () => shouldBe(page, el("Enabled checkbox"), "checked"));
    await session.step(26, "And Notifications checkbox should be checked", () => shouldBe(page, el("Notifications checkbox"), "checked"));
    await session.step(27, "When user unchecks Notifications checkbox", () => uncheck(page, el("Notifications checkbox")));
    await session.step(28, "Then Notifications checkbox should be unchecked", () => shouldBe(page, el("Notifications checkbox"), "unchecked"));
    await session.step(29, "And value of notifications readout should have text \"false\"", () => shouldHaveText(page, el("value of notifications readout"), "false"));
    await session.step(30, "When user toggles Enabled checkbox", () => toggle(page, el("Enabled checkbox")));
    await session.step(31, "Then Enabled checkbox should be unchecked", () => shouldBe(page, el("Enabled checkbox"), "unchecked"));
  });
  test("Numbers feed a computed readout", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(34, "Then value of \"dose * replicates\" readout should have text \"750\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "750"));
    await session.step(35, "When user enters \"4\" into Replicates input", () => enterInto(page, "4", el("Replicates input")));
    await session.step(36, "Then value of \"dose * replicates\" readout should have text \"1000\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "1000"));
    await session.step(37, "When user enters \"100\" into \"Dose, mg\" input", () => enterInto(page, "100", el("\"Dose, mg\" input")));
    await session.step(38, "Then value of \"dose * replicates\" readout should have text \"400\"", () => shouldHaveText(page, el("value of \"dose * replicates\" readout"), "400"));
  });
  test("Single and multiple choices", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(41, "Then Stage input should have value \"Discovery\"", () => shouldHaveValue(page, el("Stage input"), "Discovery"));
    await session.step(42, "When user selects \"Phase II\" in Stage input", () => selectIn(page, "Phase II", el("Stage input")));
    await session.step(43, "Then Stage input should have value \"Phase II\"", () => shouldHaveValue(page, el("Stage input"), "Phase II"));
    await session.step(44, "And Kinase checkbox in Targets input should be checked", () => shouldBe(page, el("Kinase checkbox in Targets input"), "checked"));
    await session.step(45, "When user checks Protease checkbox in Targets input", () => check(page, el("Protease checkbox in Targets input")));
    await session.step(46, "And user unchecks Kinase checkbox in Targets input", () => uncheck(page, el("Kinase checkbox in Targets input")));
    await session.step(47, "Then Protease checkbox in Targets input should be checked", () => shouldBe(page, el("Protease checkbox in Targets input"), "checked"));
    await session.step(48, "And Kinase checkbox in Targets input should be unchecked", () => shouldBe(page, el("Kinase checkbox in Targets input"), "unchecked"));
  });
  test("Validation shows on the input and aggregates into a readout", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(51, "Then error of Code input should have text \"Value is required\"", () => shouldHaveText(page, el("error of Code input"), "Value is required"));
    await session.step(52, "And value of \"code validity\" readout should have text \"Value is required\"", () => shouldHaveText(page, el("value of \"code validity\" readout"), "Value is required"));
    await session.step(53, "When user types \"ABC-123\" into Code input", () => typeInto(page, "ABC-123", el("Code input")));
    await session.step(54, "Then error of Code input should be hidden", () => shouldBe(page, el("error of Code input"), "hidden"));
    await session.step(55, "And value of \"code validity\" readout should have text \"valid\"", () => shouldHaveText(page, el("value of \"code validity\" readout"), "valid"));
    await session.step(56, "When user types \"ABCDEFGHIJKLMNOP\" into Code input", () => typeInto(page, "ABCDEFGHIJKLMNOP", el("Code input")));
    await session.step(57, "Then error of Code input should have text \"At most 10 characters\"", () => shouldHaveText(page, el("error of Code input"), "At most 10 characters"));
  });
  test("Pickers that open a popup", {tag: ["@demo", "@realizes:u2.inputs"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Basic inputs\" demo page", () => openDemoPage(page, "Basic inputs"));
    enter(page, "U2 Demo");
    await session.step(60, "When user selects \"acorn\" in Icon input", () => selectIn(page, "acorn", el("Icon input")));
    await session.step(61, "Then value of icon readout should have text \"acorn\"", () => shouldHaveText(page, el("value of icon readout"), "acorn"));
    await session.step(62, "When user types \"Ky\" into City input", () => typeInto(page, "Ky", el("City input")));
    await session.step(63, "And user selects \"Kyiv\" in City input", () => selectIn(page, "Kyiv", el("City input")));
    await session.step(64, "Then City input should have value \"Kyiv\"", () => shouldHaveValue(page, el("City input"), "Kyiv"));
    await session.step(65, "And value of city readout should have text \"Kyiv\"", () => shouldHaveText(page, el("value of city readout"), "Kyiv"));
    await session.step(66, "When user selects \"Abs\" in Scorer input", () => selectIn(page, "Abs", el("Scorer input")));
    await session.step(67, "Then value of scorer readout should contain text \"Abs\"", () => shouldContainText(page, el("value of scorer readout"), "Abs"));
  });
});
