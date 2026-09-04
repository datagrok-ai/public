/* ---
generated: features/demo/inputs/pickers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.range-slider, u2.multi-select, u2.button-group, u2.combobox]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, close, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotBe, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sliders, multi-select, button groups and comboboxes", () => {
  const session = feature(test);
  test("Range slider handles answer the keyboard", {tag: ["@demo", "@realizes:u2.range-slider", "@realizes:u2.multi-select", "@realizes:u2.button-group", "@realizes:u2.combobox"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Range slider\" demo page", () => openDemoPage(page, "Range slider"));
    enter(page, "U2 Demo");
    await test.step("Then value of range readout should have text \"20 … 70\"", () => shouldHaveText(page, el("value of range readout"), "20 … 70"));
    await test.step("When user presses ArrowRight in Minimum slider handle in first range slider", () => pressKeyIn(page, "ArrowRight", el("Minimum slider handle in first range slider")));
    await test.step("Then value of range readout should have text \"21 … 70\"", () => shouldHaveText(page, el("value of range readout"), "21 … 70"));
    await test.step("When user presses ArrowLeft in Maximum slider handle in first range slider", () => pressKeyIn(page, "ArrowLeft", el("Maximum slider handle in first range slider")));
    await test.step("Then value of range readout should have text \"21 … 69\"", () => shouldHaveText(page, el("value of range readout"), "21 … 69"));
    await test.step("When user presses ArrowRight in Minimum slider handle in second range slider", () => pressKeyIn(page, "ArrowRight", el("Minimum slider handle in second range slider")));
    await test.step("Then value of dose readout should have text \"10 … 1000 mg\"", () => shouldHaveText(page, el("value of dose readout"), "10 … 1000 mg"));
  });
  test("A multi-select shows its picks as tags", {tag: ["@demo", "@realizes:u2.range-slider", "@realizes:u2.multi-select", "@realizes:u2.button-group", "@realizes:u2.combobox"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Multi-select\" demo page", () => openDemoPage(page, "Multi-select"));
    enter(page, "U2 Demo");
    await test.step("Then value of targets readout should have text \"Kinase\"", () => shouldHaveText(page, el("value of targets readout"), "Kinase"));
    await test.step("When user selects \"Protease\" in Targets input", () => selectIn(page, "Protease", el("Targets input")));
    await test.step("And user closes Targets input", () => close(page, el("Targets input")));
    await test.step("Then value of targets readout should have text \"Kinase, Protease\"", () => shouldHaveText(page, el("value of targets readout"), "Kinase, Protease"));
    await test.step("When user clicks on \"Remove Kinase\" button in Targets input", () => clickOn(page, el("\"Remove Kinase\" button in Targets input")));
    await test.step("Then value of targets readout should have text \"Protease\"", () => shouldHaveText(page, el("value of targets readout"), "Protease"));
  });
  test("Button groups as actions, a single toggle and a multi toggle", {tag: ["@demo", "@realizes:u2.range-slider", "@realizes:u2.multi-select", "@realizes:u2.button-group", "@realizes:u2.combobox"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Multi-select\" demo page", () => openDemoPage(page, "Multi-select"));
    enter(page, "U2 Demo");
    await test.step("Then Grid button should be selected", () => shouldBe(page, el("Grid button"), "selected"));
    await test.step("When user clicks on Cards button", () => clickOn(page, el("Cards button")));
    await test.step("Then value of layout readout should have text \"cards\"", () => shouldHaveText(page, el("value of layout readout"), "cards"));
    await test.step("And Cards button should be selected", () => shouldBe(page, el("Cards button"), "selected"));
    await test.step("And Grid button should not be selected", () => shouldNotBe(page, el("Grid button"), "selected"));
    await test.step("When user clicks on Bold button", () => clickOn(page, el("Bold button")));
    await test.step("And user clicks on Italic button", () => clickOn(page, el("Italic button")));
    await test.step("Then value of style readout should have text \"bold, italic\"", () => shouldHaveText(page, el("value of style readout"), "bold, italic"));
    await test.step("And Bold button should be selected", () => shouldBe(page, el("Bold button"), "selected"));
    await test.step("When user clicks on Bold button", () => clickOn(page, el("Bold button")));
    await test.step("Then value of style readout should have text \"italic\"", () => shouldHaveText(page, el("value of style readout"), "italic"));
  });
  test("Comboboxes over local and async items", {tag: ["@demo", "@realizes:u2.range-slider", "@realizes:u2.multi-select", "@realizes:u2.button-group", "@realizes:u2.combobox"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Async\" demo page", () => openDemoPage(page, "Async"));
    enter(page, "U2 Demo");
    await test.step("When user selects \"Kyiv\" in first combobox", () => selectIn(page, "Kyiv", el("first combobox")));
    await test.step("Then value of local readout should have text \"Kyiv\"", () => shouldHaveText(page, el("value of local readout"), "Kyiv"));
    await test.step("When user types \"Osl\" into second combobox", () => typeInto(page, "Osl", el("second combobox")));
    await test.step("And user selects \"Oslo\" in second combobox", () => selectIn(page, "Oslo", el("second combobox")));
    await test.step("Then value of remote readout should have text \"Oslo\"", () => shouldHaveText(page, el("value of remote readout"), "Oslo"));
  });
  test("An async view refreshes on its filter", {tag: ["@demo", "@realizes:u2.range-slider", "@realizes:u2.multi-select", "@realizes:u2.button-group", "@realizes:u2.combobox"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Async\" demo page", () => openDemoPage(page, "Async"));
    enter(page, "U2 Demo");
    await test.step("Then async view should contain text \"Berlin\"", () => shouldContainText(page, el("async view"), "Berlin"));
    await test.step("When user types \"Par\" into \"Filter the async view…\" input", () => typeInto(page, "Par", el("\"Filter the async view…\" input")));
    await test.step("Then async view should contain text \"Paris\"", () => shouldContainText(page, el("async view"), "Paris"));
    await test.step("And async view should not contain text \"Berlin\"", () => shouldNotContainText(page, el("async view"), "Berlin"));
  });
});
