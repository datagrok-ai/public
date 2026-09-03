/* ---
generated: features/demo/display/sections-and-wizard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.section, u2.wizard]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {clickOn, close, collapse, expand, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sections and the wizard", () => {
  const session = feature(test);
  test("A collapsible section keeps its content", {tag: ["@demo", "@realizes:u2.section", "@realizes:u2.wizard"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Sections & wizard\" demo page", () => openDemoPage(page, "Sections & wizard"));
    enter(page, "U2 Demo");
    await test.step("Then Advanced section should be collapsed", () => shouldBe(page, el("Advanced section"), "collapsed"));
    await test.step("And Threshold input should be hidden", () => shouldBe(page, el("Threshold input"), "hidden"));
    await test.step("When user expands Advanced section", () => expand(page, el("Advanced section")));
    await test.step("Then Threshold input should have value \"0.85\"", () => shouldHaveValue(page, el("Threshold input"), "0.85"));
    await test.step("And value of expanded readout should have text \"true\"", () => shouldHaveText(page, el("value of expanded readout"), "true"));
    await test.step("When user collapses Advanced section", () => collapse(page, el("Advanced section")));
    await test.step("Then Threshold input should be hidden", () => shouldBe(page, el("Threshold input"), "hidden"));
  });
  test("The wizard gates its steps and keeps their state", {tag: ["@demo", "@realizes:u2.section", "@realizes:u2.wizard"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Sections & wizard\" demo page", () => openDemoPage(page, "Sections & wizard"));
    enter(page, "U2 Demo");
    await test.step("Then Name wizard step should be selected", () => shouldBe(page, el("Name wizard step"), "selected"));
    await test.step("And NEXT button in wizard should be disabled", () => shouldBe(page, el("NEXT button in wizard"), "disabled"));
    await test.step("And wizard should contain text \"Enter a project name\"", () => shouldContainText(page, el("wizard"), "Enter a project name"));
    await test.step("When user types \"Atlas\" into Project input in wizard", () => typeInto(page, "Atlas", el("Project input in wizard")));
    await test.step("Then NEXT button in wizard should be enabled", () => shouldBe(page, el("NEXT button in wizard"), "enabled"));
    await test.step("When user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await test.step("Then Data wizard step should be selected", () => shouldBe(page, el("Data wizard step"), "selected"));
    await test.step("And value of step readout should have text \"data\"", () => shouldHaveText(page, el("value of step readout"), "data"));
    await test.step("When user selects \"Protease\" in Targets input in wizard", () => selectIn(page, "Protease", el("Targets input in wizard")));
    await test.step("And user closes Targets input in wizard", () => close(page, el("Targets input in wizard")));
    await test.step("And user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await test.step("Then wizard should contain text \"Creating \\\"Atlas\\\" for Kinase, Protease.\"", () => shouldContainText(page, el("wizard"), "Creating \"Atlas\" for Kinase, Protease."));
    await test.step("When user clicks on BACK button in wizard", () => clickOn(page, el("BACK button in wizard")));
    await test.step("Then Data wizard step should be selected", () => shouldBe(page, el("Data wizard step"), "selected"));
    await test.step("And Targets input in wizard should contain text \"Protease\"", () => shouldContainText(page, el("Targets input in wizard"), "Protease"));
    await test.step("When user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await test.step("And user clicks on FINISH button in wizard", () => clickOn(page, el("FINISH button in wizard")));
    await test.step("Then value of wizard readout should have text \"finished\"", () => shouldHaveText(page, el("value of wizard readout"), "finished"));
  });
});
