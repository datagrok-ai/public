/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/display/sections-and-wizard.feature", import.meta.url);
  test("A collapsible section keeps its content", {tag: ["@demo", "@realizes:u2.section", "@realizes:u2.wizard"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Sections & wizard\" demo page", () => openDemoPage(page, "Sections & wizard"));
    enter(page, "U2 Demo");
    await session.step(8, "Then Advanced section should be collapsed", () => shouldBe(page, el("Advanced section"), "collapsed"));
    await session.step(9, "And Threshold input should be hidden", () => shouldBe(page, el("Threshold input"), "hidden"));
    await session.step(10, "When user expands Advanced section", () => expand(page, el("Advanced section")));
    await session.step(11, "Then Threshold input should have value \"0.85\"", () => shouldHaveValue(page, el("Threshold input"), "0.85"));
    await session.step(12, "And value of expanded readout should have text \"true\"", () => shouldHaveText(page, el("value of expanded readout"), "true"));
    await session.step(13, "When user collapses Advanced section", () => collapse(page, el("Advanced section")));
    await session.step(14, "Then Threshold input should be hidden", () => shouldBe(page, el("Threshold input"), "hidden"));
  });
  test("The wizard gates its steps and keeps their state", {tag: ["@demo", "@realizes:u2.section", "@realizes:u2.wizard"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Sections & wizard\" demo page", () => openDemoPage(page, "Sections & wizard"));
    enter(page, "U2 Demo");
    await session.step(17, "Then Name wizard step should be selected", () => shouldBe(page, el("Name wizard step"), "selected"));
    await session.step(18, "And NEXT button in wizard should be disabled", () => shouldBe(page, el("NEXT button in wizard"), "disabled"));
    await session.step(19, "And wizard should contain text \"Enter a project name\"", () => shouldContainText(page, el("wizard"), "Enter a project name"));
    await session.step(20, "When user types \"Atlas\" into Project input in wizard", () => typeInto(page, "Atlas", el("Project input in wizard")));
    await session.step(21, "Then NEXT button in wizard should be enabled", () => shouldBe(page, el("NEXT button in wizard"), "enabled"));
    await session.step(22, "When user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await session.step(23, "Then Data wizard step should be selected", () => shouldBe(page, el("Data wizard step"), "selected"));
    await session.step(24, "And value of step readout should have text \"data\"", () => shouldHaveText(page, el("value of step readout"), "data"));
    await session.step(25, "When user selects \"Protease\" in Targets input in wizard", () => selectIn(page, "Protease", el("Targets input in wizard")));
    await session.step(26, "And user closes Targets input in wizard", () => close(page, el("Targets input in wizard")));
    await session.step(27, "And user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await session.step(28, "Then wizard should contain text \"Creating \\\"Atlas\\\" for Kinase, Protease.\"", () => shouldContainText(page, el("wizard"), "Creating \"Atlas\" for Kinase, Protease."));
    await session.step(29, "When user clicks on BACK button in wizard", () => clickOn(page, el("BACK button in wizard")));
    await session.step(30, "Then Data wizard step should be selected", () => shouldBe(page, el("Data wizard step"), "selected"));
    await session.step(31, "And Targets input in wizard should contain text \"Protease\"", () => shouldContainText(page, el("Targets input in wizard"), "Protease"));
    await session.step(32, "When user clicks on NEXT button in wizard", () => clickOn(page, el("NEXT button in wizard")));
    await session.step(33, "And user clicks on FINISH button in wizard", () => clickOn(page, el("FINISH button in wizard")));
    await session.step(34, "Then value of wizard readout should have text \"finished\"", () => shouldHaveText(page, el("value of wizard readout"), "finished"));
  });
});
