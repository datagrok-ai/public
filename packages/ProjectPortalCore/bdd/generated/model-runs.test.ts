/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/model-runs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {exampleRunSaved, openApp, runExampleModel, scratchProjectAtStage, scratchWorkspace} from '../bindings/steps.js';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Saved runs on the program page, and the model card", () => {
  const session = feature(test, "features/model-runs.feature", import.meta.url);
  test("A saved run opens its model's card", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-RUNS\" at stage \"D5\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-RUNS", "D5", "EX-BIOC"));
    await session.step(13, "And a run of the example model saved to the scratch project \"BDD-RUNS\" in the \"EX-BIOC\" portfolio", () => exampleRunSaved(page, "BDD-RUNS", "EX-BIOC"));
    await session.step(14, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(15, "And user clicks on scratch runs project card", () => clickOn(page, el("scratch runs project card")));
    await session.step(16, "Then saved runs should contain text \"Fed-batch titer (example) · BDD-RUNS\"", () => shouldContainText(page, el("saved runs"), "Fed-batch titer (example) · BDD-RUNS"));
    await session.step(17, "And saved runs should contain text \"Production bioreactor\"", () => shouldContainText(page, el("saved runs"), "Production bioreactor"));
    await session.step(18, "When user clicks on example model link", () => clickOn(page, el("example model link")));
    await session.step(19, "Then model card question should contain text \"What final titer\"", () => shouldContainText(page, el("model card question"), "What final titer"));
    await session.step(20, "And model card inputs should contain text \"33 – 37.5\"", () => shouldContainText(page, el("model card inputs"), "33 – 37.5"));
  });
  test("A run is saved to a program from the dialog that follows it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(23, "Given a scratch project \"BDD-SAVE\" at stage \"D5\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-SAVE", "D5", "EX-BIOC"));
    await session.step(24, "And the workspace of the scratch project \"BDD-SAVE\" in the \"EX-BIOC\" portfolio is built", () => scratchWorkspace(page, "BDD-SAVE", "EX-BIOC"));
    await session.step(25, "When user runs the example model", () => runExampleModel(page));
    await session.step(26, "And user clicks on save run link", () => clickOn(page, el("save run link")));
    await session.step(27, "Then save run dialog should be visible", () => shouldBe(page, el("save run dialog"), "visible"));
    await session.step(28, "And run portfolio input should have value \"Example Biologics Development\"", () => shouldHaveValue(page, el("run portfolio input"), "Example Biologics Development"));
    await session.step(29, "And run unit operation input should have value \"Production bioreactor\"", () => shouldHaveValue(page, el("run unit operation input"), "Production bioreactor"));
    await session.step(30, "And run type input should have value \"Simulated\"", () => shouldHaveValue(page, el("run type input"), "Simulated"));
    await session.step(31, "When user selects \"BDD-SAVE\" in run program input", () => selectIn(page, "BDD-SAVE", el("run program input")));
    await session.step(32, "And user selects \"characterization\" in run purpose input", () => selectIn(page, "characterization", el("run purpose input")));
    await session.step(33, "And user clicks on save run button", () => clickOn(page, el("save run button")));
    await session.step(34, "Then save run dialog should not be visible", () => shouldNotBe(page, el("save run dialog"), "visible"));
    await session.step(35, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(36, "And user clicks on scratch save project card", () => clickOn(page, el("scratch save project card")));
    await session.step(37, "Then saved runs should contain text \"Fed-batch titer (example) · BDD-SAVE\"", () => shouldContainText(page, el("saved runs"), "Fed-batch titer (example) · BDD-SAVE"));
    await session.step(38, "And saved runs should contain text \"characterization\"", () => shouldContainText(page, el("saved runs"), "characterization"));
  });
});
