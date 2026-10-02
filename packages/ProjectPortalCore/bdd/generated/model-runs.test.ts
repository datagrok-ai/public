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
import {exampleRunSaved, openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
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
});
