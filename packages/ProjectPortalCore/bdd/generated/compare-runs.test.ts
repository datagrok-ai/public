/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compare-runs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, openTableStarting, scratchProjectAtStage, titledRunSaved} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Compare a criterion's runs in Compute2", () => {
  const session = feature(test, "features/compare-runs.feature", import.meta.url);
  test("The process description's two candidate runs, listed for Compare Runs", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-CMP\" at stage \"D4\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-CMP", "D4", "EX-BIOC"));
    await session.step(12, "And a run of the example model at feed rate 2 titled \"Feed 2\" saved to the scratch project \"BDD-CMP\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 2, "Feed 2", "BDD-CMP", "EX-BIOC"));
    await session.step(13, "And a run of the example model at feed rate 4 titled \"Feed 4\" saved to the scratch project \"BDD-CMP\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 4, "Feed 4", "BDD-CMP", "EX-BIOC"));
    await session.step(14, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(15, "And user clicks on scratch compare project card", () => clickOn(page, el("scratch compare project card")));
    await session.step(16, "And user clicks on compare process description runs button", () => clickOn(page, el("compare process description runs button")));
    await session.step(17, "Then compare runs dialog should be visible", () => shouldBe(page, el("compare runs dialog"), "visible"));
    await session.step(18, "And compare runs summary should contain text \"2 runs of 1 model\"", () => shouldContainText(page, el("compare runs summary"), "2 runs of 1 model"));
    await session.step(19, "And compare runs groups should contain text \"Fed-batch titer (example)\"", () => shouldContainText(page, el("compare runs groups"), "Fed-batch titer (example)"));
    await session.step(20, "And compare runs groups should contain text \"Feed 4\"", () => shouldContainText(page, el("compare runs groups"), "Feed 4"));
    await session.step(21, "When user clicks on open compare runs as table button", () => clickOn(page, el("open compare runs as table button")));
    await session.step(22, "Then the open table should have 2 rows and a name starting with \"BDD-CMP · process_description runs\"", () => openTableStarting(page, 2, "BDD-CMP · process_description runs"));
  });
});
