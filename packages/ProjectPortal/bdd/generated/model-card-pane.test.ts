/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/model-card-pane.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {selectInModelHub} from '../bindings/steps.js';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A catalogued model's card in the context panel", () => {
  const session = feature(test, "features/model-card-pane.feature", import.meta.url);
  test("Select the example model in Model Hub, read its card in the context panel", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "When user selects \"Fed-batch titer (example)\" in Model Hub", () => selectInModelHub(page, "Fed-batch titer (example)"));
    await session.step(8, "Then model card pane portfolios should contain text \"In the catalogue of EX-BIOC\"", () => shouldContainText(page, el("model card pane portfolios"), "In the catalogue of EX-BIOC"));
    await session.step(9, "And model card pane should contain text \"What final titer\"", () => shouldContainText(page, el("model card pane"), "What final titer"));
    await session.step(10, "And model card pane should contain text \"draft\"", () => shouldContainText(page, el("model card pane"), "draft"));
  });
});
