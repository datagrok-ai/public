/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/workspace-model-card.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {modelPinned, openApp, selectInWorkspace} from '../bindings/steps.js';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A catalogued model's card beside its preview in the Workspace on Home", () => {
  const session = feature(test, "features/workspace-model-card.feature", import.meta.url);
  test("Pin the example model, select it in the Workspace, read its card", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given the user has pinned the \"ProjectPortalCore:FedBatchTiter\" model", () => modelPinned(page, "ProjectPortalCore:FedBatchTiter"));
    await session.step(12, "When user selects the \"ProjectPortalCore:FedBatchTiter\" model in the Workspace on Home", () => selectInWorkspace(page, "ProjectPortalCore:FedBatchTiter"));
    await session.step(13, "Then workspace model card should contain text \"What final titer\"", () => shouldContainText(page, el("workspace model card"), "What final titer"));
    await session.step(14, "And workspace model card should contain text \"draft\"", () => shouldContainText(page, el("workspace model card"), "draft"));
    await session.step(15, "And workspace model card should contain text \"EX-BIOC\"", () => shouldContainText(page, el("workspace model card"), "EX-BIOC"));
  });
});
