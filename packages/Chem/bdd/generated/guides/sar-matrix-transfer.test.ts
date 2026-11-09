/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-transfer.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {matricesBuilt} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Transfer tab", () => {
  const session = feature(test, "features/guides/sar-matrix-transfer.feature", import.meta.url);
  test("Switch a transfer to its other target", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens the \"SarMatrixDemo\" project", () => openProject(page, "SarMatrixDemo"));
    await session.step(16, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(18, "When user clicks on SAR Transfer tab", () => clickOn(page, el("SAR Transfer tab")), undefined, "Open the SAR Transfer tab: each card is a source core");
    await session.step(19, "Then \"Transfer to\" input should have value \"Series 33 · Core 1 · R1\"", () => shouldHaveValue(page, el("\"Transfer to\" input"), "Series 33 · Core 1 · R1"));
    await session.step(21, "When user selects \"Series 98 · Core 1 · R1\" in \"Transfer to\" input", () => selectIn(page, "Series 98 · Core 1 · R1", el("\"Transfer to\" input")), undefined, "Switch the source to its other target");
    await session.step(22, "Then transfer title should have text \"Series 69 · Core 3 → Series 98 · Core 1\"", () => shouldHaveText(page, el("transfer title"), "Series 69 · Core 3 → Series 98 · Core 1"));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
});
