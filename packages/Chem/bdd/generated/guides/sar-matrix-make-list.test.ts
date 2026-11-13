/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-make-list.feature
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
import {clickOn, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Make list", () => {
  const session = feature(test, "features/guides/sar-matrix-make-list.feature", import.meta.url);
  test("Collect an analog into the Make list", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens the \"SarMatrixDemo\" project", () => openProject(page, "SarMatrixDemo"));
    await session.step(15, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(17, "When user clicks on cart of first analog worth making row", () => clickOn(page, el("cart of first analog worth making row")), undefined, "Add an analog worth making with its cart");
    await session.step(19, "And user clicks on Make list tab", () => clickOn(page, el("Make list tab")), undefined, "Open the Make list");
    await session.step(21, "Then Make list title should have text \"1 compound\"", () => shouldHaveText(page, el("Make list title"), "1 compound"), undefined, "The analog, with its potency and where it came from");
    await session.step(22, "And the \"text of cell 1 of Status\" reading of Make list grid viewer should be \"Virtual\"", () => readingReads(page, "text of cell 1 of Status", el("Make list grid viewer"), "Virtual"));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
});
