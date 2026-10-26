/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-tab.feature
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
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, openProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix tab", () => {
  const session = feature(test, "features/guides/sar-matrix-tab.feature", import.meta.url);
  test("Open a series and one of its cells", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(16, "And user opens the \"SarMatrixDemo\" project", () => openProject(page, "SarMatrixDemo"));
    await session.step(17, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(18, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
    await session.step(19, "Then Series 69 series card should be visible", () => shouldBe(page, el("Series 69 series card"), "visible"));
    await session.step(21, "When user clicks on Series 1 series card", () => clickOn(page, el("Series 1 series card")), undefined, "Click a card to open that series");
    await session.step(22, "Then the \"series\" reading of SAR Matrix Viewer viewer should be \"Series 1\"", () => readingReads(page, "series", el("SAR Matrix Viewer viewer"), "Series 1"));
    await session.step(24, "When user clicks on the \"cell 1 of c0\" area of grid viewer", () => clickArea(page, "cell 1 of c0", el("grid viewer")), undefined, "Click a cell to open its compound in the Context Panel");
    await session.step(25, "Then context panel should contain text \"SAR analysis\"", () => shouldContainText(page, el("context panel"), "SAR analysis"));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
});
