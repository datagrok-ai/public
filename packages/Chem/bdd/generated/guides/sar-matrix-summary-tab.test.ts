/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-summary-tab.feature
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
import {clickOn, scrollIntoView, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary tab", () => {
  const session = feature(test, "features/guides/sar-matrix-summary-tab.feature", import.meta.url);
  test("Read the Summary", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens the \"SarMatrixDemo\" project", () => openProject(page, "SarMatrixDemo"));
    await session.step(15, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(17, "When user clicks on first key series row", () => clickOn(page, el("first key series row")), undefined, "Open a series to see its best core");
    await session.step(18, "Then key series detail should be visible", () => shouldBe(page, el("key series detail"), "visible"));
    await session.step(20, "When user scrolls What improves potency answer tile into view", () => scrollIntoView(page, el("What improves potency answer tile")), undefined, "Scroll down to the answers");
    await session.step(22, "And user clicks on What improves potency answer tile", () => clickOn(page, el("What improves potency answer tile")), undefined, "Open an answer to see the evidence behind it");
    await session.step(23, "Then Summary evidence should contain text \"What buys potency\"", () => shouldContainText(page, el("Summary evidence"), "What buys potency"));
    await session.step(24, "And first evidence row should contain text \"weakest pair\"", () => shouldContainText(page, el("first evidence row"), "weakest pair"));
    await session.step(26, "When user clicks on first analog worth making row", () => clickOn(page, el("first analog worth making row")), undefined, "Open an analog to see what its prediction adds up from");
    await session.step(27, "Then analog worth making detail should be visible", () => shouldBe(page, el("analog worth making detail"), "visible"));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
});
