/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix.feature
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
import {clickOn, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix on a compound set", () => {
  const session = feature(test, "features/guides/sar-matrix.feature", import.meta.url);
  test("Run a SAR Matrix analysis", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens sar-matrix-demo dataset", () => openDataset(page, ds("sar-matrix-demo")));
    await session.step(14, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(15, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(16, "When user selects \"smiles\" in Molecules input in \"SAR Matrix\" dialog", () => selectIn(page, "smiles", el("Molecules input in \"SAR Matrix\" dialog")));
    await session.step(17, "And user selects \"hERG_pIC50\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "hERG_pIC50", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(19, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")), undefined, "hERG pIC50 is already a log value, so no scaling");
    await session.step(20, "Then Scaling input in \"SAR Matrix\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"SAR Matrix\" dialog"), "none"));
    await session.step(21, "When user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
    await session.step(22, "And user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(23, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(25, "And Key series card should be visible", () => shouldBe(page, el("Key series card"), "visible"), undefined, "The analysis opens on the Summary");
    await session.step(26, "And Worth making card should be visible", () => shouldBe(page, el("Worth making card"), "visible"));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
});
