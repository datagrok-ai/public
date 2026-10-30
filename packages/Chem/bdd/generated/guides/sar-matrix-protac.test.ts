/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-protac.feature
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
import {check, clickOn, scrollIntoView, selectIn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix on a PROTAC set", () => {
  const session = feature(test, "features/guides/sar-matrix-protac.feature", import.meta.url);
  test("Warhead, linker and E3 ligand in one run", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(17, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(18, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(19, "When user selects \"Compound\" in Molecules input in \"SAR Matrix\" dialog", () => selectIn(page, "Compound", el("Molecules input in \"SAR Matrix\" dialog")));
    await session.step(20, "And user selects \"Solubility logS (pred)\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Solubility logS (pred)", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(22, "Then Scaling input in \"SAR Matrix\" dialog should be disabled", () => shouldBe(page, el("Scaling input in \"SAR Matrix\" dialog"), "disabled"), undefined, "logS is negative, so it cannot be log-scaled");
    await session.step(23, "And Scaling input in \"SAR Matrix\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"SAR Matrix\" dialog"), "none"));
    await session.step(24, "When user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
    await session.step(25, "And user checks \"Use existing R-groups\" input in \"SAR Matrix\" dialog", () => check(page, el("\"Use existing R-groups\" input in \"SAR Matrix\" dialog")));
    await session.step(26, "Then Core input in \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("Core input in \"SAR Matrix\" dialog"), "visible"));
    await session.step(27, "When user selects \"Linker\" in Core input in \"SAR Matrix\" dialog", () => selectIn(page, "Linker", el("Core input in \"SAR Matrix\" dialog")));
    await session.step(28, "And user clicks on editor of R-groups input in \"SAR Matrix\" dialog", () => clickOn(page, el("editor of R-groups input in \"SAR Matrix\" dialog")));
    await session.step(29, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(30, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
    await session.step(31, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(32, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(34, "And user selects \"Warhead\" in \"Matrix columns\" input in \"SAR Matrix\" dialog", () => selectIn(page, "Warhead", el("\"Matrix columns\" input in \"SAR Matrix\" dialog")));
    await session.step(35, "And user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(36, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(38, "When user scrolls E3 Ligand component row into view", () => scrollIntoView(page, el("E3 Ligand component row")), undefined, "Scroll down to What to change");
    await session.step(40, "Then E3 Ligand component row should be visible", () => shouldBe(page, el("E3 Ligand component row"), "visible"), undefined, "What to change ranks the warhead, the linker and the E3 ligand in one fit");
    await session.step(41, "And What to change title should have text \"What to change · 3 components\"", () => shouldHaveText(page, el("What to change title"), "What to change · 3 components"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
});
