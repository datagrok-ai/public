/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-protac.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {matricesBuilt} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix on a PROTAC set", () => {
  const session = feature(test, "features/guides/sar-matrix-protac.feature", import.meta.url);
  test("warhead, linker and E3 ligand in one run", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(21, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(22, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(23, "When user selects \"Compound\" in Molecules input in \"SAR Matrix\" dialog", () => selectIn(page, "Compound", el("Molecules input in \"SAR Matrix\" dialog")));
    await session.step(24, "And user selects \"Solubility logS (pred)\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Solubility logS (pred)", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(25, "Then Scaling input in \"SAR Matrix\" dialog should be disabled", () => shouldBe(page, el("Scaling input in \"SAR Matrix\" dialog"), "disabled"));
    await session.step(26, "And Scaling input in \"SAR Matrix\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"SAR Matrix\" dialog"), "none"));
    await session.step(27, "When user selects \"Higher is better\" in Direction input in \"SAR Matrix\" dialog", () => selectIn(page, "Higher is better", el("Direction input in \"SAR Matrix\" dialog")));
    await session.step(28, "And user checks \"Use existing R-groups\" input in \"SAR Matrix\" dialog", () => check(page, el("\"Use existing R-groups\" input in \"SAR Matrix\" dialog")));
    await session.step(29, "Then Core input in \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("Core input in \"SAR Matrix\" dialog"), "visible"));
    await session.step(30, "When user selects \"Linker\" in Core input in \"SAR Matrix\" dialog", () => selectIn(page, "Linker", el("Core input in \"SAR Matrix\" dialog")));
    await session.step(31, "And user clicks on editor of R-groups input in \"SAR Matrix\" dialog", () => clickOn(page, el("editor of R-groups input in \"SAR Matrix\" dialog")));
    await session.step(32, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(33, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Warhead\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "Warhead"));
    await session.step(34, "And the \"text of cell 3 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"E3 Ligand\"", () => readingReads(page, "text of cell 3 of __name", el("grid viewer in \"Select columns...\" dialog"), "E3 Ligand"));
    await session.step(35, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
    await session.step(36, "And user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(37, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(39, "And user selects \"Warhead\" in \"Matrix columns\" input in \"SAR Matrix\" dialog", () => selectIn(page, "Warhead", el("\"Matrix columns\" input in \"SAR Matrix\" dialog")));
    await session.step(40, "Then \"SAR Matrix\" dialog should contain text \"Columns: Warhead\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "Columns: Warhead"));
    await session.step(41, "When user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(42, "Then SAR Matrix Viewer viewer should be visible", () => shouldBe(page, el("SAR Matrix Viewer viewer"), "visible"));
    await session.step(43, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(44, "And the \"source\" reading of SAR Matrix Viewer viewer should be \"R-group columns\"", () => readingReads(page, "source", el("SAR Matrix Viewer viewer"), "R-group columns"));
    await session.step(45, "And \"Summary\" tab should be visible", () => shouldBe(page, el("\"Summary\" tab"), "visible"));
    await session.step(46, "And tab panel should contain text \"core: Linker\"", () => shouldContainText(page, el("tab panel"), "core: Linker"));
    await session.step(47, "And tab panel should contain text \"across the matrix columns: Warhead\"", () => shouldContainText(page, el("tab panel"), "across the matrix columns: Warhead"));
    await session.step(48, "And tab panel should contain text \"folded into the row: E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "folded into the row: E3 Ligand"));
    await session.step(49, "And tab panel should contain text \"What to change\"", () => shouldContainText(page, el("tab panel"), "What to change"));
    await session.step(50, "And tab panel should contain text \"Best measured swap\"", () => shouldContainText(page, el("tab panel"), "Best measured swap"));
    await session.step(51, "When user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(52, "Then tab panel should contain text \"Warhead — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "Warhead — offsets from the additive fit"));
    await session.step(53, "When user clicks on \"E3 Ligand\" effects tab", () => clickOn(page, el("\"E3 Ligand\" effects tab")));
    await session.step(54, "Then tab panel should contain text \"E3 Ligand — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "E3 Ligand — offsets from the additive fit"));
    await session.step(55, "And tab panel should contain text \"The SAR Matrix columns enumerate Warhead\"", () => shouldContainText(page, el("tab panel"), "The SAR Matrix columns enumerate Warhead"));
    await session.step(56, "When user clicks on \"Measured in series\" effects tab", () => clickOn(page, el("\"Measured in series\" effects tab")));
    await session.step(57, "Then tab panel should contain text \"Swapping Warhead\"", () => shouldContainText(page, el("tab panel"), "Swapping Warhead"));
    await session.step(58, "And tab panel should contain text \"Swapping Linker\"", () => shouldContainText(page, el("tab panel"), "Swapping Linker"));
    await session.step(59, "And tab panel should contain text \"Swapping E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "Swapping E3 Ligand"));
    await session.step(60, "When user clicks on \"Worth making\" summary segment", () => clickOn(page, el("\"Worth making\" summary segment")));
    await session.step(61, "Then tab panel should contain text \"Nothing clears the trust gate\"", () => shouldContainText(page, el("tab panel"), "Nothing clears the trust gate"));
    await session.step(62, "When user clicks on \"Overview\" summary segment", () => clickOn(page, el("\"Overview\" summary segment")));
    await session.step(63, "And user clicks on \"SAR Matrix\" tab", () => clickOn(page, el("\"SAR Matrix\" tab")));
    await session.step(64, "Then \"SAR Matrix\" tab should be visible", () => shouldBe(page, el("\"SAR Matrix\" tab"), "visible"));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
