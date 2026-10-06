/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-matrix-protac.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {callWith} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix on a PROTAC set", () => {
  const session = feature(test, "features/guides/sar-matrix-protac.feature", import.meta.url);
  test("warhead, linker and E3 ligand in one run", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(21, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(22, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:Compound"],["activity","column:Solubility logS (pred)"],["scaling","none"],["activityDirection","Higher is better"],["predictVirtual","true"],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","Warhead"]]), [["table","table"],["molecules","column:Compound"],["activity","column:Solubility logS (pred)"],["scaling","none"],["activityDirection","Higher is better"],["predictVirtual","true"],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","Warhead"]]);
    await session.step(32, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(33, "And \"Summary\" tab should be visible", () => shouldBe(page, el("\"Summary\" tab"), "visible"));
    await session.step(34, "Then tab panel should contain text \"core: Linker\"", () => shouldContainText(page, el("tab panel"), "core: Linker"));
    await session.step(35, "And tab panel should contain text \"across the matrix columns: Warhead\"", () => shouldContainText(page, el("tab panel"), "across the matrix columns: Warhead"));
    await session.step(36, "And tab panel should contain text \"folded into the row: E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "folded into the row: E3 Ligand"));
    await session.step(37, "And tab panel should contain text \"What to change\"", () => shouldContainText(page, el("tab panel"), "What to change"));
    await session.step(38, "And tab panel should contain text \"Best measured swap\"", () => shouldContainText(page, el("tab panel"), "Best measured swap"));
    await session.step(39, "When user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(40, "Then tab panel should contain text \"Warhead — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "Warhead — offsets from the additive fit"));
    await session.step(41, "When user clicks on \"E3 Ligand\" effects tab", () => clickOn(page, el("\"E3 Ligand\" effects tab")));
    await session.step(42, "Then tab panel should contain text \"E3 Ligand — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "E3 Ligand — offsets from the additive fit"));
    await session.step(43, "And tab panel should contain text \"The SAR Matrix columns enumerate Warhead\"", () => shouldContainText(page, el("tab panel"), "The SAR Matrix columns enumerate Warhead"));
    await session.step(44, "When user clicks on \"Measured in series\" effects tab", () => clickOn(page, el("\"Measured in series\" effects tab")));
    await session.step(45, "Then tab panel should contain text \"Swapping Warhead\"", () => shouldContainText(page, el("tab panel"), "Swapping Warhead"));
    await session.step(46, "And tab panel should contain text \"Swapping Linker\"", () => shouldContainText(page, el("tab panel"), "Swapping Linker"));
    await session.step(47, "And tab panel should contain text \"Swapping E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "Swapping E3 Ligand"));
    await session.step(48, "When user clicks on \"Worth making\" summary segment", () => clickOn(page, el("\"Worth making\" summary segment")));
    await session.step(49, "Then tab panel should contain text \"Nothing clears the trust gate\"", () => shouldContainText(page, el("tab panel"), "Nothing clears the trust gate"));
    await session.step(50, "When user clicks on \"Overview\" summary segment", () => clickOn(page, el("\"Overview\" summary segment")));
    await session.step(51, "And user clicks on \"SAR Matrix\" tab", () => clickOn(page, el("\"SAR Matrix\" tab")));
    await session.step(52, "Then \"SAR Matrix\" tab should be visible", () => shouldBe(page, el("\"SAR Matrix\" tab"), "visible"));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
