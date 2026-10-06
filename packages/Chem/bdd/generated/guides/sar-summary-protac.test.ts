/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/sar-summary-protac.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {callWith} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary on a PROTAC set", () => {
  const session = feature(test, "features/guides/sar-summary-protac.feature", import.meta.url);
  test("E3 ligands over linkers", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(20, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(21, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","E3 Ligand"]]), [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","E3 Ligand"]]);
    await session.step(35, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(36, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(37, "Then tab panel should contain text \"What to change\"", () => shouldContainText(page, el("tab panel"), "What to change"));
    await session.step(38, "And tab panel should contain text \"Best measured swap\"", () => shouldContainText(page, el("tab panel"), "Best measured swap"));
    await session.step(39, "When user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(40, "And user clicks on \"Linker\" effects tab", () => clickOn(page, el("\"Linker\" effects tab")));
    await session.step(41, "Then tab panel should contain text \"Linker — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "Linker — offsets from the additive fit"));
    await session.step(42, "When user clicks on \"Measured in series\" effects tab", () => clickOn(page, el("\"Measured in series\" effects tab")));
    await session.step(43, "Then tab panel should contain text \"Swapping E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "Swapping E3 Ligand"));
    await session.step(44, "When user clicks on \"Method\" summary segment", () => clickOn(page, el("\"Method\" summary segment")));
    await session.step(45, "Then tab panel should contain text \"Trust gate\"", () => shouldContainText(page, el("tab panel"), "Trust gate"));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("every component in one run", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(49, "Given user is logged in", () => loggedIn(page));
    await session.step(50, "And simple mode is off", () => simpleModeOff(page));
    await session.step(51, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(52, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(53, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","E3 Ligand"]]), [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","E3 Ligand"]]);
    await session.step(67, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(68, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(69, "And user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(70, "Then tab panel should contain text \"moves Solubility logS (pred) most: E3 Ligand\"", () => shouldContainText(page, el("tab panel"), "moves Solubility logS (pred) most: E3 Ligand"));
    await session.step(71, "And tab panel should contain text \"E3 Ligand — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "E3 Ligand — offsets from the additive fit"));
    await session.step(72, "And tab panel should contain text \"Bimodal:\"", () => shouldContainText(page, el("tab panel"), "Bimodal:"));
    await session.step(73, "And tab panel should contain text \"Cross-validated R²\"", () => shouldContainText(page, el("tab panel"), "Cross-validated R²"));
    await session.step(74, "When user clicks on \"Warhead\" effects tab", () => clickOn(page, el("\"Warhead\" effects tab")));
    await session.step(75, "Then tab panel should contain text \"Warhead — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "Warhead — offsets from the additive fit"));
    await session.step(76, "When user clicks on \"Linker\" effects tab", () => clickOn(page, el("\"Linker\" effects tab")));
    await session.step(77, "Then tab panel should contain text \"Linker — offsets from the additive fit\"", () => shouldContainText(page, el("tab panel"), "Linker — offsets from the additive fit"));
    await session.step(78, "When user clicks on \"Method\" summary segment", () => clickOn(page, el("\"Method\" summary segment")));
    await session.step(79, "Then tab panel should contain text \"Every component is ranked on Effects\"", () => shouldContainText(page, el("tab panel"), "Every component is ranked on Effects"));
    await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("warheads over linkers", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(83, "Given user is logged in", () => loggedIn(page));
    await session.step(84, "And simple mode is off", () => simpleModeOff(page));
    await session.step(85, "And user opens protac-degraders dataset", () => openDataset(page, ds("protac-degraders")));
    await session.step(86, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(87, "When user calls \"Chem:SarMatrixAnalysis\" function with:", () => callWith(page, "Chem:SarMatrixAnalysis", [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","Warhead"]]), [["table","table"],["molecules","column:Compound"],["activity","column:Caco-2 Permeability (pred)"],["scaling","none"],["activityDirection","Higher is better"],["fragmentCutoff","0.4"],["fragmentationLevels","3"],["predictVirtual","true"],["useMcsAnchors","false"],["seriesColumn",""],["coreColumn","column:Linker"],["fragmentColumns","columns:Warhead,E3 Ligand"],["columnAxis","Warhead"]]);
    await session.step(101, "Then the task bar should have finished \"Building SAR matrices\"", () => taskBarFinished(page, "Building SAR matrices"));
    await session.step(102, "When user clicks on \"Summary\" tab", () => clickOn(page, el("\"Summary\" tab")));
    await session.step(103, "Then tab panel should contain text \"What to change\"", () => shouldContainText(page, el("tab panel"), "What to change"));
    await session.step(104, "When user clicks on \"Effects\" summary segment", () => clickOn(page, el("\"Effects\" summary segment")));
    await session.step(105, "And user clicks on \"Measured in series\" effects tab", () => clickOn(page, el("\"Measured in series\" effects tab")));
    await session.step(106, "Then tab panel should contain text \"Swapping Warhead\"", () => shouldContainText(page, el("tab panel"), "Swapping Warhead"));
    await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
