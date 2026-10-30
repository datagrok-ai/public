/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/pepti-hit/helm-edit.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openByLink, savedCellMatching, valueMatching} from '../../bindings/campaign.js';
import {helmEditorOk, helmEditorOpen, peptiHitCampaign, writeHelm} from '../../bindings/pepti-hit.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits, columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, menuLists, noBalloons, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A HELM edited in a PeptiHit campaign becomes its molecule and gets its V-iD", () => {
  const session = feature(test, "features/pepti-hit/helm-edit.feature", import.meta.url);
  test("The campaign opens with a HELM column, and a Helm cell's menu has the package's row actions", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the PeptiHit campaign \"BDD Helm\" with the stages \"Idea, Design\" and a table of:", () => peptiHitCampaign(page, "BDD Helm", "Idea, Design", [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]), [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]);
    await session.step(20, "When user opens the campaign \"BDD Helm\" by its link", () => openByLink(page, "BDD Helm"));
    await session.step(23, "Then the table should have the columns \"Helm, Molecule, Stage, V-iD\"", () => columnsExactly(page, "Helm, Molecule, Stage, V-iD"));
    await session.step(24, "And \"Helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Helm", "Macromolecule"));
    await session.step(25, "And \"Helm\" column should have units \"helm\"", () => columnUnits(page, "Helm", "helm"));
    await session.step(26, "And \"PeptiHit\" link should be visible", () => shouldBe(page, el("\"PeptiHit\" link"), "visible"));
    await session.step(27, "When user right-clicks on the \"cell 1 of Helm\" area of grid", () => rightClickArea(page, "cell 1 of Helm", el("grid")));
    await session.step(28, "Then the open menu should list \"Add new row\"", () => menuLists(page, "Add new row"));
    await session.step(29, "And the open menu should list \"Duplicate Sequence\"", () => menuLists(page, "Duplicate Sequence"));
    await session.step(30, "And the open menu should list \"Re-Run Calculations\"", () => menuLists(page, "Re-Run Calculations"));
  });
  test("A HELM committed in the HELM editor fills the Molecule, the template's property and the V-iD, and is saved", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the PeptiHit campaign \"BDD Helm\" with the stages \"Idea, Design\" and a table of:", () => peptiHitCampaign(page, "BDD Helm", "Idea, Design", [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]), [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]);
    await session.step(20, "When user opens the campaign \"BDD Helm\" by its link", () => openByLink(page, "BDD Helm"));
    await session.step(33, "When user double-clicks on the \"cell 1 of Helm\" area of grid", () => doubleClickArea(page, "cell 1 of Helm", el("grid")));
    await session.step(34, "Then the HELM editor should be open", () => helmEditorOpen(page));
    await session.step(35, "When user writes the HELM \"PEPTIDE1{W.Y.W.Y.W.Y.K}$$$$V2.0\" in the HELM editor", () => writeHelm(page, "PEPTIDE1{W.Y.W.Y.W.Y.K}$$$$V2.0"));
    await session.step(36, "And user clicks on OK in the HELM editor", () => helmEditorOk(page));
    await session.step(37, "Then the value of \"Helm\" column in row 1 should match \"^PEPTIDE1[{]W[.]Y[.]W[.]Y[.]W[.]Y[.]K[}]\"", () => valueMatching(page, "Helm", 1, "^PEPTIDE1[{]W[.]Y[.]W[.]Y[.]W[.]Y[.]K[}]"));
    await session.step(38, "And the campaign \"BDD Helm\" should be saved with a value matching \"^1194[.]3\" in row 1 of the \"MW\" column", () => savedCellMatching(page, "BDD Helm", "^1194[.]3", 1, "MW"));
    await session.step(39, "And the value of \"MW\" column in row 1 should match \"^1194[.]3\"", () => valueMatching(page, "MW", 1, "^1194[.]3"));
    await session.step(40, "And the value of \"Molecule\" column in row 1 should match \"M  END\"", () => valueMatching(page, "Molecule", 1, "M  END"));
    await session.step(41, "And the value of \"V-iD\" column in row 1 should match \"^V[0-9]{6}$\"", () => valueMatching(page, "V-iD", 1, "^V[0-9]{6}$"));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The campaign view is named after the campaign from the first open", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the PeptiHit campaign \"BDD Helm\" with the stages \"Idea, Design\" and a table of:", () => peptiHitCampaign(page, "BDD Helm", "Idea, Design", [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]), [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","","Idea",""]]);
    await session.step(20, "When user opens the campaign \"BDD Helm\" by its link", () => openByLink(page, "BDD Helm"));
    await session.step(45, "Then \"PeptiHit\" breadcrumb should be visible", () => shouldBe(page, el("\"PeptiHit\" breadcrumb"), "visible"));
    await session.step(46, "And the \"BDD Helm\" view should be current", () => viewIsCurrent(page, "BDD Helm"));
    await session.step(47, "And \"BDD Helm\" breadcrumb should be visible", () => shouldBe(page, el("\"BDD Helm\" breadcrumb"), "visible"));
  });
});
