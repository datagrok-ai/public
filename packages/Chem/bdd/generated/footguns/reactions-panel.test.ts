/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/reactions-panel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {checksWords, fullScreenStep, mappingAnswers, mappingConfidence, mappingNoWarn, mappingNotShow, mappingWarn, mappingWarn2, openFullScreen, paneStep, paneSteps, reactionTable} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in reaction cells, the Reaction panel and its atom mapping (Crux mode)", () => {
  const session = feature(test, "features/footguns/reactions-panel.feature", import.meta.url);
  test("A reaction Crux reads no reaction from, though RDKit does, is a block in the Checks section", {tag: ["@footguns", "@batch-2", "@catalogue-RXN-01", "@catalogue-RXN-10"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "Given user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|"]]), [["reaction"],["CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|"]]);
    await session.step(24, "When user clicks on the \"cell 1 of reaction\" area of grid", () => clickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(25, "Then the context panel's Checks section lists a block mentioning \"f:\" and \"Chem reads no reaction\"", () => checksWords(page, "block", "f:", "Chem reads no reaction"));
  });
  test("The panel's steps are the value's own steps, written as the cell writes them", {tag: ["@footguns", "@batch-2", "@catalogue-RXN-13"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(29, "Given user opens a table \"routes\" with a reaction column \"reaction\" of units \"route\":", () => reactionTable(page, "routes", "reaction", "route", [["reaction"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]--**--[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"]]), [["reaction"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]--**--[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"]]);
    await session.step(32, "When user clicks on the \"cell 1 of reaction\" area of grid", () => clickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(33, "Then the Reaction panel shows 2 steps", () => paneSteps(page, 2));
    await session.step(34, "And step 1 of the Reaction panel draws \"[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]\"", () => paneStep(page, 1, "[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]"));
    await session.step(35, "And step 2 of the Reaction panel draws \"[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C\"", () => paneStep(page, 2, "[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"));
    await session.step(36, "When user opens the Reaction panel's full-screen view", () => openFullScreen(page));
    await session.step(37, "Then the full-screen view's step 1 draws \"[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]\"", () => fullScreenStep(page, 1, "[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]"));
  });
  test("A mapping under 0.5 confidence says to check it; one above does not", {tag: ["@footguns", "@batch-2", "@catalogue-MAP-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(41, "Given user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["CC(C)(C)OC(=O)N1CCNCC1.O=C(O)C1CC1>CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F.CCN(C(C)C)C(C)C.CN(C)C=O>CC(C)(C)OC(=O)N1CCN(C(=O)C2CC2)CC1 |f:2.3|"],["COCCN.O=Cc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1>CC(=O)O[BH-](OC(C)=O)OC(C)=O.[Na+].CC(=O)O.ClCCCl>COCCNCc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1 |f:2.3|"]]), [["reaction"],["CC(C)(C)OC(=O)N1CCNCC1.O=C(O)C1CC1>CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F.CCN(C(C)C)C(C)C.CN(C)C=O>CC(C)(C)OC(=O)N1CCN(C(=O)C2CC2)CC1 |f:2.3|"],["COCCN.O=Cc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1>CC(=O)O[BH-](OC(C)=O)OC(C)=O.[Na+].CC(=O)O.ClCCCl>COCCNCc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1 |f:2.3|"]]);
    await session.step(45, "When user clicks on the \"cell 1 of reaction\" area of grid", () => clickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(46, "And the Reaction panel's mapping answers", () => mappingAnswers(page));
    await session.step(47, "Then the Reaction panel's Mapping section shows the confidence \"0.30\"", () => mappingConfidence(page, "0.30"));
    await session.step(48, "And the Reaction panel's Mapping section shows a warn mentioning \"check this mapping\"", () => mappingWarn(page, "check this mapping"));
    await session.step(49, "When user clicks on the \"cell 2 of reaction\" area of grid", () => clickArea(page, "cell 2 of reaction", el("grid")));
    await session.step(50, "And the Reaction panel's mapping answers", () => mappingAnswers(page));
    await session.step(51, "Then the Reaction panel's Mapping section shows the confidence \"0.94\"", () => mappingConfidence(page, "0.94"));
    await session.step(52, "And the Reaction panel's Mapping section shows no warn", () => mappingNoWarn(page));
  });
  test("Mapping a template says it is a template, never a bare SMILES error [template=[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]]", {tag: ["@footguns", "@batch-2", "@catalogue-MAP-02"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(56, "Given user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]"]]), [["template"],["[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]"]]);
    await session.step(59, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(60, "Then the Reaction panel's Mapping section shows a warn mentioning \"template\" and \"map numbers\"", () => mappingWarn2(page, "template", "map numbers"));
    await session.step(61, "When the Reaction panel's mapping answers", () => mappingAnswers(page));
    await session.step(62, "Then the Reaction panel's Mapping section does not show \"is not a valid SMILES string\"", () => mappingNotShow(page, "is not a valid SMILES string"));
  });
});
