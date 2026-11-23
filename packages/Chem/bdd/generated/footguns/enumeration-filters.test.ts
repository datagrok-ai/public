/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/enumeration-filters.feature
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
import {changeBlocking, enumerateDisabled, enumerateEnabled, enumerateNoWords, enumerateWordsNotName, enumerateWordsUnreadable, enumerateWordsWarnNone, enumeratorBuildingBlocks, enumeratorExclusion, enumeratorTemplates, openEnumerator, pageEngine} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {autostartsCompleted} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in the reaction enumerator's filters", () => {
  const session = feature(test, "features/footguns/enumeration-filters.feature", import.meta.url);
  test("A blocking SMARTS no engine reads blocks Enumerate and says which", {tag: ["@footguns", "@batch-5", "@catalogue-TRF-05"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(16, "And user opens the Reaction Enumerator from \"Chem | Reactions\"", () => openEnumerator(page, "Chem | Reactions"));
    await session.step(20, "Given the enumerator's templates are:", () => enumeratorTemplates(page, [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[N;H2","Amide"]]), [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[N;H2","Amide"]]);
    await session.step(23, "And the enumerator's building blocks are \"CC(=O)O, NCCN, CCN\"", () => enumeratorBuildingBlocks(page, "CC(=O)O, NCCN, CCN"));
    await session.step(24, "Then the Enumerate button is disabled", () => enumerateDisabled(page));
    await session.step(25, "And the words beside the Enumerate button name \"[N;H2\" and say it cannot be read", () => enumerateWordsUnreadable(page, "[N;H2"));
    await session.step(27, "When user changes the blocking SMARTS of template 1 to \"[N&H2]\"", () => changeBlocking(page, 1, "[N&H2]"));
    await session.step(28, "Then the Enumerate button is enabled", () => enumerateEnabled(page));
  });
  test("An exclusion SMARTS no engine reads blocks Enumerate and says which", {tag: ["@footguns", "@batch-5", "@catalogue-TRF-06"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(16, "And user opens the Reaction Enumerator from \"Chem | Reactions\"", () => openEnumerator(page, "Chem | Reactions"));
    await session.step(32, "Given the enumerator's templates are:", () => enumeratorTemplates(page, [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","","Amide"]]), [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","","Amide"]]);
    await session.step(35, "And the enumerator's building blocks are \"CC(=O)O, CCN\"", () => enumeratorBuildingBlocks(page, "CC(=O)O, CCN"));
    await session.step(36, "And the enumerator's exclusion SMARTS are \"C(=O)Cl], S(=O)(=O)[#6]\"", () => enumeratorExclusion(page, "C(=O)Cl], S(=O)(=O)[#6]"));
    await session.step(37, "Then the Enumerate button is disabled", () => enumerateDisabled(page));
    await session.step(38, "And the words beside the Enumerate button name \"C(=O)Cl]\" and say it cannot be read", () => enumerateWordsUnreadable(page, "C(=O)Cl]"));
    await session.step(39, "And the words beside the Enumerate button do not name \"S(=O)(=O)[#6]\"", () => enumerateWordsNotName(page, "S(=O)(=O)[#6]"));
  });
  test("A blocking SMARTS that matches none of the building blocks is a warn before the run", {tag: ["@footguns", "@batch-5", "@catalogue-TRF-07"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(16, "And user opens the Reaction Enumerator from \"Chem | Reactions\"", () => openEnumerator(page, "Chem | Reactions"));
    await session.step(43, "Given the enumerator's templates are:", () => enumeratorTemplates(page, [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[NX3](=O)=O","Amide"]]), [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[NX3](=O)=O","Amide"]]);
    await session.step(46, "And the enumerator's building blocks are \"O=C(O)c1ccc([N+](=O)[O-])cc1, CCN\"", () => enumeratorBuildingBlocks(page, "O=C(O)c1ccc([N+](=O)[O-])cc1, CCN"));
    await session.step(47, "Then the Enumerate button is enabled", () => enumerateEnabled(page));
    await session.step(48, "And the words beside the Enumerate button warn that \"[NX3](=O)=O\" matches none of the building blocks", () => enumerateWordsWarnNone(page, "[NX3](=O)=O"));
  });
  test("RDKit mode: the enumerator behaves as before", {tag: ["@footguns", "@batch-5", "@catalogue-UI-03", "@catalogue-TRF-05", "@catalogue-TRF-07"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(16, "And user opens the Reaction Enumerator from \"Chem | Reactions\"", () => openEnumerator(page, "Chem | Reactions"));
    await session.step(52, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(53, "And the enumerator's templates are:", () => enumeratorTemplates(page, [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[N;H2","Amide"]]), [["reaction_smarts","blocking_fg","reaction_name"],["[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]","[N;H2","Amide"]]);
    await session.step(56, "And the enumerator's building blocks are \"CC(=O)O, NCCN, CCN\"", () => enumeratorBuildingBlocks(page, "CC(=O)O, NCCN, CCN"));
    await session.step(57, "Then the Enumerate button is enabled", () => enumerateEnabled(page));
    await session.step(58, "And no foot-gun words show beside the Enumerate button", () => enumerateNoWords(page));
  });
});
