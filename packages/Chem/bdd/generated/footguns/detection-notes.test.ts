/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/detection-notes.feature
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
import {checksWord, checksWords, detectedTable} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Columns Chem's detector reads otherwise than their author meant, and the notes for usual oddities", () => {
  const session = feature(test, "features/footguns/detection-notes.feature", import.meta.url);
  test("A column of templates with one route in it: the route is flagged and the templates are still checked", {tag: ["@footguns", "@batch-3", "@catalogue-RXN-12"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(20, "Given user opens a CSV \"templates\" with a column \"reaction\" whose semantic types are detected:", () => detectedTable(page, "templates", "reaction", [["reaction"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"],["[c:1][Br].[#5]([OH])([OH])[c:2]>>[c:1][c:2]"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"]]), [["reaction"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"],["[c:1][Br].[#5]([OH])([OH])[c:2]>>[c:1][c:2]"],["[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"]]);
    await session.step(25, "Then \"reaction\" column should have semantic type \"ChemicalReaction\"", () => columnSemType(page, "reaction", "ChemicalReaction"));
    await session.step(26, "When user clicks on the \"cell 3 of reaction\" area of grid", () => clickArea(page, "cell 3 of reaction", el("grid")));
    await session.step(27, "Then the context panel's Checks section lists a warn mentioning \"route\"", () => checksWord(page, "warn", "route"));
    await session.step(28, "When user clicks on the \"cell 1 of reaction\" area of grid", () => clickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(29, "Then the context panel's Checks section lists a warn mentioning \"map number 1\" and \"twice\"", () => checksWords(page, "warn", "map number 1", "twice"));
  });
});
