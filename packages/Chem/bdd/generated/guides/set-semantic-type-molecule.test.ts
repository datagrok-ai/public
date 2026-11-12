/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/set-semantic-type-molecule.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, pressKey, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openTableOf, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {areaColors, clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Draw an undetected SMILES column as molecules", () => {
  const session = feature(test, "features/guides/set-semantic-type-molecule.feature", import.meta.url);
  test("Set the semantic type of a SMILES column to Molecule", {tag: ["@guide", "@help:govern/catalog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"Reagents\" with:", () => openTableOf(page, "Reagents", [["Reagent","Equivalents"],["CCO","1.0"],["CC#N","2.0"],["CC=O","1.5"],["OCCO","0.5"],["CCCl","1.2"],["NCCN","3.0"]]), [["Reagent","Equivalents"],["CCO","1.0"],["CC#N","2.0"],["CC=O","1.5"],["OCCO","0.5"],["CCCl","1.2"],["NCCN","3.0"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page), undefined, "The SMILES are shown as plain text");
    await session.step(26, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "Then \"Reagent\" column should have semantic type \"\"", () => columnSemType(page, "Reagent", ""), undefined, "Reagent has no semantic type");
    await session.step(30, "When user clicks on the \"header Reagent\" area of grid", () => clickArea(page, "header Reagent", el("grid")), undefined, "Click the Reagent column header");
    await session.step(31, "Then the context panel should show \"Reagent\"", () => contextPanelShows(page, "Reagent"));
    await session.step(33, "When user expands Details accordion header in context panel", () => expand(page, el("Details accordion header in context panel")), undefined, "Expand Details in the Context Panel");
    await session.step(35, "When user types \"Molecule\" into \"Search or enter semantic type\" input in Details pane in context panel", () => typeInto(page, "Molecule", el("\"Search or enter semantic type\" input in Details pane in context panel")), undefined, "Type Molecule into the Semantic type field");
    await session.step(36, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(38, "Then \"Reagent\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Reagent", "Molecule"), undefined, "Reagent is now a Molecule column");
    await session.step(39, "And \"Reagent\" column should have tag \"cell.renderer\" equal to \"Molecule\"", () => columnTag(page, "Reagent", "cell.renderer", "Molecule"));
    await session.step(41, "And the \"cell 1 of Reagent\" area of grid should be painted in at least 1 colors", () => areaColors(page, "cell 1 of Reagent", el("grid"), 1), undefined, "The grid draws the molecules");
  });
});
