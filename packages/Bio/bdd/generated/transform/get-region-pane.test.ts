/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/get-region-pane.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.calculate.extract-region]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, selectIn, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, everyValueMatches, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Get Region pane cuts a region out of a sequence column from the context panel", () => {
  const session = feature(test, "features/transform/get-region-pane.feature", import.meta.url);
  test("The pane proposes the whole range, and Get Region adds the chosen region as a sequence column", {tag: ["@realizes:bio.calculate.extract-region"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "When user clicks on the \"header fasta\" area of grid", () => clickArea(page, "header fasta", el("grid")));
    await session.step(24, "Then the context panel should show \"fasta\"", () => contextPanelShows(page, "fasta"));
    await session.step(25, "When user expands Bioinformatics accordion header in context panel", () => expand(page, el("Bioinformatics accordion header in context panel")));
    await session.step(26, "And user expands \"Get Region\" accordion header in context panel", () => expand(page, el("\"Get Region\" accordion header in context panel")));
    await session.step(27, "Then Start input in \"Get Region\" pane in context panel should have value \"1\"", () => shouldHaveValue(page, el("Start input in \"Get Region\" pane in context panel"), "1"));
    await session.step(28, "And End input in \"Get Region\" pane in context panel should have value \"38\"", () => shouldHaveValue(page, el("End input in \"Get Region\" pane in context panel"), "38"));
    await session.step(29, "And \"Column name\" input in \"Get Region\" pane in context panel should have value \"fasta: (1-38)\"", () => shouldHaveValue(page, el("\"Column name\" input in \"Get Region\" pane in context panel"), "fasta: (1-38)"));
    await session.step(30, "When user selects \"3\" in Start input in \"Get Region\" pane in context panel", () => selectIn(page, "3", el("Start input in \"Get Region\" pane in context panel")));
    await session.step(31, "And user selects \"6\" in End input in \"Get Region\" pane in context panel", () => selectIn(page, "6", el("End input in \"Get Region\" pane in context panel")));
    await session.step(32, "Then \"Column name\" input in \"Get Region\" pane in context panel should have value \"fasta: (3-6)\"", () => shouldHaveValue(page, el("\"Column name\" input in \"Get Region\" pane in context panel"), "fasta: (3-6)"));
    await session.step(33, "When user clicks on last \"Get Region\" button in \"Get Region\" pane in context panel", () => clickOn(page, el("last \"Get Region\" button in \"Get Region\" pane in context panel")));
    await session.step(34, "Then the table should have a column \"fasta: (3-6)\"", () => hasColumn(page, "fasta: (3-6)"));
    await session.step(35, "And \"fasta: (3-6)\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta: (3-6)", "Macromolecule"));
    await session.step(36, "And the value of \"fasta: (3-6)\" column in row 1 should be \"YKET\"", () => valueInRow(page, "fasta: (3-6)", 1, "YKET"));
    await session.step(37, "And every value of \"fasta: (3-6)\" column should match \"^[A-Z]{4}$\"", () => everyValueMatches(page, "fasta: (3-6)", "^[A-Z]{4}$"));
    await session.step(38, "And the table should have 9 rows", () => rowCount(page, 9));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
});
