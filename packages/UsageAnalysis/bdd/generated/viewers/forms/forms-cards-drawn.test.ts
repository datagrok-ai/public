/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-cards-drawn.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer — every selected row's card draws its own molecule and curve", () => {
  const session = feature(test, "features/viewers/forms/forms-cards-drawn.feature", import.meta.url);
  test("With four rows selected, each card draws its own molecule and its own curve", {tag: ["@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(16, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(17, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["fieldsColumnNames","smiles, multiple prefit"]]), [["fieldsColumnNames","smiles, multiple prefit"]]);
    await session.step(19, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(20, "When user selects the first 4 rows", () => selectFirstRows(page, 4));
    await session.step(21, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(22, "And the \"card kind of card 1\" reading of forms viewer should be \"current\"", () => readingReads(page, "card kind of card 1", el("forms viewer"), "current"));
    await session.step(23, "And the \"card kind of card 2\" reading of forms viewer should be \"mouse-over\"", () => readingReads(page, "card kind of card 2", el("forms viewer"), "mouse-over"));
    await session.step(24, "And the \"record of card 2\" reading of forms viewer should be \"\"", () => readingReads(page, "record of card 2", el("forms viewer"), ""));
    await session.step(25, "And the \"record of card 3\" reading of forms viewer should be \"1\"", () => readingReads(page, "record of card 3", el("forms viewer"), "1"));
    await session.step(26, "And the \"record of card 4\" reading of forms viewer should be \"2\"", () => readingReads(page, "record of card 4", el("forms viewer"), "2"));
    await session.step(27, "And the \"record of card 5\" reading of forms viewer should be \"3\"", () => readingReads(page, "record of card 5", el("forms viewer"), "3"));
    await session.step(28, "And the \"record of card 6\" reading of forms viewer should be \"4\"", () => readingReads(page, "record of card 6", el("forms viewer"), "4"));
    await session.step(29, "And the \"smiles of card 3\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "smiles of card 3", el("forms viewer"), "canvas"));
    await session.step(30, "And the \"multiple prefit of card 3\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "multiple prefit of card 3", el("forms viewer"), "canvas"));
    await session.step(31, "And the \"smiles of card 4\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "smiles of card 4", el("forms viewer"), "canvas"));
    await session.step(32, "And the \"multiple prefit of card 4\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "multiple prefit of card 4", el("forms viewer"), "canvas"));
    await session.step(33, "And the \"smiles of card 5\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "smiles of card 5", el("forms viewer"), "canvas"));
    await session.step(34, "And the \"multiple prefit of card 5\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "multiple prefit of card 5", el("forms viewer"), "canvas"));
    await session.step(35, "And the \"smiles of card 6\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "smiles of card 6", el("forms viewer"), "canvas"));
    await session.step(36, "And the \"multiple prefit of card 6\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "multiple prefit of card 6", el("forms viewer"), "canvas"));
    await session.step(37, "And forms viewer should have a \"field multiple prefit of card 6\" area", () => hasArea(page, el("forms viewer"), "field multiple prefit of card 6"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
