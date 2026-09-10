/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-renderers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {everyRecordCardShows, recordCardRows, recordCardsAreSelection} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCategoricalFilter, clearSelection, filterPasses, filterPassesFewer, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, hasArea, noBalloons, noErrors, propertyShouldBe, readingAtLeast, readingHigher, readingIs, readingLower, readingReads, readingsEqual, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer renderers, renderer size and the twenty-field cap", () => {
  const session = feature(test, "features/viewers/forms/forms-renderers.feature", import.meta.url);
  test("Forms viewer renderers, renderer size and the twenty-field cap", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(22, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(23, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(24, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(25, "And \"Structure\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Structure", "Molecule"));
    await session.step(26, "And \"Core\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Core", "Molecule"));
    await run.scenario("The field set stops at twenty columns, in table order, with no message", async () => {
      await session.step(29, "Then the \"fields shown\" reading of forms viewer should be 20", () => readingIs(page, "fields shown", el("forms viewer"), 20));
      await session.step(30, "And the \"fields\" reading of forms viewer should be \"Id, Structure, CAST Idea ID, Last Published Date, Chemist, Lab Notebook, Stereo Category, Series, Scaffold Names, Primary Series Name, Primary Scaffold Name, Has Unlabeled R-Groups, Core, R1, R2, R3, R100, R101, Chemical Space X, Chemical Space Y\"", () => readingReads(page, "fields", el("forms viewer"), "Id, Structure, CAST Idea ID, Last Published Date, Chemist, Lab Notebook, Stereo Category, Series, Scaffold Names, Primary Series Name, Primary Scaffold Name, Has Unlabeled R-Groups, Core, R1, R2, R3, R100, R101, Chemical Space X, Chemical Space Y"));
      await session.step(31, "And the \"fields\" and \"header labels\" readings of forms viewer should be the same", () => readingsEqual(page, "fields", "header labels", el("forms viewer")));
      await session.step(32, "And forms viewer should report no error", () => reportsNoError(page, el("forms viewer")));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A molecule column is drawn on a canvas and a plain string is not", async () => {
      await session.step(37, "When user sets \"fieldsColumnNames\" property of forms viewer to \"Structure, Core, Primary Series Name\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "Structure, Core, Primary Series Name"));
      await session.step(38, "Then the \"field kind of Structure\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "field kind of Structure", el("forms viewer"), "canvas"));
      await session.step(39, "And the \"field kind of Core\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "field kind of Core", el("forms viewer"), "canvas"));
      await session.step(40, "And the \"field kind of Primary Series Name\" reading of forms viewer should be \"input\"", () => readingReads(page, "field kind of Primary Series Name", el("forms viewer"), "input"));
      await session.step(41, "And the \"Structure of card 1\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "Structure of card 1", el("forms viewer"), "canvas"));
      await session.step(42, "And the \"Primary Series Name of card 1\" reading of forms viewer should be \"Pyrrolidines\"", () => readingReads(page, "Primary Series Name of card 1", el("forms viewer"), "Pyrrolidines"));
      await session.step(43, "And forms viewer should have a \"field Structure of card 1\" area", () => hasArea(page, el("forms viewer"), "field Structure of card 1"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Renderer Size grows the canvas from small through normal to large", async () => {
      await session.step(47, "Then \"rendererSize\" property of forms viewer should be \"small\"", () => propertyShouldBe(page, "rendererSize", el("forms viewer"), "small"));
      await session.step(48, "And the \"width of Structure of card 1\" reading of forms viewer should be at least 1", () => readingAtLeast(page, "width of Structure of card 1", el("forms viewer"), 1));
      await session.step(49, "When user sets \"rendererSize\" property of forms viewer to \"normal\"", () => setProperty(page, "rendererSize", el("forms viewer"), "normal"));
      await session.step(50, "Then the \"width of Structure of card 1\" reading of forms viewer should be higher than before", () => readingHigher(page, "width of Structure of card 1", el("forms viewer")));
      await session.step(51, "And the \"height of Structure of card 1\" reading of forms viewer should be higher than before", () => readingHigher(page, "height of Structure of card 1", el("forms viewer")));
      await session.step(52, "When user sets \"rendererSize\" property of forms viewer to \"large\"", () => setProperty(page, "rendererSize", el("forms viewer"), "large"));
      await session.step(53, "Then the \"width of Structure of card 1\" reading of forms viewer should be higher than before", () => readingHigher(page, "width of Structure of card 1", el("forms viewer")));
      await session.step(54, "And the \"height of Structure of card 1\" reading of forms viewer should be higher than before", () => readingHigher(page, "height of Structure of card 1", el("forms viewer")));
      await session.step(55, "And the \"field kind of Structure\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "field kind of Structure", el("forms viewer"), "canvas"));
      await session.step(56, "When user sets \"rendererSize\" property of forms viewer to \"small\"", () => setProperty(page, "rendererSize", el("forms viewer"), "small"));
      await session.step(57, "Then the \"width of Structure of card 1\" reading of forms viewer should be lower than before", () => readingLower(page, "width of Structure of card 1", el("forms viewer")));
      await session.step(58, "And the \"height of Structure of card 1\" reading of forms viewer should be lower than before", () => readingLower(page, "height of Structure of card 1", el("forms viewer")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every card of a selection carries its own molecule canvas, and a filter keeps them", async () => {
      await session.step(66, "When user sets \"fieldsColumnNames\" property of forms viewer to \"Structure, Core, Primary Series Name\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "Structure, Core, Primary Series Name"));
      await session.step(67, "And user sets \"rendererSize\" property of forms viewer to \"normal\"", () => setProperty(page, "rendererSize", el("forms viewer"), "normal"));
      await session.step(68, "And user selects the first 3 rows", () => selectFirstRows(page, 3));
      await session.step(69, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(70, "And the record cards of forms viewer should show rows \"1, 2, 3\"", () => recordCardRows(page, el("forms viewer"), "1, 2, 3"));
      await session.step(71, "And the record cards of forms viewer should be exactly the selected rows that pass the filter", () => recordCardsAreSelection(page, el("forms viewer")));
      await session.step(72, "And every record card of forms viewer should show \"canvas\" in \"Structure\"", () => everyRecordCardShows(page, el("forms viewer"), "canvas", "Structure"));
      await session.step(73, "And every record card of forms viewer should show \"canvas\" in \"Core\"", () => everyRecordCardShows(page, el("forms viewer"), "canvas", "Core"));
      await session.step(74, "When user adds a categorical filter on \"Primary Series Name\" keeping \"Pyrrolidines\"", () => addCategoricalFilter(page, "Primary Series Name", "Pyrrolidines"));
      await session.step(75, "Then fewer than 100 rows should pass the filter", () => filterPassesFewer(page, 100));
      await session.step(76, "And the record cards of forms viewer should be exactly the selected rows that pass the filter", () => recordCardsAreSelection(page, el("forms viewer")));
      await session.step(77, "And every record card of forms viewer should show \"canvas\" in \"Structure\"", () => everyRecordCardShows(page, el("forms viewer"), "canvas", "Structure"));
      await session.step(78, "And every record card of forms viewer should show \"Pyrrolidines\" in \"Primary Series Name\"", () => everyRecordCardShows(page, el("forms viewer"), "Pyrrolidines", "Primary Series Name"));
      await session.step(79, "When user hovers over \"Primary Series Name\" filter card", () => hoverOver(page, el("\"Primary Series Name\" filter card")));
      await session.step(80, "And user clicks on close of \"Primary Series Name\" filter card", () => clickOn(page, el("close of \"Primary Series Name\" filter card")));
      await session.step(81, "And user clears the row selection", () => clearSelection(page));
      await session.step(82, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A fit column promotes Renderer Size to normal without anyone setting it", async () => {
      await session.step(86, "Given user closes all views", () => closeAllViews(page));
      await session.step(87, "And user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(88, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["fieldsColumnNames","smiles, multiple prefit"]]));
      await session.step(90, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(91, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(92, "And \"multiple prefit\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple prefit", "fit"));
      await session.step(93, "And \"rendererSize\" property of forms viewer should be \"normal\"", () => propertyShouldBe(page, "rendererSize", el("forms viewer"), "normal"));
      await session.step(94, "And the \"field kind of multiple prefit\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "field kind of multiple prefit", el("forms viewer"), "canvas"));
      await session.step(95, "And the \"field kind of smiles\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "field kind of smiles", el("forms viewer"), "canvas"));
      await session.step(96, "And the \"multiple prefit of card 1\" reading of forms viewer should be \"canvas\"", () => readingReads(page, "multiple prefit of card 1", el("forms viewer"), "canvas"));
      await session.step(97, "And the \"width of multiple prefit of card 1\" reading of forms viewer should be 200", () => readingIs(page, "width of multiple prefit of card 1", el("forms viewer"), 200));
      await session.step(98, "And the \"height of multiple prefit of card 1\" reading of forms viewer should be 100", () => readingIs(page, "height of multiple prefit of card 1", el("forms viewer"), 100));
      await session.step(99, "And forms viewer should report no error", () => reportsNoError(page, el("forms viewer")));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
