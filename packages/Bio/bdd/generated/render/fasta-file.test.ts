/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/fasta-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.import.fasta, bio.export.fasta, GROK-18616]
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
import {clickOn, doubleClickOn, downloadThrough, downloadedContains, isExpanded, shouldBe, shouldHaveText, uploadDownloaded, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnSemType, columnUnits, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {renameColumn, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, openProjectWithTable, refreshBrowse, saveAsProject, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, fixtureInUsersFiles, noTableViewOpen, openTableViewsExactly, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {areaColors, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("FASTA files: open, export, reopen", () => {
  const session = feature(test, "features/render/fasta-file.feature", import.meta.url);
  test("FASTA files: open, export, reopen", {tag: ["@journey", "@realizes:bio.import.fasta", "@realizes:bio.export.fasta", "@realizes:GROK-18616"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("A FASTA file opened from the computer is a detected sequence table", async () => {
      await session.step(24, "When user uploads \"fixtures/bdd-sample.fasta\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/bdd-sample.fasta", el("\"Open local file\" icon inside browse toolbar")));
      await session.step(25, "Then the \"bdd-sample\" view should be current", () => viewIsCurrent(page, "bdd-sample"));
      await session.step(26, "And the table should have 6 rows", () => rowCount(page, 6));
      await session.step(27, "And the table should have 2 columns", () => columnCount(page, 2));
      await session.step(28, "And the value of \"description\" column in row 1 should be \"UPI0000000595:31\"", () => valueInRow(page, "description", 1, "UPI0000000595:31"));
      await session.step(29, "And the value of \"sequence\" column in row 1 should be \"MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW\"", () => valueInRow(page, "sequence", 1, "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"));
      await session.step(30, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
      await session.step(31, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
      await session.step(32, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(33, "And the \"cell 1 of sequence\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of sequence", el("grid"), 3));
      await session.step(34, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(35, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(36, "And editor of Column input in \"Sequence Space\" dialog should have text \"sequence\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "sequence"));
      await session.step(37, "When user clicks on CANCEL button in \"Sequence Space\" dialog", () => clickOn(page, el("CANCEL button in \"Sequence Space\" dialog")));
      await session.step(38, "Then \"Sequence Space\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Space\" dialog"), "hidden"));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Download as FASTA writes the table back, and the file opens as the same sequences", async () => {
      await session.step(43, "When user renames \"description\" column to \"seq id\"", () => renameColumn(page, "description", "seq id"));
      await session.step(44, "And user clicks on \"arrow to bottom\" icon", () => clickOn(page, el("\"arrow to bottom\" icon")));
      await session.step(45, "And user clicks on \"As FASTA...\" label", () => clickOn(page, el("\"As FASTA...\" label")));
      await session.step(46, "Then \"Save as FASTA\" dialog should be visible", () => shouldBe(page, el("\"Save as FASTA\" dialog"), "visible"));
      await session.step(47, "When user downloads a file through OK button in \"Save as FASTA\" dialog", () => downloadThrough(page, el("OK button in \"Save as FASTA\" dialog")));
      await session.step(48, "Then the downloaded file should contain \">UPI0000000595:31\"", () => downloadedContains(page, ">UPI0000000595:31"));
      await session.step(49, "And the downloaded file should contain \"MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW\"", () => downloadedContains(page, "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"));
      await session.step(50, "And the downloaded file should contain \">UPI000005175B:38\"", () => downloadedContains(page, ">UPI000005175B:38"));
      await session.step(51, "When user clicks on \"Browse\" view", () => clickOn(page, el("\"Browse\" view")));
      await session.step(52, "When user uploads the downloaded file through \"Open local file\" icon inside browse toolbar", () => uploadDownloaded(page, el("\"Open local file\" icon inside browse toolbar")));
      await session.step(53, "Then the table should have a column \"description\"", () => hasColumn(page, "description"));
      await session.step(54, "And the table should not have a column \"seq id\"", () => hasNoColumn(page, "seq id"));
      await session.step(55, "And the table should have 6 rows", () => rowCount(page, 6));
      await session.step(56, "And the value of \"description\" column in row 1 should be \"UPI0000000595:31\"", () => valueInRow(page, "description", 1, "UPI0000000595:31"));
      await session.step(57, "And the value of \"sequence\" column in row 1 should be \"MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW\"", () => valueInRow(page, "sequence", 1, "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"));
      await session.step(58, "And the value of \"sequence\" column in row 6 should be \"MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF\"", () => valueInRow(page, "sequence", 6, "MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF"));
      await session.step(59, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
      await session.step(60, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A FASTA file dropped onto the window is a detected sequence table", async () => {
      await session.step(65, "When user closes all views", () => closeAllViews(page));
      await session.step(66, "Then no table view is open", () => noTableViewOpen(page));
      await session.step(67, "When user drops the \"fixtures/bdd-sample.fasta\" file of the project onto browse tree", () => dropFile(page, "fixtures/bdd-sample.fasta", el("browse tree")));
      await session.step(68, "Then the open table views should be exactly \"bdd-sample\"", () => openTableViewsExactly(page, "bdd-sample"));
      await session.step(69, "And the table should have 6 rows", () => rowCount(page, 6));
      await session.step(70, "And the value of \"sequence\" column in row 1 should be \"MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW\"", () => valueInRow(page, "sequence", 1, "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"));
      await session.step(71, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
      await session.step(72, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
      await session.step(73, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(74, "And the \"cell 1 of sequence\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of sequence", el("grid"), 3));
      await session.step(75, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A FASTA-imported table keeps its type, units and renderer through a project", async () => {
      await session.step(79, "When user saves the current view as project \"bdd-bio-fasta-{run}\"", () => saveAsProject(page, session.text("bdd-bio-fasta-{run}")));
      await session.step(80, "And user closes all views", () => closeAllViews(page));
      await session.step(81, "And user opens the \"bdd-bio-fasta-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("bdd-bio-fasta-{run}")));
      await session.step(82, "Then the table should have 6 rows", () => rowCount(page, 6));
      await session.step(83, "And the value of \"sequence\" column in row 6 should be \"MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF\"", () => valueInRow(page, "sequence", 6, "MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF"));
      await session.step(84, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
      await session.step(85, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
      await session.step(86, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(87, "And the \"cell 1 of sequence\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of sequence", el("grid"), 3));
      await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A FASTA file double-clicked in Browse > Files is a detected sequence table", async () => {
      await session.step(92, "Given the \"fixtures/bdd-sample.fasta\" file of the project is in the user's files as \"bdd-fasta-{time}.fasta\"", () => fixtureInUsersFiles(page, "fixtures/bdd-sample.fasta", session.text("bdd-fasta-{time}.fasta")));
      await session.step(93, "When user closes all views", () => closeAllViews(page));
      await session.step(94, "And user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(95, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(96, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(97, "And user double-clicks on Files---My-files---bdd-fasta-{time}.fasta tree node inside browse tree", () => doubleClickOn(page, el(session.text("Files---My-files---bdd-fasta-{time}.fasta tree node inside browse tree"))));
      await session.step(98, "Then the \"bdd-fasta-{time}\" table view should open with 6 rows", () => tableViewOpened(page, session.text("bdd-fasta-{time}"), 6));
      await session.step(99, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
      await session.step(100, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
      await session.step(101, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(102, "And the \"cell 1 of sequence\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of sequence", el("grid"), 3));
      await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
