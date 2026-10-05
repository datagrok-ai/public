/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/integration/bio-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.cell.helm]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnUnits, distinctValues, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPanelHas, removeColumn, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, hasArea, noErrors, painted, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {walkToColumn} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Bio menu on the HELM showcase leaves the Helm renderer intact", () => {
  const session = feature(test, "features/integration/bio-menu.feature", import.meta.url);
  test("The Bio menu on the HELM showcase leaves the Helm renderer intact", {tag: ["@journey", "@realizes:helm.cell.helm"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(23, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(24, "Then the table should have 53 rows", () => rowCount(page, 53));
    await session.step(25, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await run.scenario("Composition docks a WebLogo over the showcase column", async () => {
      await session.step(28, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
      await session.step(29, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(30, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
      await session.step(31, "And \"Sequence Column Name\" property of WebLogo viewer should be \"HELM\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "HELM"));
      await session.step(32, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
      await session.step(33, "And WebLogo viewer should have a \"position 1\" area", () => hasArea(page, el("WebLogo viewer"), "position 1"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
      await session.step(35, "When user clicks on close icon of WebLogo viewer", () => clickOn(page, el("close icon of WebLogo viewer")));
      await session.step(36, "Then WebLogo viewer should be absent", () => shouldBe(page, el("WebLogo viewer"), "absent"));
    });
    await run.scenario("Sequence Space embeds the showcase column", async () => {
      await session.step(39, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(40, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(41, "And editor of Column input in \"Sequence Space\" dialog should have text \"HELM\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "HELM"));
      await session.step(42, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
      await session.step(43, "Then the \"Sequence Space\" dialog should close", () => dialogCloses(page, "Sequence Space"));
      await session.step(44, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(45, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
      await session.step(46, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
      await session.step(47, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
      await session.step(48, "And \"Embed_X_1\" column should have at least 10 distinct values", () => distinctValues(page, "Embed_X_1", 10));
      await session.step(49, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(50, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(51, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(54, "And user removes \"Embed_X_1\" column", () => removeColumn(page, "Embed_X_1"));
      await session.step(55, "And user removes \"Embed_Y_1\" column", () => removeColumn(page, "Embed_Y_1"));
      await session.step(56, "Then scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
    });
    await run.scenario("Hierarchical Clustering attaches a tree with a leaf for every row", async () => {
      await session.step(59, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(60, "When user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
      await session.step(61, "Then \"Hierarchical Clustering\" dialog should be visible", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "visible"));
      await session.step(62, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(63, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
      await session.step(64, "And the task bar should have finished \"Creating dendrogram\"", () => taskBarFinished(page, "Creating dendrogram"));
      await session.step(65, "And the \"tree leaves\" reading of grid should be 53", () => readingIs(page, "tree leaves", el("grid"), 53));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Convert Sequence Notation writes the column in the notation asked for", async () => {
      await session.step(69, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
      await session.step(70, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
      await session.step(71, "And \"Convert Sequence Notation\" dialog should contain text \"Current notation: helm\"", () => shouldContainText(page, el("\"Convert Sequence Notation\" dialog"), "Current notation: helm"));
      await session.step(72, "When user selects \"separator\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "separator", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
      await session.step(73, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
      await session.step(74, "Then 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(75, "And a new column \"separator(HELM)\" should have been added", () => newColumnNamed(page, "separator(HELM)"));
      await session.step(76, "And \"separator(HELM)\" column should have units \"separator\"", () => columnUnits(page, "separator(HELM)", "separator"));
      await session.step(77, "And the value of \"separator(HELM)\" column in row 2 should be \"A-C-D-E-F-G-H-I-K-L\"", () => valueInRow(page, "separator(HELM)", 2, "A-C-D-E-F-G-H-I-K-L"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "When user removes \"separator(HELM)\" column", () => removeColumn(page, "separator(HELM)"));
    });
    await run.scenario("Extract Region cuts a HELM region", async () => {
      await session.step(82, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
      await session.step(83, "Then \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "visible"));
      await session.step(84, "When user selects \"1\" in Start input in \"Get Sequence Region\" dialog", () => selectIn(page, "1", el("Start input in \"Get Sequence Region\" dialog")));
      await session.step(85, "And user selects \"2\" in End input in \"Get Sequence Region\" dialog", () => selectIn(page, "2", el("End input in \"Get Sequence Region\" dialog")));
      await session.step(86, "And user enters \"region 1-2\" into \"Column name\" input in \"Get Sequence Region\" dialog", () => enterInto(page, "region 1-2", el("\"Column name\" input in \"Get Sequence Region\" dialog")));
      await session.step(87, "And user clicks on OK button in \"Get Sequence Region\" dialog", () => clickOn(page, el("OK button in \"Get Sequence Region\" dialog")));
      await session.step(88, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(89, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(90, "And \"region 1-2\" column should have units \"helm\"", () => columnUnits(page, "region 1-2", "helm"));
      await session.step(91, "And the value of \"region 1-2\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "region 1-2", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(92, "And the value of \"region 1-2\" column in row 2 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "region 1-2", 2, "PEPTIDE1{A.C}$$$$"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
      await session.step(94, "When user removes \"region 1-2\" column", () => removeColumn(page, "region 1-2"));
    });
    await run.scenario("Scan Liabilities annotates the showcase column", async () => {
      await session.step(97, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
      await session.step(98, "Then \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "visible"));
      await session.step(99, "When user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
      await session.step(100, "Then a new column \"~HELM_annotations\" should have been added", () => newColumnNamed(page, "~HELM_annotations"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Similarity Search lists the neighbours of a HELM row", async () => {
      await session.step(104, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(105, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(106, "And \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
      await session.step(107, "And the \"source column\" reading of \"Sequence Similarity Search\" viewer should be \"HELM\"", () => readingReads(page, "source column", el("\"Sequence Similarity Search\" viewer"), "HELM"));
      await session.step(108, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 11", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 11));
      await session.step(109, "And \"Sequence Similarity Search\" viewer should be painted", () => painted(page, el("\"Sequence Similarity Search\" viewer")));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
      await session.step(111, "When user clicks on close icon of \"Sequence Similarity Search\" viewer", () => clickOn(page, el("close icon of \"Sequence Similarity Search\" viewer")));
      await session.step(112, "Then \"Sequence Similarity Search\" viewer should be absent", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "absent"));
    });
    await run.scenario("Diversity Search picks a varied HELM subset", async () => {
      await session.step(115, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(116, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(117, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(118, "And the \"source column\" reading of \"Sequence Diversity Search\" viewer should be \"HELM\"", () => readingReads(page, "source column", el("\"Sequence Diversity Search\" viewer"), "HELM"));
      await session.step(119, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 10", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 10));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
      await session.step(121, "When user clicks on close icon of \"Sequence Diversity Search\" viewer", () => clickOn(page, el("close icon of \"Sequence Diversity Search\" viewer")));
      await session.step(122, "Then \"Sequence Diversity Search\" viewer should be absent", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "absent"));
    });
    await run.scenario("Subsequence Search adds a filter on the HELM column", async () => {
      await session.step(125, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
      await session.step(126, "Then filters viewer should be visible", () => shouldBe(page, el("filters viewer"), "visible"));
      await session.step(127, "And the filter panel should have a filter on \"HELM\" column", () => filterPanelHas(page, "HELM"));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split to Monomers gives a Monomer column per position", async () => {
      await session.step(131, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(132, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
      await session.step(133, "And editor of Sequence input in \"Split to Monomers\" dialog should have text \"HELM\"", () => shouldHaveText(page, el("editor of Sequence input in \"Split to Monomers\" dialog"), "HELM"));
      await session.step(134, "When user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
      await session.step(135, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(136, "And \"1\" column should have semantic type \"Monomer\"", () => columnSemType(page, "1", "Monomer"));
      await session.step(137, "And the value of \"1\" column in row 2 should be \"A\"", () => valueInRow(page, "1", 2, "A"));
      await session.step(138, "And the value of \"10\" column in row 2 should be \"L\"", () => valueInRow(page, "10", 2, "L"));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After the sweep the column is still a HELM column painted by Helm", async () => {
      await session.step(142, "When user moves the current cell of grid to the \"HELM\" column", () => walkToColumn(page, el("grid"), "HELM"));
      await session.step(143, "Then the table should have 53 rows", () => rowCount(page, 53));
      await session.step(144, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
      await session.step(145, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
      await session.step(146, "And \"HELM\" column should have tag \"cell.renderer\" equal to \"helm\"", () => columnTag(page, "HELM", "cell.renderer", "helm"));
      await session.step(147, "And the \"cell type of HELM\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM", el("grid"), "helm"));
      await session.step(148, "And the \"cell 2 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 2 of HELM", el("grid"), 3));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
