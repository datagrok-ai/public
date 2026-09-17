/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/clustering/chem-dialog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.cp.hier-clustering-chem-dialog-end-to-end]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {noSuchReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hierarchical clustering from the Chem menu", () => {
  const session = feature(test, "features/clustering/chem-dialog.feature", import.meta.url);
  test("Hierarchical clustering from the Chem menu", {tag: ["@journey", "@realizes:dendrogram.cp.hier-clustering-chem-dialog-end-to-end", "@known-failure", "@GROK-19595"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
    await run.scenario("The dialog opens on the molecule column with every distance and linkage", async () => {
      await session.step(19, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(20, "Then \"Hierarchical Clustering\" dialog should be visible", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "visible"));
      await session.step(21, "And Table input in \"Hierarchical Clustering\" dialog should have value \"mol1K\"", () => shouldHaveValue(page, el("Table input in \"Hierarchical Clustering\" dialog"), "mol1K"));
      await session.step(22, "And editor of Features input in \"Hierarchical Clustering\" dialog should have text \"(1) molecule\"", () => shouldHaveText(page, el("editor of Features input in \"Hierarchical Clustering\" dialog"), "(1) molecule"));
      await session.step(23, "And Distance input in \"Hierarchical Clustering\" dialog should have value \"euclidean\"", () => shouldHaveValue(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "euclidean"));
      await session.step(24, "And Linkage input in \"Hierarchical Clustering\" dialog should have value \"ward\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "ward"));
      await session.step(25, "And Distance input in \"Hierarchical Clustering\" dialog should offer \"euclidean, manhattan\"", () => shouldOffer(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "euclidean, manhattan"));
      await session.step(26, "And Linkage input in \"Hierarchical Clustering\" dialog should offer \"single, complete, average, weighted, centroid, median, ward\"", () => shouldOffer(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "single, complete, average, weighted, centroid, median, ward"));
      await session.step(27, "When user clicks on CANCEL button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("CANCEL button in \"Hierarchical Clustering\" dialog")));
      await session.step(28, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Euclidean distance and ward linkage attach a tree with a leaf for every molecule", async () => {
      await session.step(32, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(33, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(34, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(35, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
      await session.step(36, "And the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(37, "And the task bar should have shown \"Creating dendrogram\"", () => taskBarShown(page, "Creating dendrogram"));
      await session.step(38, "And \"Assign Clusters\" icon should be visible", () => shouldBe(page, el("\"Assign Clusters\" icon"), "visible"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the tree, then manhattan distance and single linkage attach a tree of another height", async () => {
      await session.step(42, "When user remembers the \"tree height\" reading of grid", () => rememberReading(page, "tree height", el("grid")));
      await session.step(43, "And user clicks on \"Remove Dendrogram\" icon", () => clickOn(page, el("\"Remove Dendrogram\" icon")));
      await session.step(44, "Then \"Assign Clusters\" icon should be absent", () => shouldBe(page, el("\"Assign Clusters\" icon"), "absent"));
      await session.step(45, "And grid should not report a \"tree leaves\" reading", () => noSuchReading(page, el("grid"), "tree leaves"));
      await session.step(46, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(47, "And user selects \"manhattan\" in Distance input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "manhattan", el("Distance input in \"Hierarchical Clustering\" dialog")));
      await session.step(48, "And user selects \"single\" in Linkage input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "single", el("Linkage input in \"Hierarchical Clustering\" dialog")));
      await session.step(49, "Then Distance input in \"Hierarchical Clustering\" dialog should have value \"manhattan\"", () => shouldHaveValue(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "manhattan"));
      await session.step(50, "And Linkage input in \"Hierarchical Clustering\" dialog should have value \"single\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "single"));
      await session.step(51, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(52, "Then the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(53, "And the \"tree height\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "tree height", el("grid")));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Numeric columns with median linkage attach a tree with a leaf for every row", async () => {
      await session.step(58, "When user remembers the \"tree height\" reading of grid", () => rememberReading(page, "tree height", el("grid")));
      await session.step(59, "And user clicks on \"Remove Dendrogram\" icon", () => clickOn(page, el("\"Remove Dendrogram\" icon")));
      await session.step(60, "And user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(61, "And user clicks on editor of Features input in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("editor of Features input in \"Hierarchical Clustering\" dialog")));
      await session.step(62, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(63, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(64, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(65, "And user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(66, "Then the \"text of cell 4 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"pIC50_HIV_Integrase\"", () => readingReads(page, "text of cell 4 of __name", el("grid viewer in \"Select columns...\" dialog"), "pIC50_HIV_Integrase"));
      await session.step(67, "And the \"text of cell 5 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Q\"", () => readingReads(page, "text of cell 5 of __name", el("grid viewer in \"Select columns...\" dialog"), "Q"));
      await session.step(68, "And \"Select columns...\" dialog should contain text \"2 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "2 checked"));
      await session.step(69, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(70, "Then editor of Features input in \"Hierarchical Clustering\" dialog should contain text \"(2)\"", () => shouldContainText(page, el("editor of Features input in \"Hierarchical Clustering\" dialog"), "(2)"));
      await session.step(71, "When user selects \"median\" in Linkage input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "median", el("Linkage input in \"Hierarchical Clustering\" dialog")));
      await session.step(72, "Then Linkage input in \"Hierarchical Clustering\" dialog should have value \"median\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "median"));
      await session.step(73, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(74, "Then the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(75, "And the \"tree height\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "tree height", el("grid")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Numeric columns with centroid linkage are run", async () => {
      await session.step(79, "When user clicks on \"Remove Dendrogram\" icon", () => clickOn(page, el("\"Remove Dendrogram\" icon")));
      await session.step(80, "And user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
      await session.step(81, "And user clicks on editor of Features input in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("editor of Features input in \"Hierarchical Clustering\" dialog")));
      await session.step(82, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(83, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(84, "And user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(85, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(86, "And user selects \"centroid\" in Linkage input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "centroid", el("Linkage input in \"Hierarchical Clustering\" dialog")));
      await session.step(87, "Then Linkage input in \"Hierarchical Clustering\" dialog should have value \"centroid\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "centroid"));
      await session.step(88, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(89, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
    });
    await run.scenario("Numeric columns with centroid linkage attach a tree with a leaf for every row", async () => {
      await session.step(93, "Then the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
