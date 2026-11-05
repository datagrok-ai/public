/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/clustering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [public#4473]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, newestMatchingDistinct, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("ML > Cluster runs DBSCAN on two columns and Markov clustering on several", () => {
  const session = feature(test, "features/analyze/clustering.feature", import.meta.url);
  test("ML > Cluster runs DBSCAN on two columns and Markov clustering on several", {tag: ["@eda", "@journey", "@known-failure:100", "@realizes:public#4473"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await run.scenario("DBSCAN labels the two dense groups apart and the far point as noise", async () => {
      await session.step(27, "Given user opens a table \"points\" with:", () => openTableOf(page, "points", [["x","y"],["0","0"],["0.01","0"],["0","0.01"],["0.01","0.01"],["5","5"],["5.01","5"],["5","5.01"],["5.01","5.01"],["20","20"]]), [["x","y"],["0","0"],["0.01","0"],["0","0.01"],["0.01","0.01"],["5","5"],["5.01","5"],["5","5.01"],["5.01","5.01"],["20","20"]]);
      await session.step(38, "When user picks \"ML > Cluster > DBSCAN...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > DBSCAN..."));
      await session.step(39, "Then \"DBSCAN\" dialog should be visible", () => shouldBe(page, el("\"DBSCAN\" dialog"), "visible"));
      await session.step(40, "And Epsilon input in \"DBSCAN\" dialog should have value \"0.02\"", () => shouldHaveValue(page, el("Epsilon input in \"DBSCAN\" dialog"), "0.02"));
      await session.step(41, "And \"Minimum points\" input in \"DBSCAN\" dialog should have value \"4\"", () => shouldHaveValue(page, el("\"Minimum points\" input in \"DBSCAN\" dialog"), "4"));
      await session.step(42, "When user selects \"x\" in \"X Col\" input in \"DBSCAN\" dialog", () => selectIn(page, "x", el("\"X Col\" input in \"DBSCAN\" dialog")));
      await session.step(43, "And user selects \"y\" in \"Y Col\" input in \"DBSCAN\" dialog", () => selectIn(page, "y", el("\"Y Col\" input in \"DBSCAN\" dialog")));
      await session.step(44, "And user clicks on OK button in \"DBSCAN\" dialog", () => clickOn(page, el("OK button in \"DBSCAN\" dialog")));
      await session.step(45, "Then \"DBSCAN\" dialog should be hidden", () => shouldBe(page, el("\"DBSCAN\" dialog"), "hidden"));
      await session.step(46, "And a new column \"Cluster (DBSCAN)\" should have been added", () => newColumnNamed(page, "Cluster (DBSCAN)"));
      await session.step(47, "And the newest column matching \"^Cluster \\(DBSCAN\\)\" should have 3 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(DBSCAN\\)", 3));
      await session.step(48, "And the value of \"Cluster (DBSCAN)\" column in row 1 should be \"1\"", () => valueInRow(page, "Cluster (DBSCAN)", 1, "1"));
      await session.step(49, "And the value of \"Cluster (DBSCAN)\" column in row 4 should be \"1\"", () => valueInRow(page, "Cluster (DBSCAN)", 4, "1"));
      await session.step(50, "And the value of \"Cluster (DBSCAN)\" column in row 5 should be \"2\"", () => valueInRow(page, "Cluster (DBSCAN)", 5, "2"));
      await session.step(51, "And the value of \"Cluster (DBSCAN)\" column in row 8 should be \"2\"", () => valueInRow(page, "Cluster (DBSCAN)", 8, "2"));
      await session.step(52, "And the value of \"Cluster (DBSCAN)\" column in row 9 should be \"-1\"", () => valueInRow(page, "Cluster (DBSCAN)", 9, "-1"));
      await session.step(53, "And the table should have the columns \"x, y, Cluster (DBSCAN)\"", () => columnsExactly(page, "x, y, Cluster (DBSCAN)"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MCL on a text and a number column docks its viewer and adds the cluster columns", async () => {
      await session.step(57, "When user closes all views", () => closeAllViews(page));
      await session.step(58, "And user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(59, "When user picks \"ML > Cluster > MCL...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > MCL..."));
      await session.step(60, "Then \"Markov clustering\" dialog should be visible", () => shouldBe(page, el("\"Markov clustering\" dialog"), "visible"));
      await session.step(61, "And Columns input in \"Markov clustering\" dialog should contain text \"(0)\"", () => shouldContainText(page, el("Columns input in \"Markov clustering\" dialog"), "(0)"));
      await session.step(62, "And \"Similarity Threshold\" input in \"Markov clustering\" dialog should have value \"80\"", () => shouldHaveValue(page, el("\"Similarity Threshold\" input in \"Markov clustering\" dialog"), "80"));
      await session.step(63, "And \"Inflation Factor\" input in \"Markov clustering\" dialog should have value \"2\"", () => shouldHaveValue(page, el("\"Inflation Factor\" input in \"Markov clustering\" dialog"), "2"));
      await session.step(64, "And \"Max Iterations\" input in \"Markov clustering\" dialog should have value \"5\"", () => shouldHaveValue(page, el("\"Max Iterations\" input in \"Markov clustering\" dialog"), "5"));
      await session.step(65, "And \"Min Cluster Size\" input in \"Markov clustering\" dialog should have value \"5\"", () => shouldHaveValue(page, el("\"Min Cluster Size\" input in \"Markov clustering\" dialog"), "5"));
      await session.step(66, "When user clicks on Columns input in \"Markov clustering\" dialog", () => clickOn(page, el("Columns input in \"Markov clustering\" dialog")));
      await session.step(67, "And user toggles the \"Species\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Species", el("\"Select columns...\" dialog")));
      await session.step(68, "And user toggles the \"Petal.Length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Petal.Length", el("\"Select columns...\" dialog")));
      await session.step(69, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(70, "Then Columns input in \"Markov clustering\" dialog should contain text \"(2)\"", () => shouldContainText(page, el("Columns input in \"Markov clustering\" dialog"), "(2)"));
      await session.step(71, "When user clicks on OK button in \"Markov clustering\" dialog", () => clickOn(page, el("OK button in \"Markov clustering\" dialog")));
      await session.step(72, "Then \"Markov clustering\" dialog should be hidden", () => shouldBe(page, el("\"Markov clustering\" dialog"), "hidden"));
      await session.step(73, "And 5 new columns should have been added", () => newColumnsCount(page, 5));
      await session.step(74, "And a new column \"Cluster (MCL)\" should have been added", () => newColumnNamed(page, "Cluster (MCL)"));
      await session.step(75, "And a new column \"EmbedX (MCL)\" should have been added", () => newColumnNamed(page, "EmbedX (MCL)"));
      await session.step(76, "And a new column \"EmbedY (MCL)\" should have been added", () => newColumnNamed(page, "EmbedY (MCL)"));
      await session.step(77, "And a new column \"Cluster size (MCL)\" should have been added", () => newColumnNamed(page, "Cluster size (MCL)"));
      await session.step(78, "And a new column \"Connectivity (MCL)\" should have been added", () => newColumnNamed(page, "Connectivity (MCL)"));
      await session.step(79, "And the newest column matching \"^Cluster \\(MCL\\)\" should have no missing values", () => newestMatchingFilled(page, "^Cluster \\(MCL\\)"));
      await session.step(80, "And the newest column matching \"^EmbedX \\(MCL\\)\" should have no missing values", () => newestMatchingFilled(page, "^EmbedX \\(MCL\\)"));
      await session.step(81, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(82, "And the table should have 150 rows", () => rowCount(page, 150));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MCL on four number columns adds the cluster columns", async () => {
      await session.step(87, "When user closes all views", () => closeAllViews(page));
      await session.step(88, "And user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(89, "When user picks \"ML > Cluster > MCL...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > MCL..."));
      await session.step(90, "And user clicks on Columns input in \"Markov clustering\" dialog", () => clickOn(page, el("Columns input in \"Markov clustering\" dialog")));
      await session.step(91, "And user toggles the \"Sepal.Length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Sepal.Length", el("\"Select columns...\" dialog")));
      await session.step(92, "And user toggles the \"Sepal.Width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Sepal.Width", el("\"Select columns...\" dialog")));
      await session.step(93, "And user toggles the \"Petal.Length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Petal.Length", el("\"Select columns...\" dialog")));
      await session.step(94, "And user toggles the \"Petal.Width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Petal.Width", el("\"Select columns...\" dialog")));
      await session.step(95, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(96, "Then Columns input in \"Markov clustering\" dialog should contain text \"(4)\"", () => shouldContainText(page, el("Columns input in \"Markov clustering\" dialog"), "(4)"));
      await session.step(97, "When user clicks on OK button in \"Markov clustering\" dialog", () => clickOn(page, el("OK button in \"Markov clustering\" dialog")));
      await session.step(98, "Then \"Markov clustering\" dialog should be hidden", () => shouldBe(page, el("\"Markov clustering\" dialog"), "hidden"));
      await session.step(99, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(100, "And a new column \"Cluster (MCL)\" should have been added", () => newColumnNamed(page, "Cluster (MCL)"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true, line: 100});
    run.finish();
  });
});
