/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/pca.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.pca]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, shouldBe, shouldContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, distinctValues, everyValueBetween, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Principal component analysis and Reduce Dimensionality over cars", () => {
  const session = feature(test, "features/analyze/pca.feature", import.meta.url);
  test("Principal component analysis and Reduce Dimensionality over cars", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.pca"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(21, "Then table \"cars\" should have 30 rows", () => tableRows(page, "cars", 30));
    await run.scenario("Three components over every column add PC1 to PC3", async () => {
      await session.step(24, "When user picks \"ML > Analyze > PCA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PCA..."));
      await session.step(25, "Then \"PCA\" dialog should be visible", () => shouldBe(page, el("\"PCA\" dialog"), "visible"));
      await session.step(26, "When user clicks on editor of Features input in \"PCA\" dialog", () => clickOn(page, el("editor of Features input in \"PCA\" dialog")));
      await session.step(27, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(28, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(29, "Then \"Select columns...\" dialog should contain text \"16 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "16 checked"));
      await session.step(30, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(31, "Then editor of Features input in \"PCA\" dialog should contain text \"(16)\"", () => shouldContainText(page, el("editor of Features input in \"PCA\" dialog"), "(16)"));
      await session.step(32, "When user enters \"3\" into Components input in \"PCA\" dialog", () => enterInto(page, "3", el("Components input in \"PCA\" dialog")));
      await session.step(33, "And user unchecks Center input in \"PCA\" dialog", () => uncheck(page, el("Center input in \"PCA\" dialog")));
      await session.step(34, "And user unchecks Scale input in \"PCA\" dialog", () => uncheck(page, el("Scale input in \"PCA\" dialog")));
      await session.step(35, "Then Center input in \"PCA\" dialog should be unchecked", () => shouldBe(page, el("Center input in \"PCA\" dialog"), "unchecked"));
      await session.step(36, "And Scale input in \"PCA\" dialog should be unchecked", () => shouldBe(page, el("Scale input in \"PCA\" dialog"), "unchecked"));
      await session.step(37, "When user clicks on OK button in \"PCA\" dialog", () => clickOn(page, el("OK button in \"PCA\" dialog")));
      await session.step(38, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(39, "And \"PCA\" dialog should be hidden", () => shouldBe(page, el("\"PCA\" dialog"), "hidden"));
      await session.step(40, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(41, "And the table should have a column \"PC1\"", () => hasColumn(page, "PC1"));
      await session.step(42, "And the table should have a column \"PC2\"", () => hasColumn(page, "PC2"));
      await session.step(43, "And the table should have a column \"PC3\"", () => hasColumn(page, "PC3"));
      await session.step(44, "And \"PC1\" column should have no missing values", () => columnComplete(page, "PC1"));
      await session.step(45, "And \"PC2\" column should have no missing values", () => columnComplete(page, "PC2"));
      await session.step(46, "And \"PC3\" column should have no missing values", () => columnComplete(page, "PC3"));
      await session.step(47, "And \"PC1\" column should have at least 2 distinct values", () => distinctValues(page, "PC1", 2));
      await session.step(48, "And \"PC2\" column should have at least 2 distinct values", () => distinctValues(page, "PC2", 2));
      await session.step(49, "And \"PC3\" column should have at least 2 distinct values", () => distinctValues(page, "PC3", 2));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Center and Scale add a second set of components beside the first", async () => {
      await session.step(54, "When user picks \"ML > Analyze > PCA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PCA..."));
      await session.step(55, "Then \"PCA\" dialog should be visible", () => shouldBe(page, el("\"PCA\" dialog"), "visible"));
      await session.step(56, "When user clicks on editor of Features input in \"PCA\" dialog", () => clickOn(page, el("editor of Features input in \"PCA\" dialog")));
      await session.step(57, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(58, "Then \"Select columns...\" dialog should contain text \"19 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "19 checked"));
      await session.step(59, "When user types \"PC\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "PC", el("Search input in \"Select columns...\" dialog")));
      await session.step(60, "Then the \"rows shown\" reading of grid viewer in \"Select columns...\" dialog should be 3", () => readingIs(page, "rows shown", el("grid viewer in \"Select columns...\" dialog"), 3));
      await session.step(61, "And the \"text of cell 17 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"PC1\"", () => readingReads(page, "text of cell 17 of __name", el("grid viewer in \"Select columns...\" dialog"), "PC1"));
      await session.step(62, "And the \"text of cell 18 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"PC2\"", () => readingReads(page, "text of cell 18 of __name", el("grid viewer in \"Select columns...\" dialog"), "PC2"));
      await session.step(63, "And the \"text of cell 19 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"PC3\"", () => readingReads(page, "text of cell 19 of __name", el("grid viewer in \"Select columns...\" dialog"), "PC3"));
      await session.step(64, "When user clicks on the \"cell 17 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 17 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(65, "And user clicks on the \"cell 18 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 18 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(66, "And user clicks on the \"cell 19 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 19 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(67, "Then the \"text of cell 17 of x\" reading of grid viewer in \"Select columns...\" dialog should be \"false\"", () => readingReads(page, "text of cell 17 of x", el("grid viewer in \"Select columns...\" dialog"), "false"));
      await session.step(68, "And the \"text of cell 18 of x\" reading of grid viewer in \"Select columns...\" dialog should be \"false\"", () => readingReads(page, "text of cell 18 of x", el("grid viewer in \"Select columns...\" dialog"), "false"));
      await session.step(69, "And the \"text of cell 19 of x\" reading of grid viewer in \"Select columns...\" dialog should be \"false\"", () => readingReads(page, "text of cell 19 of x", el("grid viewer in \"Select columns...\" dialog"), "false"));
      await session.step(70, "And \"Select columns...\" dialog should contain text \"16 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "16 checked"));
      await session.step(71, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(72, "Then editor of Features input in \"PCA\" dialog should contain text \"(16)\"", () => shouldContainText(page, el("editor of Features input in \"PCA\" dialog"), "(16)"));
      await session.step(73, "When user enters \"3\" into Components input in \"PCA\" dialog", () => enterInto(page, "3", el("Components input in \"PCA\" dialog")));
      await session.step(74, "And user checks Center input in \"PCA\" dialog", () => check(page, el("Center input in \"PCA\" dialog")));
      await session.step(75, "And user checks Scale input in \"PCA\" dialog", () => check(page, el("Scale input in \"PCA\" dialog")));
      await session.step(76, "Then Center input in \"PCA\" dialog should be checked", () => shouldBe(page, el("Center input in \"PCA\" dialog"), "checked"));
      await session.step(77, "And Scale input in \"PCA\" dialog should be checked", () => shouldBe(page, el("Scale input in \"PCA\" dialog"), "checked"));
      await session.step(78, "When user clicks on OK button in \"PCA\" dialog", () => clickOn(page, el("OK button in \"PCA\" dialog")));
      await session.step(79, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(80, "And \"PCA\" dialog should be hidden", () => shouldBe(page, el("\"PCA\" dialog"), "hidden"));
      await session.step(81, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(82, "And the table should have a column \"PC1 (2)\"", () => hasColumn(page, "PC1 (2)"));
      await session.step(83, "And the table should have a column \"PC2 (2)\"", () => hasColumn(page, "PC2 (2)"));
      await session.step(84, "And the table should have a column \"PC3 (2)\"", () => hasColumn(page, "PC3 (2)"));
      await session.step(85, "And \"PC1 (2)\" column should have no missing values", () => columnComplete(page, "PC1 (2)"));
      await session.step(86, "And \"PC2 (2)\" column should have no missing values", () => columnComplete(page, "PC2 (2)"));
      await session.step(87, "And \"PC3 (2)\" column should have no missing values", () => columnComplete(page, "PC3 (2)"));
      await session.step(88, "And \"PC1 (2)\" column should have at least 2 distinct values", () => distinctValues(page, "PC1 (2)", 2));
      await session.step(89, "And \"PC2 (2)\" column should have at least 2 distinct values", () => distinctValues(page, "PC2 (2)", 2));
      await session.step(90, "And \"PC3 (2)\" column should have at least 2 distinct values", () => distinctValues(page, "PC3 (2)", 2));
      await session.step(91, "And every value of \"PC1 (2)\" column should lie between -6 and 6", () => everyValueBetween(page, "PC1 (2)", -6, 6));
      await session.step(92, "And every value of \"PC2 (2)\" column should lie between -6 and 6", () => everyValueBetween(page, "PC2 (2)", -6, 6));
      await session.step(93, "And every value of \"PC3 (2)\" column should lie between -6 and 6", () => everyValueBetween(page, "PC3 (2)", -6, 6));
      await session.step(94, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reduce Dimensionality on three columns adds the embedding and its clusters, and plots them", async () => {
      await session.step(98, "When user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."));
      await session.step(99, "Then \"Reduce Dimensionality\" dialog should be visible", () => shouldBe(page, el("\"Reduce Dimensionality\" dialog"), "visible"));
      await session.step(100, "When user clicks on Columns input in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("Columns input in \"Reduce Dimensionality\" dialog")));
      await session.step(101, "And user toggles the \"wheel.base\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "wheel.base", el("\"Select columns...\" dialog")));
      await session.step(102, "And user toggles the \"length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "length", el("\"Select columns...\" dialog")));
      await session.step(103, "And user toggles the \"width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "width", el("\"Select columns...\" dialog")));
      await session.step(104, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(105, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(106, "When user clicks on OK button in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("OK button in \"Reduce Dimensionality\" dialog")));
      await session.step(107, "Then the task bar should have finished \"dimensionality reduction\"", () => taskBarFinished(page, "dimensionality reduction"));
      await session.step(108, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(109, "And a new column \"Cluster (DBSCAN)\" should have been added", () => newColumnNamed(page, "Cluster (DBSCAN)"));
      await session.step(110, "And \"title\" property of scatter plot viewer should be \"UMAP (wheel.base, length, width)\"", () => propertyShouldBe(page, "title", el("scatter plot viewer"), "UMAP (wheel.base, length, width)"));
      await session.step(111, "And \"xColumnName\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(112, "And \"yColumnName\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(113, "And \"colorColumnName\" property of scatter plot viewer should be \"Cluster (DBSCAN)\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Cluster (DBSCAN)"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
