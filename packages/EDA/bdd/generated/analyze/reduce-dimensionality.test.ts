/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/reduce-dimensionality.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, painted, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("ML > Reduce Dimensionality embeds the chosen columns in two dimensions and plots them", () => {
  const session = feature(test, "features/analyze/reduce-dimensionality.feature", import.meta.url);
  test("The dialog opens with no column picked, and OK on three columns adds the embedding, its plot and clusters", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(26, "When user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."));
    await session.step(27, "Then \"Reduce Dimensionality\" dialog should be visible", () => shouldBe(page, el("\"Reduce Dimensionality\" dialog"), "visible"));
    await session.step(28, "And Table input in \"Reduce Dimensionality\" dialog should have value \"cars\"", () => shouldHaveValue(page, el("Table input in \"Reduce Dimensionality\" dialog"), "cars"));
    await session.step(29, "And Columns input in \"Reduce Dimensionality\" dialog should contain text \"(0)\"", () => shouldContainText(page, el("Columns input in \"Reduce Dimensionality\" dialog"), "(0)"));
    await session.step(30, "And Method input in \"Reduce Dimensionality\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Reduce Dimensionality\" dialog"), "UMAP"));
    await session.step(31, "And Method input in \"Reduce Dimensionality\" dialog should offer \"UMAP, t-SNE\"", () => shouldOffer(page, el("Method input in \"Reduce Dimensionality\" dialog"), "UMAP, t-SNE"));
    await session.step(32, "And \"Plot Embeddings\" input in \"Reduce Dimensionality\" dialog should be checked", () => shouldBe(page, el("\"Plot Embeddings\" input in \"Reduce Dimensionality\" dialog"), "checked"));
    await session.step(33, "And Postprocessing input in \"Reduce Dimensionality\" dialog should have value \"DBSCAN clustering\"", () => shouldHaveValue(page, el("Postprocessing input in \"Reduce Dimensionality\" dialog"), "DBSCAN clustering"));
    await session.step(34, "When user clicks on Columns input in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("Columns input in \"Reduce Dimensionality\" dialog")));
    await session.step(35, "And user toggles the \"wheel.base\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "wheel.base", el("\"Select columns...\" dialog")));
    await session.step(36, "And user toggles the \"length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "length", el("\"Select columns...\" dialog")));
    await session.step(37, "And user toggles the \"width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "width", el("\"Select columns...\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(39, "Then Columns input in \"Reduce Dimensionality\" dialog should contain text \"(3)\"", () => shouldContainText(page, el("Columns input in \"Reduce Dimensionality\" dialog"), "(3)"));
    await session.step(40, "And OK button in \"Reduce Dimensionality\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Reduce Dimensionality\" dialog"), "enabled"));
    await session.step(41, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(42, "When user clicks on OK button in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("OK button in \"Reduce Dimensionality\" dialog")));
    await session.step(43, "Then \"Reduce Dimensionality\" dialog should be hidden", () => shouldBe(page, el("\"Reduce Dimensionality\" dialog"), "hidden"));
    await session.step(44, "And the task bar should have finished \"dimensionality reduction\"", () => taskBarFinished(page, "dimensionality reduction"));
    await session.step(45, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
    await session.step(46, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(47, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
    await session.step(48, "And a new column \"Cluster (DBSCAN)\" should have been added", () => newColumnNamed(page, "Cluster (DBSCAN)"));
    await session.step(49, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
    await session.step(50, "And the newest column matching \"^Embed_Y_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_Y_"));
    await session.step(51, "And the newest column matching \"^Cluster \\(DBSCAN\\)\" should have no missing values", () => newestMatchingFilled(page, "^Cluster \\(DBSCAN\\)"));
    await session.step(52, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(53, "And \"title\" property of scatter plot viewer should be \"UMAP (wheel.base, length, width)\"", () => propertyShouldBe(page, "title", el("scatter plot viewer"), "UMAP (wheel.base, length, width)"));
    await session.step(54, "And \"xColumnName\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "Embed_X_1"));
    await session.step(55, "And \"yColumnName\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "Embed_Y_1"));
    await session.step(56, "And \"colorColumnName\" property of scatter plot viewer should be \"Cluster (DBSCAN)\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Cluster (DBSCAN)"));
    await session.step(57, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
  });
  test("With no column picked OK is disabled", {tag: ["@eda", "@known-failure:64"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await knownFailure(async () => {
      await session.step(62, "When user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."));
      await session.step(63, "Then Columns input in \"Reduce Dimensionality\" dialog should contain text \"(0)\"", () => shouldContainText(page, el("Columns input in \"Reduce Dimensionality\" dialog"), "(0)"));
      await session.step(64, "And OK button in \"Reduce Dimensionality\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Reduce Dimensionality\" dialog"), "disabled"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    }, 64);
  });
  test("t-SNE with no postprocessing adds only the embedding and plots it under the method's name", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(68, "When user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."));
    await session.step(69, "And user clicks on Columns input in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("Columns input in \"Reduce Dimensionality\" dialog")));
    await session.step(70, "And user toggles the \"wheel.base\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "wheel.base", el("\"Select columns...\" dialog")));
    await session.step(71, "And user toggles the \"length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "length", el("\"Select columns...\" dialog")));
    await session.step(72, "And user toggles the \"width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "width", el("\"Select columns...\" dialog")));
    await session.step(73, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(74, "And user selects \"t-SNE\" in Method input in \"Reduce Dimensionality\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Reduce Dimensionality\" dialog")));
    await session.step(75, "And user selects \"None\" in Postprocessing input in \"Reduce Dimensionality\" dialog", () => selectIn(page, "None", el("Postprocessing input in \"Reduce Dimensionality\" dialog")));
    await session.step(76, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(77, "When user clicks on OK button in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("OK button in \"Reduce Dimensionality\" dialog")));
    await session.step(78, "Then the task bar should have finished \"dimensionality reduction\"", () => taskBarFinished(page, "dimensionality reduction"));
    await session.step(79, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
    await session.step(80, "And 2 new columns should have been added", () => newColumnsCount(page, 2));
    await session.step(81, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(82, "And \"title\" property of scatter plot viewer should be \"t-SNE (wheel.base, length, width)\"", () => propertyShouldBe(page, "title", el("scatter plot viewer"), "t-SNE (wheel.base, length, width)"));
    await session.step(83, "And the table should not have a column \"Cluster (DBSCAN)\"", () => hasNoColumn(page, "Cluster (DBSCAN)"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
  });
  test("Without Plot Embeddings the columns are added and no scatter plot is docked", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(87, "When user picks \"ML > Reduce Dimensionality...\" from the top menu", () => pickFromTopMenu(page, "ML > Reduce Dimensionality..."));
    await session.step(88, "And user clicks on Columns input in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("Columns input in \"Reduce Dimensionality\" dialog")));
    await session.step(89, "And user toggles the \"wheel.base\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "wheel.base", el("\"Select columns...\" dialog")));
    await session.step(90, "And user toggles the \"length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "length", el("\"Select columns...\" dialog")));
    await session.step(91, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(92, "And user unchecks \"Plot Embeddings\" input in \"Reduce Dimensionality\" dialog", () => uncheck(page, el("\"Plot Embeddings\" input in \"Reduce Dimensionality\" dialog")));
    await session.step(93, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(94, "When user clicks on OK button in \"Reduce Dimensionality\" dialog", () => clickOn(page, el("OK button in \"Reduce Dimensionality\" dialog")));
    await session.step(95, "Then the task bar should have finished \"dimensionality reduction\"", () => taskBarFinished(page, "dimensionality reduction"));
    await session.step(96, "And the newest column matching \"^Embed_Y_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_Y_"));
    await session.step(97, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(98, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
});
