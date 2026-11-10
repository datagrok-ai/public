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
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, newColumnsCount, newestMatchingDistinct, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("ML > Cluster runs DBSCAN on two columns and Markov clustering on several", () => {
  const session = feature(test, "features/analyze/clustering.feature", import.meta.url);
  test("ML > Cluster runs DBSCAN on two columns and Markov clustering on several", {tag: ["@eda", "@journey", "@realizes:public#4473"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await run.scenario("DBSCAN labels the two dense groups apart and the far point as noise", async () => {
      await session.step(20, "Given user opens a table \"points\" with:", () => openTableOf(page, "points", [["x","y"],["0","0"],["0.01","0"],["0","0.01"],["0.01","0.01"],["5","5"],["5.01","5"],["5","5.01"],["5.01","5.01"],["20","20"]]), [["x","y"],["0","0"],["0.01","0"],["0","0.01"],["0.01","0.01"],["5","5"],["5.01","5"],["5","5.01"],["5.01","5.01"],["20","20"]]);
      await session.step(31, "When user picks \"ML > Cluster > DBSCAN...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > DBSCAN..."));
      await session.step(32, "Then \"DBSCAN\" dialog should be visible", () => shouldBe(page, el("\"DBSCAN\" dialog"), "visible"));
      await session.step(33, "When user selects \"x\" in \"X Col\" input in \"DBSCAN\" dialog", () => selectIn(page, "x", el("\"X Col\" input in \"DBSCAN\" dialog")));
      await session.step(34, "And user selects \"y\" in \"Y Col\" input in \"DBSCAN\" dialog", () => selectIn(page, "y", el("\"Y Col\" input in \"DBSCAN\" dialog")));
      await session.step(35, "And user clicks on OK button in \"DBSCAN\" dialog", () => clickOn(page, el("OK button in \"DBSCAN\" dialog")));
      await session.step(36, "Then a new column \"Cluster (DBSCAN)\" should have been added", () => newColumnNamed(page, "Cluster (DBSCAN)"));
      await session.step(37, "And the newest column matching \"^Cluster \\(DBSCAN\\)\" should have 3 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(DBSCAN\\)", 3));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MCL on two number columns docks its viewer and adds the cluster columns", async () => {
      await session.step(42, "Given user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(43, "When user picks \"ML > Cluster > MCL...\" from the top menu", () => pickFromTopMenu(page, "ML > Cluster > MCL..."));
      await session.step(44, "Then \"Markov clustering\" dialog should be visible", () => shouldBe(page, el("\"Markov clustering\" dialog"), "visible"));
      await session.step(45, "When user clicks on Columns input in \"Markov clustering\" dialog", () => clickOn(page, el("Columns input in \"Markov clustering\" dialog")));
      await session.step(46, "And user toggles the \"Petal.Length\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Petal.Length", el("\"Select columns...\" dialog")));
      await session.step(47, "And user toggles the \"Petal.Width\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Petal.Width", el("\"Select columns...\" dialog")));
      await session.step(48, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Markov clustering\" dialog", () => clickOn(page, el("OK button in \"Markov clustering\" dialog")));
      await session.step(50, "Then 5 new columns should have been added", () => newColumnsCount(page, 5));
      await session.step(51, "And a new column \"Cluster (MCL)\" should have been added", () => newColumnNamed(page, "Cluster (MCL)"));
      await session.step(52, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
