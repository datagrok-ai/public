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
import {check, clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Principal component analysis", () => {
  const session = feature(test, "features/analyze/pca.feature", import.meta.url);
  test("Principal component analysis", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.pca"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("Three components over every column add PC1 to PC3", async () => {
      await session.step(16, "When user picks \"ML > Analyze > PCA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PCA..."));
      await session.step(17, "Then \"PCA\" dialog should be visible", () => shouldBe(page, el("\"PCA\" dialog"), "visible"));
      await session.step(18, "When user clicks on editor of Features input in \"PCA\" dialog", () => clickOn(page, el("editor of Features input in \"PCA\" dialog")));
      await session.step(19, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(20, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(21, "Then \"Select columns...\" dialog should contain text \"16 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "16 checked"));
      await session.step(22, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(23, "Then editor of Features input in \"PCA\" dialog should contain text \"(16)\"", () => shouldContainText(page, el("editor of Features input in \"PCA\" dialog"), "(16)"));
      await session.step(24, "When user enters \"3\" into Components input in \"PCA\" dialog", () => enterInto(page, "3", el("Components input in \"PCA\" dialog")));
      await session.step(25, "And user clicks on OK button in \"PCA\" dialog", () => clickOn(page, el("OK button in \"PCA\" dialog")));
      await session.step(26, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(27, "And \"PCA\" dialog should be hidden", () => shouldBe(page, el("\"PCA\" dialog"), "hidden"));
      await session.step(28, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(29, "And the table should have a column \"PC1\"", () => hasColumn(page, "PC1"));
      await session.step(30, "And the table should have a column \"PC2\"", () => hasColumn(page, "PC2"));
      await session.step(31, "And the table should have a column \"PC3\"", () => hasColumn(page, "PC3"));
      await session.step(32, "And \"PC1\" column should have no missing values", () => columnComplete(page, "PC1"));
      await session.step(33, "And \"PC3\" column should have no missing values", () => columnComplete(page, "PC3"));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Center and Scale add a second set of components beside the first", async () => {
      await session.step(38, "When user picks \"ML > Analyze > PCA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PCA..."));
      await session.step(39, "Then \"PCA\" dialog should be visible", () => shouldBe(page, el("\"PCA\" dialog"), "visible"));
      await session.step(40, "When user clicks on editor of Features input in \"PCA\" dialog", () => clickOn(page, el("editor of Features input in \"PCA\" dialog")));
      await session.step(41, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(42, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(43, "And user enters \"3\" into Components input in \"PCA\" dialog", () => enterInto(page, "3", el("Components input in \"PCA\" dialog")));
      await session.step(44, "And user checks Center input in \"PCA\" dialog", () => check(page, el("Center input in \"PCA\" dialog")));
      await session.step(45, "And user checks Scale input in \"PCA\" dialog", () => check(page, el("Scale input in \"PCA\" dialog")));
      await session.step(46, "Then Center input in \"PCA\" dialog should be checked", () => shouldBe(page, el("Center input in \"PCA\" dialog"), "checked"));
      await session.step(47, "And Scale input in \"PCA\" dialog should be checked", () => shouldBe(page, el("Scale input in \"PCA\" dialog"), "checked"));
      await session.step(48, "When user clicks on OK button in \"PCA\" dialog", () => clickOn(page, el("OK button in \"PCA\" dialog")));
      await session.step(49, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(50, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(51, "And the table should have a column \"PC1 (2)\"", () => hasColumn(page, "PC1 (2)"));
      await session.step(52, "And the table should have a column \"PC2 (2)\"", () => hasColumn(page, "PC2 (2)"));
      await session.step(53, "And the table should have a column \"PC3 (2)\"", () => hasColumn(page, "PC3 (2)"));
      await session.step(54, "And \"PC1 (2)\" column should have no missing values", () => columnComplete(page, "PC1 (2)"));
      await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
