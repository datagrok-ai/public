/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/fitting.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.model.bioreactor]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBeSwitchedOn, shouldHaveValue, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {nestedTableRows, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {loadTable, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Fitting a model to data", () => {
  const session = feature(test, "features/fitting.feature", import.meta.url);
  test("Fitting a model to data", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.model.bioreactor"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens the \"Bioreactor\" model of the Diff Studio library", () => openLibraryModel(page, "Bioreactor"));
    await run.scenario("Fit opens a view of its own", async () => {
      await session.step(24, "When user clicks on Fit ribbon item", () => clickOn(page, el("Fit ribbon item")));
      await session.step(25, "Then the \"Bioreactor - fitting\" view should be current", () => viewIsCurrent(page, "Bioreactor - fitting"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A parameter is varied, and its bounds appear with it", async () => {
      await session.step(29, "When user switches on \"FFox\" input", () => switchOn(page, el("\"FFox\" input")));
      await session.step(30, "Then \"FFox (min)\" input should be switched on", () => shouldBeSwitchedOn(page, el("\"FFox (min)\" input")));
      await session.step(31, "And \"FFox (min)\" input should have value \"0.15\"", () => shouldHaveValue(page, el("\"FFox (min)\" input"), "0.15"));
      await session.step(32, "And \"FFox (max)\" input should have value \"0.25\"", () => shouldHaveValue(page, el("\"FFox (max)\" input"), "0.25"));
      await session.step(33, "When user enters \"1.0\" into \"FFox (max)\" input", () => enterInto(page, "1.0", el("\"FFox (max)\" input")));
      await session.step(34, "Then \"FFox (max)\" input should have value \"1.0\"", () => shouldHaveValue(page, el("\"FFox (max)\" input"), "1.0"));
    });
    await run.scenario("The fit needs an output to aim at and a table to aim with", async () => {
      await session.step(37, "Given the \"System:AppData/DiffStudio/library/bioreactor-experiment.csv\" file is loaded as a table", () => loadTable(page, "System:AppData/DiffStudio/library/bioreactor-experiment.csv"));
      await session.step(38, "When user switches on \"Bioreactor\" input", () => switchOn(page, el("\"Bioreactor\" input")));
      await session.step(39, "And user selects \"bioreactor-experiment\" in \"Bioreactor\" input", () => selectIn(page, "bioreactor-experiment", el("\"Bioreactor\" input")));
      await session.step(40, "Then \"Bioreactor\" input should have value \"bioreactor-experiment\"", () => shouldHaveValue(page, el("\"Bioreactor\" input"), "bioreactor-experiment"));
      await session.step(41, "And \"argument\" input should have value \"t\"", () => shouldHaveValue(page, el("\"argument\" input"), "t"));
    });
    await run.scenario("Running the fit gives its result and the loss of each iteration", async () => {
      await session.step(44, "When user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
      await session.step(45, "Then the table should have a column \"RMSE by iterations\"", () => hasColumn(page, "RMSE by iterations"));
      await session.step(46, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(47, "And the \"RMSE by iterations\" table should have at least 2 rows", () => nestedTableRows(page, "RMSE by iterations", 2));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
