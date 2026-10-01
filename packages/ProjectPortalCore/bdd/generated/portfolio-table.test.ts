/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/portfolio-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addressEndsWith, matrixFollows, openApp, tableViewers} from '../bindings/steps.js';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterIsExactlyCategory} from '@datagrok-libraries/bdd/bindings/platform/data';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The portfolio as a table, with filters and analytics", () => {
  const session = feature(test, "features/portfolio-table.feature", import.meta.url);
  test("The development portfolio as a table, filtered to one stage", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(13, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(14, "And user clicks on table button", () => clickOn(page, el("table button")));
    await session.step(15, "Then the \"Example Biologics Development\" table view should open with 7 rows", () => tableViewOpened(page, "Example Biologics Development", 7));
    await session.step(16, "And the address should end with \"/Portal/EX-BIOC/table\"", () => addressEndsWith(page, "/Portal/EX-BIOC/table"));
    await session.step(17, "And the portfolio table should show the viewers \"Programs by stage, Days in stage, Readiness for the next gate, Compounds by role\"", () => tableViewers(page, "Programs by stage, Days in stage, Readiness for the next gate, Compounds by role"));
    await session.step(18, "When user adds a categorical filter on \"stage\" keeping \"D4\"", () => addCategoricalFilter(page, "stage", "D4"));
    await session.step(19, "Then the filter should pass exactly the rows where \"stage\" is \"D4\"", () => filterIsExactlyCategory(page, "stage", "D4"));
    await session.step(20, "And the readiness matrix should show the programs that pass the filter", () => matrixFollows(page));
  });
});
