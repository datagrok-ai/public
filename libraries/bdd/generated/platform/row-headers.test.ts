/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/platform/row-headers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {rowsRangeSelected, selectRowsByHeaders} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Selecting rows by their row headers", () => {
  const session = feature(test, "features/platform/row-headers.feature", import.meta.url);
  test("A click and a Shift+click on the row headers select the rows between them", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(10, "When user selects rows 3 to 7 by their row headers", () => selectRowsByHeaders(page, 3, 7));
    await session.step(11, "Then rows 3 to 7 should be selected", () => rowsRangeSelected(page, 3, 7));
  });
});
