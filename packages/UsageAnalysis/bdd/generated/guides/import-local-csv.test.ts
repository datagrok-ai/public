/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/import-local-csv.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Open a file from your computer", () => {
  const session = feature(test, "features/guides/import-local-csv.feature", import.meta.url);
  test("Open a local CSV file as a table", {tag: ["@guide", "@help:access/files"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And simple mode is off", () => simpleModeOff(page));
    await session.step(10, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(11, "When user uploads \"fixtures/browse-import.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/browse-import.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(12, "Then the \"browse-import\" view should be current", () => viewIsCurrent(page, "browse-import"));
    await session.step(13, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(14, "And the table should have 3 columns", () => columnCount(page, 3));
  });
});
