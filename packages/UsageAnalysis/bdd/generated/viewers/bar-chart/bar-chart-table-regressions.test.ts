/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-table-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, public#405]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doesNotOfferColumn, offersColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart regressions — the Table property switched in the context panel", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-table-regressions.feature", import.meta.url);
  test("After a Table switch the column properties offer the new table's columns (public#405)", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@realizes:public#405"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "Given user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(20, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(21, "Then the context panel should show \"Bar chart\"", () => contextPanelShows(page, "Bar chart"));
    await session.step(22, "When user selects \"cars\" in \"Table\" property in context panel", () => selectIn(page, "cars", el("\"Table\" property in context panel")));
    await session.step(23, "Then bar chart viewer should be bound to table \"cars\"", () => boundTable(page, el("bar chart viewer"), "cars"));
    await session.step(24, "Given \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
    await session.step(25, "Then \"Split\" property in context panel should offer the columns \"model\"", () => offersColumns(page, el("\"Split\" property in context panel"), "model"));
    await session.step(26, "And \"Split\" property in context panel should not offer the column \"RACE\"", () => doesNotOfferColumn(page, el("\"Split\" property in context panel"), "RACE"));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
});
