/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/app-address-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("App address edge cases: an address naming no app", () => {
  const session = feature(test, "features/apps/app-address-edge-cases.feature", import.meta.url);
  test("An address naming no app says \"Application not found\" and opens the Apps view", {tag: ["@apps", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(13, "When user opens the address \"/apps/BddNoSuchApp{time}\"", () => openAddress(page, session.text("/apps/BddNoSuchApp{time}")));
    await session.step(14, "Then an error balloon containing \"Application not found\" should have been shown", () => errorBalloonText(page, "Application not found"));
    await session.step(15, "And the \"Apps\" view should be current", () => viewIsCurrent(page, "Apps"));
  });
});
