/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/catalog-looks.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.image-input, u2.filter-builder]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../bindings/demo.js';
import {rightOf, topAligned} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Catalog looks", () => {
  const session = feature(test, "features/demo/catalog-looks.feature", import.meta.url);
  test("The image preview sits beside its caption", {tag: ["@demo", "@realizes:u2.image-input", "@realizes:u2.filter-builder"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Catalog\" demo page", () => openDemoPage(page, "Catalog"));
    enter(page, "U2 Demo");
    await session.step(11, "Then editor of \"image-input-default\" image input should be right of label of \"image-input-default\" image input", () => rightOf(page, el("editor of \"image-input-default\" image input"), el("label of \"image-input-default\" image input")));
    await session.step(12, "And editor of \"image-input-default\" image input should be top-aligned with label of \"image-input-default\" image input", () => topAligned(page, el("editor of \"image-input-default\" image input"), el("label of \"image-input-default\" image input")));
  });
  test("A read-only filter builder offers nothing that edits it", {tag: ["@demo", "@realizes:u2.image-input", "@realizes:u2.filter-builder"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Catalog\" demo page", () => openDemoPage(page, "Catalog"));
    enter(page, "U2 Demo");
    await session.step(15, "Then Advanced link in \"filter-builder-readOnly\" filter builder should be hidden", () => shouldBe(page, el("Advanced link in \"filter-builder-readOnly\" filter builder"), "hidden"));
    await session.step(16, "And \"Add condition\" button in \"filter-builder-readOnly\" filter builder should be hidden", () => shouldBe(page, el("\"Add condition\" button in \"filter-builder-readOnly\" filter builder"), "hidden"));
    await session.step(17, "And Advanced link in \"filter-builder-default\" filter builder should be visible", () => shouldBe(page, el("Advanced link in \"filter-builder-default\" filter builder"), "visible"));
    await session.step(18, "And \"Add condition\" button in \"filter-builder-default\" filter builder should be visible", () => shouldBe(page, el("\"Add condition\" button in \"filter-builder-default\" filter builder"), "visible"));
  });
  test("A blank required filter builder shows its error", {tag: ["@demo", "@realizes:u2.image-input", "@realizes:u2.filter-builder"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the \"Catalog\" demo page", () => openDemoPage(page, "Catalog"));
    enter(page, "U2 Demo");
    await session.step(21, "Then error of \"filter-builder-invalid\" filter builder should be visible", () => shouldBe(page, el("error of \"filter-builder-invalid\" filter builder"), "visible"));
    await session.step(22, "And error of \"filter-builder-invalid\" filter builder should have text \"Value can't be empty\"", () => shouldHaveText(page, el("error of \"filter-builder-invalid\" filter builder"), "Value can't be empty"));
  });
});
