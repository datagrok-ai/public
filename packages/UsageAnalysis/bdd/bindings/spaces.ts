/* Spaces as the Browse tree shows them. Everything general — the browse panel, the browse tree,
   the expanded state of a node — is platform vocabulary in the library; what stays here is the
   space view's own gallery and search; the spaces a scenario leaves on the server are cleaned by
   the platform's "no space named … is on the server". */
import {element} from '@datagrok-libraries/bdd';

/* The gallery a space view shows, and its search — the same two elements in the Spaces list and
   inside a space. The cards themselves are links, so "X link in space gallery" names one without
   catching the identically-classed links of an open help pane. */
// the plain "gallery" is the platform's own element now (bindings/platform/elements.ts) — the same
// selector, so "X link in gallery" and "X link in space gallery" name the same cards
element('space gallery', {selector: '.grok-gallery-grid', aliases: ['space content']});
element('space search', {selector: '.grok-gallery-search-bar .ui-input-type-ahead'});
