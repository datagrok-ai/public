/* The steps only the trellis plot needs. Everything else its features use is the library's
   `viewers` tier and the platform's data steps, read from the hit areas and readings the trellis
   reports (`cell F | Caucasian`, `cell body …`, `x plus`, `x range slider max handle`, `cells`,
   `cells drawn`, `distinct cell signatures`, `cell signature <cell>` and the rest — see
   `core/client/d4/lib/src/viewers/trellis_plot/CLAUDE.md`).

   What is left here is the four paging icons, whose `aria-disabled` is a state of an element and
   not of a hit area. The inner viewer's look, its type selector and the cells-wide-by-tall read
   were the matrix plot's too, and are in the library now
   (`bindings/tiers/viewers/widgets.ts`). */
import {element} from '@datagrok-libraries/bdd';

// --- the paging icons ---------------------------------------------------------------------------

/* The (+)/(-) icons that page one category row of an axis in and out. They are hit areas as well
   (`x plus`, so a click lands on them by name), but "disabled" is an attribute of the element the
   trellis writes it on. */
element('x plus icon', {selector: '[name="x-axis-icons"] .d4-trellis-plot-add-cat-icon'});
element('x minus icon', {selector: '[name="x-axis-icons"] .d4-trellis-plot-remove-cat-icon'});
element('y plus icon', {selector: '[name="y-axis-icons"] .d4-trellis-plot-add-cat-icon'});
element('y minus icon', {selector: '[name="y-axis-icons"] .d4-trellis-plot-remove-cat-icon'});
