// Free-text value matchers
// See the following link for examples of all matchers (numerical, string, date, bool)
// https://datagrok.ai/help/explore/data-search-patterns

let matcher = DG.ValueMatcher.numerical('> 30');

grok.shell.info(matcher.operator); // '>'
grok.shell.info(matcher.match(40)); // true
grok.shell.info(matcher.validate(20)); // "> 30" does not match "20"