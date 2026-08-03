
//name: #{NAME}
//description: Browse and edit #{DOMAIN_TABLE} rows
//meta.role: app
//input: string path {meta.url: true; optional: true}
//output: view result
export function #{NAME}(_path) {
  // The declared `path` input is what makes the app URL-addressable; the view
  // reads the deep link itself, so the parameter is unused.
  return new DomainAppView(grok.dapi.domains.table('#{DOMAIN_TABLE}'));
}
