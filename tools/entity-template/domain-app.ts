
//name: #{NAME}
//description: Browse and edit #{DOMAIN_TABLE} rows
//meta.role: app
//input: string path {meta.url: true; optional: true}
//output: view result
export async function #{NAME}(): Promise<DG.ViewBase> {
  // The declared `path` input is what makes the app URL-addressable; the page
  // reads the deep link itself (restoreFromUrl), so the function takes no parameter.
  return (await domains.table('#{DOMAIN_TABLE}')).app();
}
