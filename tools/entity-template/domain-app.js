
//name: #{NAME}
//description: Browse and edit #{DOMAIN_TABLE} rows
//tags: app
//output: view result
export async function #{NAME}() {
  return (await domains.table('#{DOMAIN_TABLE}')).app();
}
