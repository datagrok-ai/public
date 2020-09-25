
//name: #{NAME}
//description: Creates a #{NAME} view
//tags: view
//input: map params
//input: string path
//output: view result
export function #{NAME_LOWERCASE}View(params = null, path = '') {
    return new #{NAME_TITLECASE}View(params, path);
}
