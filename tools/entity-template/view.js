
//name: #{NAME}
//description: Creates #{NAME} view
//meta.role: view
//input: map params
//input: string path
//output: view result
export function _#{NAME}(params = null, path = '') {
  return new #{NAME}(params, path);
}
