// GENERATED CODE - DO NOT MODIFY BY HAND
// Generated from @ExportApi() annotated server plugins
// See core/server/datlas/docs/code-gen.md

let _root: string = '';
let _token: string = '';

export function init(root: string, token: string) {
  _root = root;
  _token = token;
}

async function _fetch(url: string, options: RequestInit = {}): Promise<any> {
  options.headers = { ...options.headers as any, 'Authorization': _token };
  const res = await fetch(`${_root}${url}`, options);
  return res.json();
}

export namespace dapi2 {
  export namespace chats {
  }

}
