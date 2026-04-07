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
    export async function getComments(id: string, options?: {page?: number, limit?: number, order?: string}): Promise<any[]> {
      let url = `/chats/comments/${id}`;
      const params = new URLSearchParams();
      if (options?.page !== undefined) params.set('page', String(options.page));
      if (options?.limit !== undefined) params.set('limit', String(options.limit));
      if (options?.order !== undefined) params.set('order', String(options.order));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function countComments(id: string): Promise<number> {
      let url = `/chats/comments/${id}/count`;
      return _fetch(url);
    }

    export async function upVoteComment(id: string): Promise<string> {
      let url = `/chats/comments/${id}/up_vote`;
      return _fetch(url, {method: 'POST'});
    }

    export async function downVoteComment(id: string): Promise<string> {
      let url = `/chats/comments/${id}/down_vote`;
      return _fetch(url, {method: 'POST'});
    }

    export async function postComment(id: string, comment: any): Promise<any> {
      let url = `/chats/comments/${id}`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(comment)});
    }

    export async function deleteComment(id: string): Promise<boolean> {
      let url = `/chats/comments/${id}`;
      return _fetch(url, {method: 'DELETE'});
    }

    export async function deleteChat(id: string): Promise<boolean> {
      let url = `/chats/${id}`;
      return _fetch(url, {method: 'DELETE'});
    }

    export async function save(chat: any): Promise<any> {
      let url = `/chats`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(chat)});
    }

    export async function view(id: string): Promise<string> {
      let url = `/chats/${id}/view`;
      return _fetch(url, {method: 'POST'});
    }

    export async function isWatched(id: string): Promise<boolean> {
      let url = `/chats/${id}/watch`;
      return _fetch(url);
    }

    export async function watch(id: string): Promise<string> {
      let url = `/chats/${id}/watch`;
      return _fetch(url, {method: 'POST'});
    }

    export async function read(id: string): Promise<string> {
      let url = `/chats/${id}/read`;
      return _fetch(url, {method: 'POST'});
    }

    export async function unWatch(id: string): Promise<string> {
      let url = `/chats/${id}/un_watch`;
      return _fetch(url, {method: 'POST'});
    }

    export async function list(options?: {entityId?: string, userId?: string, onlyPrivate?: boolean, include?: string, text?: string, id?: string, page?: number, limit?: number, order?: string}): Promise<any[]> {
      let url = `/chats`;
      const params = new URLSearchParams();
      if (options?.entityId !== undefined) params.set('entityId', String(options.entityId));
      if (options?.userId !== undefined) params.set('userId', String(options.userId));
      if (options?.onlyPrivate !== undefined) params.set('onlyPrivate', String(options.onlyPrivate));
      if (options?.include !== undefined) params.set('include', String(options.include));
      if (options?.text !== undefined) params.set('text', String(options.text));
      if (options?.id !== undefined) params.set('id', String(options.id));
      if (options?.page !== undefined) params.set('page', String(options.page));
      if (options?.limit !== undefined) params.set('limit', String(options.limit));
      if (options?.order !== undefined) params.set('order', String(options.order));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function getChatsWithGroups(options?: {ids?: string}): Promise<any[]> {
      let url = `/chats/with_groups`;
      const params = new URLSearchParams();
      if (options?.ids !== undefined) params.set('ids', String(options.ids));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function count(options?: {entityId?: string, userId?: string, onlyPrivate?: boolean, include?: string, text?: string, id?: string}): Promise<number> {
      let url = `/chats/count`;
      const params = new URLSearchParams();
      if (options?.entityId !== undefined) params.set('entityId', String(options.entityId));
      if (options?.userId !== undefined) params.set('userId', String(options.userId));
      if (options?.onlyPrivate !== undefined) params.set('onlyPrivate', String(options.onlyPrivate));
      if (options?.include !== undefined) params.set('include', String(options.include));
      if (options?.text !== undefined) params.set('text', String(options.text));
      if (options?.id !== undefined) params.set('id', String(options.id));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function countUnread(options?: {entityId?: string, userId?: string, include?: string, text?: string, id?: string}): Promise<number> {
      let url = `/chats/count_unread`;
      const params = new URLSearchParams();
      if (options?.entityId !== undefined) params.set('entityId', String(options.entityId));
      if (options?.userId !== undefined) params.set('userId', String(options.userId));
      if (options?.include !== undefined) params.set('include', String(options.include));
      if (options?.text !== undefined) params.set('text', String(options.text));
      if (options?.id !== undefined) params.set('id', String(options.id));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function get(id: string, options?: {include?: string}): Promise<any> {
      let url = `/chats/${id}`;
      const params = new URLSearchParams();
      if (options?.include !== undefined) params.set('include', String(options.include));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function exportChats(): Promise<any> {
      let url = `/chats/export`;
      return _fetch(url, {method: 'POST'});
    }

  }

}
