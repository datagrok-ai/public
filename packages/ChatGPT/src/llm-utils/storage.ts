/* eslint-disable max-len */
/* Create a new file: conversation-storage.ts */

import {OpenAI} from 'openai';
import {MessageType, UIMessage} from './panel';

export interface StoredConversation<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam, TMeta = any> {
  id: string;
  timestamp: number;
  initialPrompt: string;
  messages: T[];
  uiMessages: UIMessage[];
  meta?: TMeta;
}

export interface StoredConversationWithContext<T extends MessageType = OpenAI.Chat.ChatCompletionMessageParam, TMeta = any>
  extends StoredConversation<T, TMeta> {
  contextId: string; // can be connection ID, project ID, or whatever context is appropriate
}

export class ConversationStorage {
  private static DB_NAME = 'chatgpt-conversations';
  private static DB_VERSION = 1;
  private static STORE_NAME = 'conversations';
  private static MAX_CONVERSATIONS = 20;

  private static async openDB(): Promise<IDBDatabase> {
    return new Promise((resolve, reject) => {
      const request = indexedDB.open(this.DB_NAME, this.DB_VERSION);

      request.onerror = () => reject(request.error);
      request.onsuccess = () => resolve(request.result);

      request.onupgradeneeded = (event) => {
        const db = (event.target as IDBOpenDBRequest).result;

        if (!db.objectStoreNames.contains(this.STORE_NAME)) {
          const store = db.createObjectStore(this.STORE_NAME, {keyPath: 'id'});
          store.createIndex('timestamp', 'timestamp', {unique: false});
          store.createIndex('contextId', 'contextId', {unique: false});
        }
      };
    });
  }

  static async saveConversation<TMeta = any, MessageT extends MessageType = OpenAI.Chat.ChatCompletionMessageParam>(
    messages: MessageT[],
    uiMessages: UIMessage[],
    initialPrompt: string,
    contextId?: string,
    meta?: TMeta
  ): Promise<string> {
    const db = await this.openDB();
    const conversationId = `${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;

    const conversation: StoredConversation<MessageT> = {
      id: conversationId,
      timestamp: Date.now(),
      initialPrompt,
      messages,
      uiMessages,
      ...meta ? {meta} : {},
    };

    const conversationWithContext: StoredConversationWithContext<MessageT> | StoredConversation<MessageT> = contextId ?
      {...conversation, contextId} :
      conversation;

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readwrite');
      const store = transaction.objectStore(this.STORE_NAME);

      const request = store.add(conversationWithContext);
      request.onsuccess = () => {
        this.pruneOldConversations(db);
        resolve(conversationId);
      };
      request.onerror = () => reject(request.error);
    });
  }

  static async updateConversation<MessageT extends MessageType = OpenAI.Chat.ChatCompletionMessageParam>(
    conversationId: string,
    messages: MessageT[],
    uiMessages: UIMessage[],
    meta?: any
  ): Promise<void> {
    const db = await this.openDB();

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readwrite');
      const store = transaction.objectStore(this.STORE_NAME);

      const getRequest = store.get(conversationId);
      getRequest.onsuccess = () => {
        const conversation = getRequest.result as StoredConversation<MessageT>;
        if (conversation) {
          conversation.messages = messages;
          conversation.uiMessages = uiMessages;
          conversation.timestamp = Date.now(); // Update timestamp
          if (meta != undefined)
            conversation.meta = meta;
          const putRequest = store.put(conversation);
          putRequest.onsuccess = () => resolve();
          putRequest.onerror = () => reject(putRequest.error);
        } else
          reject(new Error('Conversation not found'));
      };
      getRequest.onerror = () => reject(getRequest.error);
    });
  }

  static async getConversation<MessageT extends MessageType = OpenAI.Chat.ChatCompletionMessageParam>(conversationId: string): Promise<StoredConversation<MessageT> | null> {
    const db = await this.openDB();

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readonly');
      const store = transaction.objectStore(this.STORE_NAME);
      const request = store.get(conversationId);

      request.onsuccess = () => resolve(request.result || null);
      request.onerror = () => reject(request.error);
    });
  }

  static async listConversations<MessageT extends MessageType = OpenAI.Chat.ChatCompletionMessageParam>(
    contextId?: string,
    limit: number = 50
  ): Promise<StoredConversation<MessageT>[]> {
    const db = await this.openDB();

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readonly');
      const store = transaction.objectStore(this.STORE_NAME);

      let request: IDBRequest;

      if (contextId) {
        const index = store.index('contextId');
        request = index.getAll(contextId);
      } else
        request = store.getAll();


      request.onsuccess = () => {
        const conversations = request.result as StoredConversation<MessageT>[];
        const sorted = conversations
          .sort((a, b) => b.timestamp - a.timestamp)
          .slice(0, limit);
        resolve(sorted);
      };
      request.onerror = () => reject(request.error);
    });
  }

  static async deleteConversation(conversationId: string): Promise<void> {
    const db = await this.openDB();

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readwrite');
      const store = transaction.objectStore(this.STORE_NAME);
      const request = store.delete(conversationId);

      request.onsuccess = () => resolve();
      request.onerror = () => reject(request.error);
    });
  }

  private static async pruneOldConversations(db: IDBDatabase): Promise<void> {
    const transaction = db.transaction([this.STORE_NAME], 'readwrite');
    const store = transaction.objectStore(this.STORE_NAME);
    const index = store.index('timestamp');

    const request = index.openCursor(null, 'prev');
    let count = 0;

    request.onsuccess = (event) => {
      const cursor = (event.target as IDBRequest).result;
      if (cursor) {
        count++;
        if (count > this.MAX_CONVERSATIONS)
          cursor.delete();

        cursor.continue();
      }
    };
  }

  static async clearAll(): Promise<void> {
    const db = await this.openDB();

    return new Promise((resolve, reject) => {
      const transaction = db.transaction([this.STORE_NAME], 'readwrite');
      const store = transaction.objectStore(this.STORE_NAME);
      const request = store.clear();

      request.onsuccess = () => resolve();
      request.onerror = () => reject(request.error);
    });
  }
}
