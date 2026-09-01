import { IDartApi } from "./api/grok_api.g";
import { toJs } from "./wrappers";
import { Entity, Func } from './entities';
import { FUNC_TYPES } from './const';
import * as rxjs from 'rxjs';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export interface AIConfigBase {
    provider: string;
    configured: boolean;
    indexEntities: boolean;
}

export interface OpenAIConfig extends AIConfigBase {
    provider: 'openai';
}

export interface AzureOpenAIConfig extends AIConfigBase {
    provider: 'azure';
    apiMode: 'openai_compat' | 'legacy';
    apiVersion?: string;
    modelToDeployment?: Record<string, string>;
}

export type AIConfig = OpenAIConfig | AzureOpenAIConfig;

export interface AIChatOptions {
    engine?: string;
    model?: string;
    systemPrompt?: string;
}

export type AIAttachment =
    | {type: 'image', data: Blob}
    | {type: 'document', data: Blob, title?: string};

export interface AIRunOptions {
    signal?: AbortSignal;
    schema?: object;
    attachments?: AIAttachment[];
}

export interface AIUsage {
    inputTokens?: number;
    outputTokens?: number;
    costUsd?: number;
}

export interface AIResult {
    text: string;
    structuredOutput?: unknown;
    usage: AIUsage;
}

export type AIStream = AsyncGenerator<string, AIResult>;

export interface AITurnEvent {
    chat: AIChat;
    prompt: string;
    options?: AIRunOptions;
}

const _events = {
    onPrompt: new rxjs.Subject<AITurnEvent>(),
    onChunk: new rxjs.Subject<AITurnEvent & {chunk: string}>(),
    onResult: new rxjs.Subject<AITurnEvent & {result: AIResult}>(),
    onCancel: new rxjs.Subject<AITurnEvent>(),
    onError: new rxjs.Subject<AITurnEvent & {error: any}>(),
};

export abstract class AIChat {
    protected abstract respond(prompt: string, options?: AIRunOptions): AIStream;

    async* stream(prompt: string, options?: AIRunOptions): AIStream {
        const event: AITurnEvent = {chat: this, prompt, options};
        const s = this.respond(prompt, options);
        let settled = false;
        _events.onPrompt.next(event);
        try {
            let step = await s.next();
            while (!step.done) {
                _events.onChunk.next({...event, chunk: step.value});
                yield step.value;
                step = await s.next();
            }
            settled = true;
            _events.onResult.next({...event, result: step.value});
            return step.value;
        }
        catch (error) {
            settled = true;
            if (options?.signal?.aborted)
                _events.onCancel.next(event);
            else
                _events.onError.next({...event, error});
            throw error;
        }
        finally {
            if (!settled) {
                await s.return(undefined as any);
                _events.onCancel.next(event);
            }
        }
    }

    async prompt(prompt: string, options?: AIRunOptions): Promise<AIResult> {
        const s = this.stream(prompt, options);
        let step = await s.next();
        while (!step.done)
            step = await s.next();
        return step.value;
    }

    close(): void {}
}

export abstract class AIEngine {
    abstract id: string;
    abstract name: string;
    abstract available(): Promise<boolean>;
    models?(): Promise<string[]>;
    abstract chat(options?: AIChatOptions): Promise<AIChat>;

    async prompt(prompt: string, options?: AIChatOptions & AIRunOptions): Promise<AIResult> {
        const chat = await this.chat(options);
        try {
            return await chat.prompt(prompt, options);
        }
        finally {
            chat.close();
        }
    }

    async* stream(prompt: string, options?: AIChatOptions & AIRunOptions): AIStream {
        const chat = await this.chat(options);
        try {
            return yield* chat.stream(prompt, options);
        }
        finally {
            chat.close();
        }
    }
}

const _registered = new Map<string, AIEngine>();

async function engineFor(options?: AIChatOptions): Promise<AIEngine> {
    if (options?.engine != null) {
        const engine = await AI.getEngine(options.engine);
        if (!engine)
            throw new Error(`AI engine not found: '${options.engine}'`);
        return engine;
    }
    for (const engine of await AI.engines())
        if (await engine.available().catch(() => false))
            return engine;
    throw new Error('No AI engine is available. Install a package that provides one ' +
        '(e.g. Grokky), or register one with grok.ai.registerEngine().');
}

/**
 * Datagrok AI integration entry point.
 */
export const AI = {
    registerEngine(engine: AIEngine): void {
        _registered.set(engine.id, engine);
    },

    async engines(): Promise<AIEngine[]> {
        const funcs = Func.find({meta: {role: FUNC_TYPES.AI_ENGINE}});
        const discovered: AIEngine[] = await Promise.all(funcs.map((f) => f.apply()));
        return [..._registered.values(), ...discovered.filter((e) => !_registered.has(e.id))];
    },

    async getEngine(id: string): Promise<AIEngine | null> {
        return (await AI.engines()).find((e) => e.id === id) ?? null;
    },

    async chat(options?: AIChatOptions): Promise<AIChat> {
        return (await engineFor(options)).chat(options);
    },

    async prompt(prompt: string, options?: AIChatOptions & AIRunOptions): Promise<AIResult> {
        return (await engineFor(options)).prompt(prompt, options);
    },

    async* stream(prompt: string, options?: AIChatOptions & AIRunOptions): AIStream {
        return yield* (await engineFor(options)).stream(prompt, options);
    },

    get onPrompt(): rxjs.Observable<AITurnEvent> { return _events.onPrompt; },
    get onChunk(): rxjs.Observable<AITurnEvent & {chunk: string}> { return _events.onChunk; },
    get onResult(): rxjs.Observable<AITurnEvent & {result: AIResult}> { return _events.onResult; },
    get onCancel(): rxjs.Observable<AITurnEvent> { return _events.onCancel; },
    get onError(): rxjs.Observable<AITurnEvent & {error: any}> { return _events.onError; },

    /**
     * Configuration sub-namespace for AI provider settings.
     */
    config: {
        get current(): AIConfig {
            return api.grok_AI_Config();
        },

        get configured(): boolean {
            return this.current.configured;
        },

        get indexEntities(): boolean {
            const c = this.current;
            return c.configured && c.indexEntities;
        },

        get proxyUrl(): string {
            return api.grok_Dapi_OpenAI_Proxy();
        },

        get proxyToken(): string {
            return api.grok_Dapi_Get_Token();
        },
    },

    /**
     * Performs a semantic search over indexed Datagrok entities.
     *
     * The input [text] is converted into a vector embedding using the
     * configured embedding model and compared against stored entity
     * embeddings using vector distance.
     *
     * Results are filtered by the provided similarity [threshold],
     * optionally restricted to specific entity [types], and limited
     * by [limit].
     *
     * @param text Text query to search for. Must be a non-empty string.
     * @param threshold Similarity threshold in the range (0, 1).
     *   Higher values require closer (more similar) matches.
     * @param limit Maximum number of entities to return.
     * @param types Optional list of entity types to restrict the search. If omitted or empty, all
     *   entity types are considered. Use {@link Entity#entityType} to get type name.
     *
     * @returns A list of entities ordered by descending similarity
     *   (1.0 = perfect match).
     *
     * @throws Error if AI indexing is disabled or no AI provider is configured.
     */
    async searchEntities(
        text: string,
        threshold: number = 0.5,
        limit: number = 100,
        types?: string[],
    ): Promise<Entity[]> {
        const entities: Entity[] =
            await api.grok_AI_SearchEntities(text, threshold, limit, types ?? null);

        return entities.map((e: Entity) => toJs(e));
    },

    async processPrompt(prompt: string): Promise<boolean> {
        return !!(await api.grok_Shell_AI_Prompt(prompt));
    },
};
