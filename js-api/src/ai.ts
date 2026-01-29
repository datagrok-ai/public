import { IDartApi } from "./api/grok_api.g";
import { toJs } from "./wrappers";
import { Entity } from './entities';

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

/**
 * Datagrok AI integration entry point.
 */
export const AI = {
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
};
