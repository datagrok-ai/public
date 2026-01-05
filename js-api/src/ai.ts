import { IDartApi } from "./api/grok_api.g";
import { toJs } from "./wrappers";
import { Entity } from './entities';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export class AiPlugin {
    /**
     * Base URL of the Datagrok OpenAI proxy endpoint.
     *
     * All OpenAI-compatible API requests made through Datagrok should be sent
     * to this URL. The proxy forwards requests to the configured AI provider and
     * injects the correct authorization token.
     */
    get openAiProxyUrl(): string {
        return api.grok_Dapi_OpenAI_Proxy();
    }

    /**
     * Authentication token for requests sent to the Datagrok OpenAI proxy.
     *
     * This value is the **Datagrok access token of the currently logged-in user**.
     * It is used to authenticate the request with the Datagrok server.
     *
     * When a request is forwarded to the configured AI provider, Datagrok
     * **replaces this token with the server-configured OpenAI (or compatible)
     * API token** before sending the request to the external provider.
     *
     * Client code should always use this token when calling [openAiProxyUrl]
     * and should not supply a provider-specific API key directly.
     */
    get openAiProxyToken(): string {
        return api.grok_Dapi_Get_Token();
    }

    /**
     * Indicates whether an OpenAI-compatible AI provider is configured
     * on the Datagrok server.
     *
     * Returns `true` if the administrator has configured an AI provider
     * (for example, OpenAI or Azure OpenAI) and the platform is ready to
     * accept embedding or completion requests.
     */
    get openAiConfigured(): boolean {
        return api.grok_AI_OpenAiConfigured();
    }

    /**
     * Indicates whether AI-based entity indexing is enabled.
     *
     * When enabled, Datagrok maintains vector embeddings for supported
     * entities (projects, packages, connections, etc.), allowing semantic
     * search via [searchEntities].
     */
    get entityIndexingEnabled(): boolean {
        return api.grok_AI_EntityIndexingEnabled();
    }

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
    }
}
