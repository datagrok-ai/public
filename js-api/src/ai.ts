import { IDartApi } from "./api/grok_api.g";
import { toJs } from "./wrappers";
import { Entity } from './entities';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export type AiProvider = "openai" | "azure";
export type AzureApiMode = "openai_compat" | "legacy";
export type AiConfig =
    {
    provider: AiProvider;
    configured: boolean;
    indexEntities: boolean;
    } & {
        apiMode?: AzureApiMode;
        apiVersion?: string;
        modelToDeployment?: Record<string, string>;
    };


/**
 * Datagrok AI integration entry point.
 *
 * Client usage depends on the configured AI provider:
 *
 * ## OpenAI (public API)
 * Use the standard OpenAI JS client:
 *
 * ```ts
 * import OpenAI from "openai";
 *
 * const client = new OpenAI({
 *   baseURL: `${grok.ai.openAiProxyUrl}/v1`,
 *   dangerouslyAllowBrowser: true,
 *   apiKey: grok.ai.openAiProxyToken
 * });
 * ```
 *
 * ## Azure OpenAI
 * Datagrok forwards requests to Azure OpenAI.
 *
 * - `apiMode === "openai_compat"` (recommended):
 *   - Use the same OpenAI client as above.
 *   - Pass **Azure deployment names** as `model`
 *     (use `modelToDeployment` mapping if provided).
 *
 * - `apiMode === "legacy"`:
 *   - Use an Azure-aware OpenAI client (deployment + api-version semantics).
 *   - Pass **Azure deployment names** as `model`
 *     (use `modelToDeployment` mapping if provided).
 *
 * In all cases, client code authenticates only with Datagrok.
 * Provider credentials are injected server-side.
 */
export class AiPlugin {
    /**
     * Current AI provider configuration.
     *
     * How to choose the client:
     *
     * - `provider: "openai"`
     *   - Use the standard OpenAI JS client.
     *   - Use normal OpenAI model IDs (e.g. `gpt-4o`, `text-embedding-3-small`).
     *
     * - `provider: "azure"`
     *   - Requests are forwarded to Azure OpenAI through Datagrok.
     *   - `apiMode: "openai_compat"`:
     *     - Use the standard OpenAI JS client.
     *     - Azure usually expects deployment names as `model`.
     *     - Use `modelToDeployment` to map OpenAI model IDs to deployment names.
     *   - `apiMode: "legacy"`:
     *     - Prefer an Azure-aware client.
     */
    get config(): AiConfig {
        return api.grok_AI_Config();
    }

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
