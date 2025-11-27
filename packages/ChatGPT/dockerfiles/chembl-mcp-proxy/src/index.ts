import express from "express";
import {McpServer} from "@modelcontextprotocol/sdk/server/mcp.js";
import {StreamableHTTPServerTransport} from "@modelcontextprotocol/sdk/server/streamableHttp.js";
import {MCPServerStdio} from "@openai/agents";
import path from "path";
import {fileURLToPath} from "url";
import {z} from "zod";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const REMOTE_MCP_PATH = path.resolve(__dirname, '../chembl-mcp-server/build/index.js');

const underlyingMcp = new MCPServerStdio({
  command: "node",
  args: [REMOTE_MCP_PATH],
});

const server = new McpServer({
  name: "proxy-server",
  version: "1.0.0",
});

function convertToZod(jsonSchema: any) {
  if (!jsonSchema || Object.keys(jsonSchema).length === 0) {
    return z.object({});
  }

  const shape: Record<string, any> = {};
  for (const key in jsonSchema.properties) {
    const prop = jsonSchema.properties[key];
    if (prop.type === "string") shape[key] = z.string();
    else if (prop.type === "number") shape[key] = z.number();
    else if (prop.type === "boolean") shape[key] = z.boolean();
    else shape[key] = z.any();
  }

  return z.object(shape);
}

async function registerUnderlyingTools() {
  if (!underlyingMcp.listTools) {
    console.warn("Underlying MCP does not support listTools()");
    return;
  }

  const tools = await underlyingMcp.listTools();
  if (!tools) return;

  for (const tool of tools) {
    server.registerTool(
      tool.name,
      {
        title: tool.title || tool.name,
        description: tool.description,
        inputSchema: convertToZod(tool.inputSchema),
        outputSchema: convertToZod(tool.outputSchema),
      },
      async (input: any) => {
        const result = await underlyingMcp.callTool(tool.name, input);
        return {
          content: [{ type: "text", text: JSON.stringify(result) }],
          structuredContent: {result},
        };
      }
    );
  }
}

const app = express();
app.use(express.json());

app.post("/mcp", async (req, res) => {
  const transport = new StreamableHTTPServerTransport({
    sessionIdGenerator: undefined,
    enableJsonResponse: true,
  });

  res.on("close", () => transport.close());

  await server.connect(transport);
  await transport.handleRequest(req, res, req.body);
});

(async () => {
  try {
    await underlyingMcp.connect();
    await registerUnderlyingTools();
    const PORT = parseInt(process.env.PORT || "3333", 10);
    app.listen(PORT, () =>
      console.log(`Proxy MCP Server running on http://localhost:${PORT}/mcp`)
    );
  } catch (err) {
    console.error("Failed to start underlying MCP server:", err);
    process.exit(1);
  }
})();

process.on("SIGINT", async () => {
  console.log("Closing underlying MCP server...");
  await underlyingMcp.close();
  process.exit();
});
