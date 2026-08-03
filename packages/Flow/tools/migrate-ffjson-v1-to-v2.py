"""Convert .ffjson v1 (LiteGraph-backed) to v2 (Rete-native).

Usage:
    py migrate-ffjson-v1-to-v2.py <input.ffjson> [output.ffjson]

If output is omitted, overwrites the input in place. Run from anywhere.

What it does:
    - Reads the v1 LiteGraph payload (`graph.nodes` array + `graph.links` tuples).
    - For each node, emits a v2 entry with:
        * id stringified
        * typeName preserved (the v1 `type` field already matches the v2 registry)
        * label from `title` (user-renamed) or the last `/`-segment of `type`
        * pos pulled from the indexed-keys "0" / "1"
        * properties stripped of `_input_*` and `_passthroughCount`
        * inputValues recovered from `_input_<name>` keys
    - For each link tuple `[id, srcNode, srcSlot, dstNode, dstSlot, type]`, looks
      up the v1 nodes' input/output arrays to translate slot indices into v2
      string keys. Pass-through outputs become `<inputName>__pt`.
"""

import json
import sys
from pathlib import Path


def slot_key_for_input(node, slot_idx):
    """v2 input key = v1 input.name at that index."""
    inputs = node.get("inputs") or []
    if slot_idx < len(inputs):
        return inputs[slot_idx]["name"]
    return f"input_{slot_idx}"


def slot_key_for_output(node, slot_idx):
    """v2 output key.

    Pass-through outputs in v1 have name '→' at indices 0..passthroughCount-1
    and pair with the input at the same index. Real outputs use their name.
    """
    pt_count = (node.get("properties") or {}).get("_passthroughCount", 0)
    outputs = node.get("outputs") or []
    if slot_idx >= len(outputs):
        return f"output_{slot_idx}"
    out = outputs[slot_idx]
    if slot_idx < pt_count:
        # pass-through — derive key from corresponding input name
        inp_name = (node.get("inputs") or [{}])[slot_idx].get("name", f"in{slot_idx}")
        return f"{inp_name}__pt"
    return out["name"]


def convert_node(v1_node):
    nid = str(v1_node["id"])
    type_name = v1_node["type"]
    title = v1_node.get("title")
    label = title if title else type_name.rsplit("/", 1)[-1]
    pos = v1_node.get("pos", {})
    x = pos.get("0", 0) if isinstance(pos, dict) else (pos[0] if len(pos) > 0 else 0)
    y = pos.get("1", 0) if isinstance(pos, dict) else (pos[1] if len(pos) > 1 else 0)

    raw_props = v1_node.get("properties") or {}
    properties = {}
    input_values = {}
    for k, v in raw_props.items():
        if k.startswith("_input_"):
            input_values[k[len("_input_"):]] = v
        elif k == "_passthroughCount":
            continue
        else:
            properties[k] = v

    return {
        "id": nid,
        "typeName": type_name,
        "label": label,
        "pos": {"x": x, "y": y},
        "properties": properties,
        "inputValues": input_values,
    }


def convert_links(v1_graph, node_index):
    """Convert v1 link tuples to v2 connection objects.

    v1 link tuple: [linkId, srcNodeId, srcSlot, dstNodeId, dstSlot, type]
    """
    out = []
    for link in v1_graph.get("links") or []:
        if not link or len(link) < 5:
            continue
        link_id, src_id, src_slot, dst_id, dst_slot = link[0], link[1], link[2], link[3], link[4]
        src_node = node_index.get(src_id)
        dst_node = node_index.get(dst_id)
        if not src_node or not dst_node:
            continue
        out.append({
            "id": str(link_id),
            "source": str(src_id),
            "sourceOutput": slot_key_for_output(src_node, src_slot),
            "target": str(dst_id),
            "targetInput": slot_key_for_input(dst_node, dst_slot),
        })
    return out


def convert(doc_v1):
    if doc_v1.get("version") != "1.0":
        raise SystemExit(f"Expected v1.0, got {doc_v1.get('version')}")
    v1_graph = doc_v1.get("graph") or {}
    v1_nodes = v1_graph.get("nodes") or []
    node_index = {n["id"]: n for n in v1_nodes}

    nodes = [convert_node(n) for n in v1_nodes]
    connections = convert_links(v1_graph, node_index)

    metadata = doc_v1.get("metadata") or {}
    settings = (metadata.get("settings") or {})
    if not settings:
        settings = {
            "scriptName": doc_v1.get("name") or "Flow",
            "scriptDescription": doc_v1.get("description") or "",
            "tags": [],
        }

    return {
        "version": "2.0",
        "name": doc_v1.get("name") or settings.get("scriptName") or "Flow",
        "description": doc_v1.get("description") or settings.get("scriptDescription") or "",
        "author": doc_v1.get("author") or "unknown",
        "created": doc_v1.get("created"),
        "modified": doc_v1.get("modified"),
        "nodes": nodes,
        "connections": connections,
        "metadata": {"settings": settings},
    }


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    inp = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else inp
    with inp.open("r", encoding="utf-8") as f:
        v1 = json.load(f)
    v2 = convert(v1)
    with out.open("w", encoding="utf-8") as f:
        json.dump(v2, f, indent=2)
    print(f"converted: {inp} -> {out}  (nodes: {len(v2['nodes'])}, connections: {len(v2['connections'])})")


if __name__ == "__main__":
    main()
