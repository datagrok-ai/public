const TYPE_UINT_8_LIST: u16 = 0xA103;

pub fn frame_with_header(header_json: &str, payload: &[u8]) -> Vec<u8> {
    let json = header_json.as_bytes();
    let mut out = Vec::with_capacity(10 + json.len() + payload.len());
    out.extend_from_slice(&TYPE_UINT_8_LIST.to_le_bytes());
    out.extend_from_slice(&(json.len() as f64).to_le_bytes());
    out.extend_from_slice(json);
    out.extend_from_slice(payload);
    out
}
