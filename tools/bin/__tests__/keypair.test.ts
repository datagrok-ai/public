import {describe, it, expect, afterEach, vi} from 'vitest';
import crypto from 'crypto';
import * as kp from '../utils/keypair';

const URL_ = 'http://stand/api';

/** The vector `core/server/datlas/test/users/user_key_auth_test.dart` verifies. */
const PUBLIC_KEY: kp.Jwk = {
  kty: 'EC', crv: 'P-256',
  x: 'RRrfdF3svy_vsD6lcV3SdgpFsRd1jCkIoZI_GEf9LCg',
  y: 'WtzxWanjHSqG3XY9J5bEDXALqtM4wmogOLOAGtNhoRo',
};
const FINGERPRINT = '5rjE44EisLl_gyjodvZjZ77v-fAOCA0mHU6ct80IEw0';

const json = (body: any) => new Response(JSON.stringify(body), {status: 200, headers: {'content-type': 'application/json'}});

afterEach(() => { vi.unstubAllGlobals(); delete process.env.GROK_PRIVATE_KEY; });

describe('fingerprint', () => {
  it('is the RFC 7638 thumbprint the server derives from the same key', () => {
    expect(kp.fingerprint(PUBLIC_KEY)).toBe(FINGERPRINT);
  });

  it('ignores member order, since the thumbprint is over a canonical form', () => {
    const shuffled: kp.Jwk = {y: PUBLIC_KEY.y, kty: 'EC', x: PUBLIC_KEY.x, crv: PUBLIC_KEY.crv};
    expect(kp.fingerprint(shuffled)).toBe(FINGERPRINT);
  });

  it('refuses a key type neither side can verify', () => {
    expect(() => kp.fingerprint({kty: 'OKP', crv: 'Ed25519', x: 'abc'})).toThrow('Unsupported key type');
  });
});

describe('generateKeyPair / publicPart', () => {
  it('generates P-256 and keeps the private scalar out of the public half', () => {
    const {publicKey, privateKey} = kp.generateKeyPair();
    expect(privateKey.crv).toBe('P-256');
    expect(privateKey.d).toBeTruthy();
    expect(publicKey.d).toBeUndefined();
    expect(kp.publicPart(privateKey)).toEqual(publicKey);
    expect(kp.fingerprint(publicKey)).toBe(kp.fingerprint(kp.publicPart(privateKey)));
  });
});

describe('signNonce', () => {
  const verifies = (publicKey: kp.Jwk, text: string, signature: string): boolean =>
    crypto.verify('sha256', Buffer.from(text),
      {key: crypto.createPublicKey({key: publicKey as any, format: 'jwk'}), dsaEncoding: 'ieee-p1363'},
      Buffer.from(signature, 'base64url'));

  it('signs prefix + audience + nonce in the raw r||s form the server expects', () => {
    const {publicKey, privateKey} = kp.generateKeyPair();
    const signature = kp.signNonce(privateKey, URL_, 'a-nonce');
    expect(Buffer.from(signature, 'base64url')).toHaveLength(64);
    expect(verifies(publicKey, `datagrok-login:${URL_}:a-nonce`, signature)).toBe(true);
  });

  it('does not verify over the bare nonce, nor for another server', () => {
    const {publicKey, privateKey} = kp.generateKeyPair();
    const signature = kp.signNonce(privateKey, URL_, 'a-nonce');
    expect(verifies(publicKey, 'a-nonce', signature)).toBe(false);
    expect(verifies(publicKey, 'datagrok-login:http://elsewhere/api:a-nonce', signature)).toBe(false);
  });
});

describe('keyLogin', () => {
  it('asks for a nonce, signs it, and returns the token', async () => {
    const {privateKey} = kp.generateKeyPair();
    const calls: any[] = [];
    vi.stubGlobal('fetch', vi.fn(async (url: string, init: any) => {
      calls.push({url, body: JSON.parse(init.body)});
      return url.endsWith('/challenge') ? json({nonce: 'n-1'}) : json({isSuccess: true, token: 'Bearer t'});
    }));
    expect(await kp.keyLogin(URL_, privateKey)).toBe('Bearer t');
    expect(calls[0].url).toBe(`${URL_}/users/login/key/challenge`);
    expect(calls[0].body.fingerprint).toBe(kp.fingerprint(kp.publicPart(privateKey)));
    expect(calls[1].body.nonce).toBe('n-1');
    expect(calls[1].body.signature).toBeTruthy();
    expect(calls[1].body.audience).toBe(URL_);
  });

  it('reports the server’s reason rather than a bare failure', async () => {
    const {privateKey} = kp.generateKeyPair();
    vi.stubGlobal('fetch', vi.fn(async (url: string) =>
      url.endsWith('/challenge') ? json({nonce: 'n-1'}) : json({isSuccess: false, comment: 'unknown key'})));
    // `comment` is where UserLoginResponse carries the refusal, not `message`.
    await expect(kp.keyLogin(URL_, privateKey)).rejects.toThrow('unknown key');
  });

  it('names the version a server needs when it has no keypair routes', async () => {
    const {privateKey} = kp.generateKeyPair();
    // A server without them answers 404, or 401 because unknown paths are refused before routing.
    for (const status of [404, 401]) {
      vi.stubGlobal('fetch', vi.fn(async () => new Response('{"message":"Invalid session"}', {status})));
      await expect(kp.keyLogin(URL_, privateKey)).rejects.toThrow('needs Datagrok 1.28 or later');
      await expect(kp.keyLogin(URL_, privateKey)).rejects.toHaveProperty('name', 'ServerTooOldError');
    }
  });

  it('does not present an HTML error page as a token', async () => {
    const {privateKey} = kp.generateKeyPair();
    vi.stubGlobal('fetch', vi.fn(async () => new Response('<html>502</html>', {status: 502})));
    await expect(kp.keyLogin(URL_, privateKey)).rejects.toThrow('Unexpected response');
  });
});

describe('keypairFor', () => {
  it('takes the key from GROK_PRIVATE_KEY, raw or base64', () => {
    const {privateKey} = kp.generateKeyPair();
    process.env.GROK_PRIVATE_KEY = JSON.stringify(privateKey);
    expect(kp.keypairFor(URL_)).toEqual(privateKey);
    process.env.GROK_PRIVATE_KEY = Buffer.from(JSON.stringify(privateKey)).toString('base64');
    expect(kp.keypairFor(URL_)).toEqual(privateKey);
  });

  it('stands aside for an explicitly supplied developer key — that is another identity', () => {
    const {privateKey} = kp.generateKeyPair();
    process.env.GROK_PRIVATE_KEY = JSON.stringify(privateKey);
    expect(kp.keypairFor(URL_, 'someone-elses-key')).toBeUndefined();
  });

  it('is undefined when nothing is configured, so the caller falls back', () => {
    expect(kp.keypairFor('http://never-configured/api')).toBeUndefined();
  });
});

describe('keyFromEnv', () => {
  it('is undefined for an unset variable', () => {
    expect(kp.keyFromEnv('GROK_PRIVATE_KEY_NOT_SET')).toBeUndefined();
  });
});
