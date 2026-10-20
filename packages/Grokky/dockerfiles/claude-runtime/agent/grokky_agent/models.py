from typing import Literal, NotRequired, TypedDict

from .config import load_yaml


class Rates(TypedDict):
    input: float
    output: float
    cache_read: float
    cache_write: float


class LongContextRates(Rates):
    above_tokens: int


class Pricing(TypedDict):
    base: Rates
    long_context: NotRequired[LongContextRates]


class Model(TypedDict):
    family: Literal['claude', 'gpt']
    context: int
    images: bool
    prices: Pricing


class Usage(TypedDict, total=False):
    input_tokens: int | None
    output_tokens: int | None
    cache_read_tokens: int | None
    cache_write_tokens: int | None


CATALOG: dict[str, Model] = load_yaml('models.yaml')


def cost(model: str, usage: Usage) -> float:
    """USD for one call to the catalog model; a prompt over the long-context threshold is billed at those rates."""
    pricing = CATALOG[model]['prices']
    inp, out = usage.get('input_tokens') or 0, usage.get('output_tokens') or 0
    read, write = usage.get('cache_read_tokens') or 0, usage.get('cache_write_tokens') or 0
    long = pricing.get('long_context')
    rates = long if long and inp + read + write > long['above_tokens'] else pricing['base']
    return (inp * rates['input'] + out * rates['output'] + read * rates['cache_read'] +
            write * rates['cache_write']) / 1_000_000
