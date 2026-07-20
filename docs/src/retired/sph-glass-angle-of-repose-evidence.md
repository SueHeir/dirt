# Retired SPH repose: external-evidence audit

## Decision

**No glass angle-of-repose number, rolling-friction coefficient, or validation
claim is admitted for current DIRT.**  This is a negative evidence decision,
not a replacement calibration.  It supplements the
[retirement record](./sph-glass-angle-of-repose.md) by recording an external
adversarial check of the citation that had been associated with the historical
22--26 degree band.

## Independent source check

On 2026-07-19, a query was sent directly to the Crossref Works API, outside
the DIRT source tree:

```text
https://api.crossref.org/works/10.1016/S0378-4371(99)00183-1
```

The returned bibliographic record identifies DOI
`10.1016/s0378-4371(99)00183-1` as Zhou, Wright, Yang, Xu, and Yu (1999),
*Rolling friction in the dynamic simulation of sandpile formation*, *Physica
A: Statistical Mechanics and its Applications*, a journal article.  A compact
copy of the metadata response, normalized with `jq -c`, is retained below;
its SHA-256 is
`3c39b2e2124f52dc9d3b38b5b0359404e1f6a6a386fb40ba1748e00429838ba8`.

```json
{"api_status":"ok","doi":"10.1016/s0378-4371(99)00183-1","title":"Rolling friction in the dynamic simulation of sandpile formation","container_title":"Physica A: Statistical Mechanics and its Applications","published":[1999,7],"type":"journal-article","authors":["Y.C. Zhou","B.D. Wright","R.Y. Yang","B.H. Xu","A.B. Yu"]}
```

This is deliberately an adversarial result: the independently supplied title
describes a **dynamic simulation of sandpile formation**, rather than a
protocol-matched experiment on the glass material and apparatus claimed by a
future calibration.  It verifies bibliographic identity and exposes a scope
mismatch; it does not establish a numerical angle, a material specification,
or a contact-law parameter.

For a second, non-identical check, Crossref's title-search endpoint returned
the target record together with a later comment and a different 2011
simulation-validation article.  The normalized three-result list has SHA-256
`8ce8999b561491b04d5fca581ad190426a9db12fe221c9831cff8ca0a1fe9023`:

```json
[{"doi":"10.1016/s0378-4371(99)00183-1","title":"Rolling friction in the dynamic simulation of sandpile formation","type":"journal-article"},{"doi":"10.1016/j.physa.2005.01.019","title":"Comment on \"Rolling friction in the dynamic simulation of sandpile formation\"","type":"journal-article"},{"doi":"10.1061/(asce)gm.1943-5622.0000156","title":"Identification and Validation of Rolling Friction Models by Dynamic Simulation of Sandpile Formation","type":"journal-article"}]
```

The search confirmation reduces the risk of a mistyped DOI or title, but it is
still bibliographic metadata, not an experimental observation.  Neither check
is coupled to a DIRT executable, a DIRT acceptance threshold, or a fitted
parameter, so neither can tautologically make this repository pass.

A third, separately operated bibliographic index, OpenAlex, resolves the DOI
to the same normalized DOI, title, 1999 year, and article type. This protects
against a transcription error or a single-service anomaly; it still says
nothing about glass material, apparatus, formation procedure, measured angle,
or a rolling-contact law. Three agreeing metadata services are therefore
stronger evidence of **identity only**, not evidence of comparability.

## Reproduction and limits

Anyone can independently repeat the three external queries and compare their
normalized outputs (metadata can legitimately change over time):

```bash
curl -LfsS 'https://api.crossref.org/works/10.1016/S0378-4371(99)00183-1' |
  jq -c '{api_status:.status,doi:.message.DOI,title:.message.title[0],container_title:.message["container-title"][0],published:.message.published["date-parts"][0],type:.message.type,authors:[.message.author[] | "\(.given) \(.family)"]}' |
  sha256sum
curl -LfsS 'https://api.crossref.org/works?query.title=Rolling%20friction%20in%20the%20dynamic%20simulation%20of%20sandpile%20formation&rows=3' |
  jq -c '[.message.items[] | {doi:.DOI,title:.title[0],type:.type}]' |
  sha256sum
curl -LfsS 'https://api.openalex.org/works/https://doi.org/10.1016/S0378-4371(99)00183-1' |
  jq -c '{doi,title,publication_year,type}'
```

This audit did **not** obtain or read the full paper, inspect a primary glass
experiment, recover its apparatus/protocol, run a solver, or assess
uncertainty.  It therefore cannot decide whether any particular future
22--26 degree target is correct.  It only demonstrates why this citation
cannot carry that burden by itself.  Admission of a future claim still
requires the independent material/protocol evidence, solver ledger,
replicates, and cross-code reconciliation listed in the retirement record;
those requirements have not been relaxed.

This document and the normalization commands were prepared with AI assistance.
The external responses were retrieved on the stated date, but no human
literature review or experimental validation has been performed.  Treat this
as a transparent provenance boundary, not as scientific authorship or a
substitute for independent expert review.
