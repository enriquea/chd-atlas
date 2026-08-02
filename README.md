# chd-atlas

A curated set of genes associated with CHD.

## Published data API

The atlas is served as static JSON built from the curated sources — no server,
no query language. Build it locally with:

```bash
uv run chd-atlas build --root . --out dist
```

The command refuses to write anything if `chd-atlas validate` would report an
error, so a published site always passed its own gate. Two builds of one commit
are byte-identical, and `manifest.json` carries a sha256 for every file served.

The shape of each file is documented in [docs/data-api.md](docs/data-api.md).
Consumers should read that before writing against the output — in particular the
note on contested genes, which must never be displayed as settled.

## Validating the corpus

```bash
uv run chd-atlas validate
```

Reports every problem it finds rather than stopping at the first, and exits 1 on
any error, 2 if `--root` is not a directory.
