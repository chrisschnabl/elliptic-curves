# P079 Lab 01

Python implementation to use Elliptic Curves to get Diffie-Hellman with Curve25519 (X25519) and Signatures with Twisted Edwards Curve25519 (Ed25519). This is not production ready code that should be used in any real application.

# Run

## Locally

Requires `uv` to be installed.

```bash
uv run src/example.py
```

## Dockerized

```bash
./run.sh
```

# Run tests

## Locally

Requires `uv` to be installed.

```bash
uv run pytest
```

## Dockerized

```bash
./run_tests.sh
```
