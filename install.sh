#!/usr/bin/env sh

if [ ${ENV} = "DEV" ]; then
    poetry export --dev --without-hashes -f requirements.txt | /venv/bin/pip install -r /dev/stdin
else
    poetry export --without-hashes -f requirements.txt | /venv/bin/pip install -r /dev/stdin
fi
