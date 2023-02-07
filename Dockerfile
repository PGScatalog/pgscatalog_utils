
FROM python:3.10 as builder

# docker build --build-arg "ENV=PROD" ...

ARG ENV

WORKDIR /app

RUN pip install poetry

RUN python -m venv /venv

COPY install.sh poetry.lock pyproject.toml /app/

RUN chmod +x install.sh && ./install.sh

COPY . . 

RUN poetry build && /venv/bin/pip install dist/*.whl

FROM python:3.10.9-slim-bullseye

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

COPY --from=builder /venv /venv

ENV PATH="/venv/bin:${PATH}"


