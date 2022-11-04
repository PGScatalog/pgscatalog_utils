
FROM python:3.10 as builder

# docker build --build-arg "ENV=PROD" ...

ARG ENV

RUN apt-get update && apt-get install -y sqlite3

WORKDIR /app

RUN pip install poetry

RUN python -m venv /venv

COPY install.sh poetry.lock pyproject.toml /app

RUN chmod +x install.sh && ./install.sh

COPY . . 

RUN poetry build && /venv/bin/pip install dist/*.whl

FROM builder as final

COPY --from=builder /venv /venv

ENV PATH="/venv/bin:${PATH}"


