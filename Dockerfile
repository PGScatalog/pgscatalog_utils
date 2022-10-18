FROM python:3.10 as builder
WORKDIR /app
COPY . /app/

RUN pip install poetry && poetry config virtualenvs.in-project true && \
    poetry install --no-ansi --no-dev
  
RUN poetry build

FROM python:3.10

WORKDIR /opt/

COPY --from=builder /app/dist/pgscatalog_utils-0.3.0-py3-none-any.whl .

RUN pip install pgscatalog_utils-0.3.0-py3-none-any.whl

RUN apt-get update && apt-get install -y sqlite3