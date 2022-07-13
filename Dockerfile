FROM python:3.10
WORKDIR /app
COPY . /app/

RUN pip install poetry
RUN poetry config virtualenvs.in-project true
RUN poetry install --no-ansi
RUN poetry build
RUN pip install dist/pgscatalog_utils-0.1.0-py3-none-any.whl