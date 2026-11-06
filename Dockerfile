FROM python:3.12
# do not run the server as root
RUN useradd --create-home --uid 10001 app
WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=app:app . .

ENV PYTHONPATH="${PYTHONPATH}:./src"
USER app
CMD ["uwsgi", "--ini", "uwsgi.ini", "--http", ":8000", "--wsgi-file", "src/matrixdb/api/main.py", "--callable", "app", "--processes", "2", "--threads", "2"]
