FROM python:3.11-slim

WORKDIR /app

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc libffi-dev && \
    rm -rf /var/lib/apt/lists/*

# Python dependencies
COPY pyproject.toml requirements.txt ./
RUN pip install --no-cache-dir -e ".[web]" 2>/dev/null || \
    pip install --no-cache-dir \
    fastapi>=0.109.0 \
    "uvicorn[standard]>=0.27.0" \
    sqlalchemy>=2.0.25 \
    "python-jose[cryptography]>=3.3.0" \
    "passlib[bcrypt]>=1.7.4" \
    python-multipart>=0.0.6 \
    pydantic-settings>=2.1.0 \
    email-validator>=2.1.0

# Copy application
COPY . .
RUN pip install --no-cache-dir -e .

# Non-root user
RUN useradd --create-home appuser && \
    mkdir -p /app/job_results /app/data && \
    chown -R appuser:appuser /app
USER appuser

EXPOSE 8000

CMD ["uvicorn", "run_server:app", "--host", "0.0.0.0", "--port", "8000"]
