# Live-analysis API for the Cancer Protein Explorer.
# Build:  docker build -t cancer-explorer-api .
# Run:    docker run -p 8100:8100 cancer-explorer-api
FROM python:3.11-slim

# ProDy/SciPy need a C/C++ toolchain to build wheels that aren't prebuilt.
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install the package (with the API extra) first for better layer caching.
COPY pyproject.toml requirements.txt ./
COPY cancer_tool ./cancer_tool
RUN pip install --no-cache-dir ".[api]"

# Committed precomputed genes are served straight from disk.
COPY data ./data
COPY api ./api

EXPOSE 8100
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8100"]
