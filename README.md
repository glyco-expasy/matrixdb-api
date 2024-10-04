# MatrixDB REST API

This repository contains the REST API for the MatrixDB project, providing access to interaction data involving components of the extracellular matrix (ECM), such as proteins, proteoglycans, and bioactive fragments. The API is written in Python and deployed using Docker with MongoDB and Solr as dependencies.

## Technologies Used

- **Python**: Backend language for building the REST API.
- **Flask**: For handling HTTP requests and routing.
- **MongoDB**: NoSQL database for storing ECM interaction data.
- **Solr**: For handling full-text search functionality.
- **Docker**: For containerized deployment and easy environment setup.
- **Docker Compose**: For orchestrating the multi-container setup (API, MongoDB, Solr).


## Installation

### Prerequisites

Ensure you have the following installed:

- [Docker](https://www.docker.com/get-started)
- [Docker Compose](https://docs.docker.com/compose/)
- [Python 3.8+](https://www.python.org/) (for local development)

### Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/MatrixDB-API.git
   cd MatrixDB-API

2. Set up environment variables:

Create a .env file in the root directory to configure environment-specific settings (MongoDB and Solr URLs, API secrets).

```
MONGO_URI=mongodb://mongo:27017/matrixdb
SOLR_URL=http://solr:8983/solr/matrixdb
```

3. Build and start the services using Docker Compose:

```
docker-compose up --build
```
This will start the API server, MongoDB, and Solr containers.

4. The API will be available at http://localhost:8000


## Usage

Once the services are running, you can interact with the API via HTTP requests. Here's an example of common API interactions:

```bash
# GET biomolecules by id 
curl http://localhost:8000/biomolecules/gag_1

# GET search query 'LOX' for basic search
curl api/search?query=LOX&mode=0 
