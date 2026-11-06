"""Input validation shared by the API endpoints.

Request bodies used to go straight into MongoDB and Solr queries with no checks on type or size, so one request could carry
a million ids (or a dict where a list was expected). Everything an endpoint accepts is checked here and rejected with a 400.
"""
import json

from flask import request

MAX_IDS = 1000          # biomolecules / associations / experiments / expressions per request
MAX_NETWORK_IDS = 500   # the network call is the expensive one
MAX_ID_LENGTH = 64


class ApiError(Exception):
    def __init__(self, message, status=400):
        super().__init__(message)
        self.message = message
        self.status = status


def json_body():
    """The request body as parsed JSON, or a 400 when it is missing or malformed."""
    try:
        return json.loads(request.get_data() or b"null")
    except ValueError:
        raise ApiError("The request body is not valid JSON")


def string_list(value, name, max_items=MAX_IDS):
    """A list of short strings (ids), at most max_items of them."""
    if not isinstance(value, list):
        raise ApiError(f"'{name}' must be a list")
    if len(value) > max_items:
        raise ApiError(f"'{name}' has more than {max_items} entries", 413)
    for v in value:
        if not isinstance(v, str) or not v or len(v) > MAX_ID_LENGTH:
            raise ApiError(f"'{name}' must contain only ids (non-empty strings of at most {MAX_ID_LENGTH} characters)")
    return value


def ids_from_body(key, max_items=MAX_IDS):
    body = json_body()
    if key is None:                      # the body itself is the list
        return string_list(body, "body", max_items)
    if not isinstance(body, dict) or key not in body:
        raise ApiError(f"'{key}' is required")
    return string_list(body[key], key, max_items)


def found(document, what="resource"):
    """404 instead of an unhandled error when a lookup returned nothing."""
    if document is None:
        raise ApiError(f"Unknown {what}", 404)
    return document
