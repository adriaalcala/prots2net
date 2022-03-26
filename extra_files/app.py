"""App configuration module."""
from flasgger import Swagger  # type: ignore
from flask import Flask, request, jsonify
from flask_httpauth import HTTPBasicAuth  # type: ignore
from werkzeug.security import check_password_hash

from logger import logger
from service.routes.pubsub import api as pubsub_api
from service.constants import VALID_USERNAME_PASSWORD_PAIRS


app = Flask(__name__, static_url_path='/', static_folder='docs/_build/')
app.register_blueprint(pubsub_api)
auth = HTTPBasicAuth()


@auth.verify_password
def verify_password(username, password):
    if username in VALID_USERNAME_PASSWORD_PAIRS and \
            check_password_hash(VALID_USERNAME_PASSWORD_PAIRS.get(username), password):
        return username

app.config["SWAGGER"] = {
    "title": "Prots2Net",
    "uiversion": 3,
}

swag = Swagger(
    app,
    decorators=[auth.login_required],
    template={
        "swagger": "2.0",
        "info": {
            "title": "Prots2Net",
            "version": "1.0",
        },
        "consumes": [
            "application/json",
        ],
        "produces": [
            "application/json",
        ],
    },
)


@app.route('/docs', methods=['GET'])
@app.route('/docs/<path:path>', methods=['GET'])
@auth.login_required
def docs(path='index.html'):
    """Library documentation endpoint.
    ---
    tags:
        - Docs
    parameters:
        - name: path
          in: path
          description: Path to documentation
    responses:
        200:
            description: Documentation
    """
    return app.send_static_file(path)


if __name__ == '__main__':
    app.run()