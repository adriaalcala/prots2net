from waitress import serve
from app import server #so "app" is the name of my Dash script I want to serve

#Run Dash on a Public IP
if __name__ == "__main__":
    serve(server)