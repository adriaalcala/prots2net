import dash
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import gi

gi.require_version("Gtk", "3.0")
gi.require_version('WebKit2', '4.0')
from gi.repository import Gtk, WebKit2

from flask import request

# ======== #
# DASH APP #
# ======== #
# Dash uses the Flask web framework
# For the Dash deployment, we need to access the Flask application instance
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"

app = dash.Dash(
    __name__,
    routes_pathname_prefix="/dash/",
    external_stylesheets=[dbc.themes.BOOTSTRAP, FONT_AWESOME]
)


### Load files
@app.callback(
    Output('B_add-files', 'value'),
    [Input('B_add-files', 'n_clicks')],
    [State('files-check', 'value')] )
def select_files(n_clicks, options):
    window = Gtk.Window()
    dialog =  Gtk.FileChooserDialog(title="Please choose a file", parent=window, action=Gtk.FileChooserAction.OPEN)
    dialog.add_buttons(
            Gtk.STOCK_CANCEL,
            Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN,
            Gtk.ResponseType.OK,
        )

    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        print("Open clicked")
        print("File selected: " + dialog.get_filename())
    elif response == Gtk.ResponseType.CANCEL:
        print("Cancel clicked")

    dialog.destroy()
    window.destroy()
    return str(n_clicks)

# ======= #
# RUN APP #
# ======= #
server = app.server
app.config['suppress_callback_exceptions'] = True

# Run Dash on a Public IP
if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=int("8080"), debug=True)