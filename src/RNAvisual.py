'''
RNA Secondary Structure Visualizations using dash
Able to choose which structure to visualize and select multiple for comparison
'''
import dash
import dash_bio as dashbio
import dash_core_components as dcc
import dash_html_components as html
from dash.exceptions import PreventUpdate
import pandas as pd
import sys


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


def set_sequences():
    print("seq")
    the_sequences = {}
    seq_data = pd.read_csv(sys.argv[1], sep="\t")
    for index, row in seq_data.iterrows():
        name = str(int(row.SNP))
        structure = {
            'sequence': row.SEQ,
            'structure': row.SEQ_DB
        }
        snp_name = name + "_SNP"
        snp_structure = {
            'sequence': row.SNP_SEQ,
            'structure': row.SNP_SEQ_DB
        }
        the_sequences.update({name: structure, snp_name: snp_structure})
    return the_sequences


sequences = set_sequences()

key = list(sequences.keys())[0]

app.layout = html.Div([
    dashbio.FornaContainer(
        id='forna'
    ),
    html.Hr(),
    html.P('Select the sequences to display below.'),
    dcc.Dropdown(
        id='forna-sequence-display',
        options=[
            {'label': name, 'value': name} for name in sequences.keys()
        ],
        multi=True,
        value=[key]
    )
])


@app.callback(
    dash.dependencies.Output('forna', 'sequences'),
    [dash.dependencies.Input('forna-sequence-display', 'value')]
)
def show_selected_sequences(value):
    if value is None:
        raise PreventUpdate
    return [
        sequences[selected_sequence]
        for selected_sequence in value
    ]


if __name__ == '__main__':
    global seq_file
    seq_file = sys.argv[1]
    print(seq_file)
    app.run_server(debug=True)
