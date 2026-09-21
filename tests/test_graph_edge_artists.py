import ast
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

# Load the rendering helper without triggering the CLI's global style setup.
source = Path(__file__).resolve().parents[1] / 'processes/graph_network/graph_network.py'
node = next(n for n in ast.parse(source.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'set_edge_artist_zorder')
namespace = {}
exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), 'exec'), namespace)
set_zorder = namespace['set_edge_artist_zorder']


class EdgeArtistTests(unittest.TestCase):
    def test_networkx_return_types_render(self):
        graphs = [nx.empty_graph(1), nx.path_graph(2), nx.DiGraph([(0, 1)]), nx.MultiGraph([(0, 1), (0, 1)])]
        for graph in graphs:
            with self.subTest(graph=type(graph).__name__, edges=graph.number_of_edges()):
                fig, ax = plt.subplots()
                try:
                    artists = nx.draw_networkx_edges(graph, nx.circular_layout(graph), ax=ax)
                    set_zorder(artists, 1)
                    for artist in artists if isinstance(artists, list) else [artists]:
                        if artist is not None:
                            self.assertEqual(artist.get_zorder(), 1)
                    with tempfile.TemporaryDirectory() as directory:
                        fig.savefig(Path(directory)/'network.pdf')
                finally:
                    plt.close(fig)
        set_zorder(None, 1)
