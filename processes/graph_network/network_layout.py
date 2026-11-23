import numpy as np
from scipy.spatial import cKDTree

def separate_markers(positions, radii, gap):
    """Relax colliding circles in drawing units, with the spring layout as the starting point."""
    xy = np.array(positions, dtype=float, copy=True)
    radii = np.asarray(radii)
    if len(xy) < 2:
        return xy
    for _ in range(600):
        pairs = cKDTree(xy).query_pairs(2 * radii.max() + gap, output_type="ndarray")
        a, b = pairs.T
        delta = xy[b] - xy[a]
        lengths = np.linalg.norm(delta, axis=1)
        overlap = radii[a] + radii[b] + gap - lengths
        keep = overlap > 0.001
        if not keep.any():
            break
        a, b, delta, lengths = a[keep], b[keep], delta[keep], lengths[keep]
        coincident = lengths == 0
        if coincident.any():
            angle = a[coincident] * 2.4 + b[coincident]
            delta[coincident] = np.column_stack((np.cos(angle), np.sin(angle)))
            lengths[coincident] = 1
        shift = delta * (overlap[keep] / lengths / 2)[:, None]
        moves = np.zeros_like(xy)
        np.add.at(moves, a, -shift)
        np.add.at(moves, b, shift)
        neighbors = np.bincount(np.concatenate((a, b)), minlength=len(xy))
        xy += moves / np.maximum(1, neighbors)[:, None]
    return xy

def spaced_network_positions(pos, radius_by_node, drawing_span, gap):
    """Fit marker spacing in points, then return coordinates in the original units."""
    nodes = list(pos)
    if len(nodes) < 2:
        return dict(pos)
    coords = np.asarray([pos[node] for node in nodes], dtype=float)
    span = np.ptp(coords, axis=0).max()
    units_per_point = (span if span > 0 else 1) / drawing_span
    radii = [radius_by_node[node] for node in nodes]
    packed = separate_markers(coords / units_per_point, radii, gap)
    return dict(zip(nodes, packed * units_per_point))
