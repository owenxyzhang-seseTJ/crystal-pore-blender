import numpy as np
from scipy.ndimage import label


def periodic_components(mask):
    if np.asarray(mask).ndim != 3:
        raise ValueError("Periodic connectivity needs a three-dimensional mask.")
    local, n = label(mask)
    parent = np.arange(n + 1)
    weight = np.zeros((n + 1, 3), int)
    cycles = {i: [] for i in range(n + 1)}

    def find(i):
        path = []
        root = i
        while parent[root] != root:
            path.append(root)
            root = int(parent[root])
        accumulated = np.zeros(3, int)
        for node in reversed(path):
            accumulated = weight[node] + accumulated
            parent[node] = root
            weight[node] = accumulated
        return root, weight[i].copy()

    for axis in range(3):
        first = np.take(local, 0, axis)
        last = np.take(local, -1, axis)
        good = (first > 0) & (last > 0)
        pairs = np.unique(np.column_stack([last[good], first[good]]), axis=0)
        shift = np.eye(3, dtype=int)[axis]
        for a, b in pairs:
            ra, wa = find(int(a))
            rb, wb = find(int(b))
            if ra == rb:
                cycle = wb - wa - shift
                if cycle.any():
                    cycles[ra].append(cycle.tolist())
            else:
                parent[rb] = ra
                weight[rb] = shift + wa - wb
                cycles[ra].extend(cycles[rb])
    roots = np.array([find(i)[0] for i in range(n + 1)])
    sizes = np.bincount(local.ravel(), minlength=n + 1)
    components = []
    percolating = []
    for root in sorted(set(roots[1:])):
        ids = np.where(roots == root)[0]
        ids = ids[ids > 0]
        winding = np.unique(np.array(cycles[root], int).reshape(-1, 3), axis=0)
        axes = []
        if len(winding):
            for axis in range(3):
                target = np.eye(3)[axis]
                solution = np.linalg.lstsq(winding.T, target, rcond=None)[0]
                if np.linalg.norm(winding.T @ solution - target) < 1e-7:
                    axes.append(axis)
            percolating.extend(ids.tolist())
        components.append(
            {
                "root": int(root),
                "local_labels": ids.tolist(),
                "voxels": int(sizes[ids].sum()),
                "winding_vectors": winding.tolist(),
                "independent_pure_lattice_axes": axes,
                "percolates": bool(len(winding)),
                "dimension": int(np.linalg.matrix_rank(winding)) if len(winding) else 0,
            }
        )
    return local, np.isin(local, percolating), components
