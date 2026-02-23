# SonarPyVista

Utilities for building PyVista meshes from processed cavern sonar data.

---

## Methods

### `wireframe_from_cavlines2(lines)`

Builds a PyVista wireframe mesh from cavern survey line geometries.

**Parameters:**

- `lines` (`pandas.DataFrame` or `geopandas.GeoDataFrame`): DataFrame with 'geometry' column with 3D LineString objects.

**Returns:**

- `mesh` (`pyvista.PolyData`): Merged wireframe mesh.

---

### `add_rind(array_3d, rind_width)`

Add a rind (border) of zeros to a 3D array.

**Parameters:**

- `array_3d` (`np.ndarray`): 3D NumPy array to be padded.
- `rind_width` (`int`): Number of zeros to pad on each side of every axis.

**Returns:**

- `padded_array` (`np.ndarray`): Padded 3D array.

---

### `remove_rind(array_3d, rind_width)`

Remove rind (border) of zeros from a 3D array.

**Parameters:**

- `array_3d` (`np.ndarray`): 3D NumPy array to remove padding.
- `rind_width` (`int`): Number of layers to remove from each side of the array.

**Returns:**

- `array_cleaned` (`np.ndarray`): Array with padding removed.

---

### `mesh_xyz(xyz, iterations=3, n_iter=200, rf=0.1)`

Takes a geopandas.GeoDataFrame of processed xyz date (UTM) and creates a surface of the exterior of the cavern.  

**Parameters:**

- `xyz` (`geopandas.GeoDataFrame`): GeoDataFrame containing processed sonar data with columns: ['x', 'y', 'z', 'dx', 'dy', 'dz'].
- `iterations` (`int`, optional): Number of dilation/erosion passes, default is 3.
- `n_iter` (`int`, optional): Number of smoothing iterations, default is 200.
- `rf` (`float`, optional): Relaxation factor for smoothing (0 disables smoothing), default is 0.1.

**Returns:**

- `envelope` (`pyvista.PolyData`): The final triangulated mesh surface of the cavern.