# SonarDXF

 Utilities for constructing structured 3D line geometries from sonar survey point data.

---

## Methods

### `process_lines(xyz_gdf)`

Converts a GeoDataFrame of 3D points into horizontal and vertical LineString segments.

**Parameters:**

- ` xyz_gdf` (`geopandas.GeoDataFrame` or `pandas.DataFrame`): DataFrame of points with columns: 'x', 'y', 'z', 'depth', 'tilt', and 'cAziN'.

**Returns:**

- `lines` (`pandas.DataFrame`): DataFrame with a 'geometry' column containing 3D LineString objects.