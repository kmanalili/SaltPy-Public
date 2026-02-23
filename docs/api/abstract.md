# Project Overview

## SaltPy: Open-Source Sonar Survey Processing for Solution-Mined Caverns

**Authors:**  
Matthew W. Bauer¹², Brandon Lampe¹, Kaitlyn Manalili¹², Baris Ates¹  

¹ Agapito Associates, LLC  
² Colorado School of Mines

---

## Citation

Bauer, M. W.; Lampe, B.; Manalili, K.; Ates, B. (2026):
SaltPy: Open-Source Sonar Survey Processing for Solution-Mined Caverns;
SMRI Spring 2026 Technical Conference, 26–29 April 2026, Edinburgh, Scotland.

## Abstract

Solution-mined caverns are designed and monitored using diverse datasets including surface surveys, well logs, and sonar surveys, often delivered in heterogeneous and vendor-specific formats. This variability complicates visualization, interpretation, and integration of cavern geometry across engineering and geospatial workflows.

SaltPy is an open-source, vendor-agnostic Python framework developed to processing, georeferencing, and visualizing sonar survey data from solution-mined caverns. The software addresses challenges associated with transforming proprietary sonar exports into spatially consistent, analysis-ready datasets suitable for engineering evaluation, GIS integration, and three-dimensional modeling.

SaltPy implements spherical-to-Cartesian coordinate transformations, enabling direct integration with common geospatial (GIS) environments. Cavern geometries are represented as point clouds using widely adopted geospatial data structures, promoting interoperability and long-term reproducibility. The framework also supports the generation of higher-order geometric products, including wireframe models, STL volumes, and volumetric cavern envelopes derived from voxel-based reconstruction with watertight repair utilities. These capabilities facilitate interactive three-dimensional visualization, interpretation, and export to common engineering formats.

SaltPy was developed to address the variability inherent in sonar survey data and common challenges in illustrating three-dimensional cavern geometries while minimizing technical debt. The software provides a reproducible and extensible foundation for salt cavern sonar interpretation and research. This discussion will provide an introduction to how SaltPy works and provide examples of how it may be implemented.

 