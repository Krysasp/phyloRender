# phyloRender – an R package for phylodynamic /phylogeographic analyses
## Overview
phyloRender contains scripts to analyzing, and visualizing phylogenetic and mutation data. The package is a component of the on-going development program (phyloAct) that enables end-to-end phylogenetic workflows, from raw data handling to advanced visualization. phyloRender performs all downstream analyses for outputs from Part I of phyloAct. The R functions were designed to manipulate, analyze, and visualize BEAST MCC phylogenetic trees with spatial and trait annotations. It supports:
 * Extracting tip coordinates and node metadata.
 * Adding discrete and continuous traits.
 * Managing Highest Posterior Density (HPD) intervals in latitude/longitude.
 * Flexible visualization, including swapping lat/lon, bounding ellipses, and animations.
 * Tree transformations (cladogram, radial, tip/branch highlighting).

All functions are **modular**, allowing them to be integrated into your workflow in main script.
________________________________________
## **Table of Contents**
* **Tree Loading and Parsing**
* **Trait Management**
* **HPD Handling**
* **Coordinate Transformations**
* **Visualization / Animation Support**
* **Examples**
* **Support**
________________________________________
### 1. Tree Loading and Parsing
 #### read_mcc_tr()
Purpose:
Load a BEAST MCC tree, parse node and tip metadata, and compute node heights and posterior probabilities.
Parameters:
| **Parameter**	|	   **Type**			   |  **Default**	|   	      **Description**		                    |
|---------------|------------------|--------------|---------------------------------------------- |
|     trName    |      string      |         -    | Path to the tree file (.tre, .nexus)          |
|      sep      |      string      |       `"\    |                   "`                          |
|     BEAST2    |      logical     |     FALSE    | Whether the tree comes from BEAST2 (adjusts parsing).|
|               |                  |              |                                               |

________________________________________
### 2. Trait Management
```bash
addDiscreteTraits(tr)
```
Adds discrete traits from BEAST MCC annotations.
Features:
  * Extracts all discrete traits listed in the tree annotations.
  * Stores trait matrices in tr$props and unique values in tr$uprops.
```bash
getFullTraitSet(tr, propIndex=1)
```
Generates a full probability matrix for a given discrete trait.
Parameters:
  * propIndex – choose which trait column to extract.
```bash
addFullTraitSet(tr, propIndex=1)
```
Adds the full probability matrix to the tree.
```bash
addContinuousTraits(tr, traitNames=NULL)
```
Extracts continuous traits from the MCC tree annotations (e.g., rate, location).
    * traitNames=NULL will automatically detect all continuous traits.
    * Stores results as matrices in tr$continuous_traits.
________________________________________
### 3. HPD Handling
```bash
addLatLonHPD(tr, latlonName="latlon", hpdNum=1)
```
Adds HPD intervals for tip latitude and longitude.
Parameters:
  * latlonName – prefix of lat/lon in metadata (default "latlon").
  * hpdNum – HPD level to extract (e.g., 1 = 80% HPD).
Future Enhancements:
 * Multi-HPD extraction (support multiple HPD levels).
 * Compute HPD centroids for visualization.
 * Fit bounding ellipses to HPDs for plotting.
```bash
fit_HPDs_to_standard(tr, npts=50, ltol=0.005)
```
Fits HPD points to a standard ellipse shape for visualization.
Parameters:
  * npts – number of points to fit per HPD.
  * ltol – tolerance threshold for degenerate points.
```bash
interpolate_discrete_element(fromSet, toSet, fractTime)
```
Interpolate discrete traits along branches for animation.
________________________________________
### 4. Coordinate Transformations
```bash
switch_latlon_advanced(tr, swapTips=TRUE, swapHPDs=TRUE, hpdLevels=NULL, plotCheck=FALSE)
```
Purpose: Swap latitude and longitude in MCC trees.
Features:
  * Swaps tip coordinates (tr$latlon).
  * Swaps HPD node coordinates (tr$lat80, tr$lon80).
  * Works with multiple HPD levels.
  * Optional plot for sanity check.
Parameters:
|  **Parameter**  |  **Type**   |  Default**  |            **Description**             |
|-----------------|-------------|-------------|----------------------------------------|
|    swapTips     |   logical   |    TRUE     |  Swap tip coordinates?                 |
|    swapHPDs     |   logical   |    TRUE     |    Swap HPDs?                          |
|    hpdLevels    |   integer   |    NULL     | Which HPD levels to swap; NULL = all.  |
|    plotCheck    |   logical   |    FALSE    |  Plot tips before/after swap.          |
________________________________________
### 5. Visualization / Animation Support
```bash
animate_phyloRender(tr, times, outdir="frames", prefix="frame", xlim=c(-180,180), ylim=c(-60,80), propIndex=0)
```
Generates a series of PNG frames of spatial phylodynamics over time.
Highlights tip points and HPD polygons.
* Supports discrete traits via propIndex.
* Can be combined with switch_latlon_advanced() for coordinate flipping.
________________________________________
6. Examples of usage for this source script as follows:
```bash
# Load tree
tr <- read_mcc_tr("example_tree.nexus")

# Add continuous traits
tr <- addContinuousTraits(tr)

# Add discrete traits
tr <- addDiscreteTraits(tr)
tr <- addFullTraitSet(tr, propIndex=1)

# Add HPD intervals
tr <- addLatLonHPD(tr, hpdNum=1)
tr <- fit_HPDs_to_standard(tr)

# Swap latitude and longitude (optional)
tr <- switch_latlon_advanced(tr, swapTips=TRUE, swapHPDs=TRUE, plotCheck=TRUE)

# Animate spatial phylodynamics
times <- seq(0, max(tr$node.times), length.out=100)
animate_phylomovie(tr, times, outdir="frames", prefix="frame", propIndex=1)

# Export discrete trait matrix to CSV
write.csv(tr$fullset, "trait_matrix.csv")
```
________________________________________
7. Support
  •	Requirements:
    o	R >= 4.4.2
    o	Packages: ape, conicfit, graphics
  •	Tips:
      o	Use plotCheck=TRUE to verify transformations.
      o	HPD ellipses allow visualization of uncertainty.
      o	Functions are compatible with main.R pipeline.
  •	Contact / Help:
      o	Questions about usage, bug reports, or feature requests can be directed to jonatasp92@gmail.com
________________________________________
