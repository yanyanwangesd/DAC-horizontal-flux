# DAC-horizontal-flux
This landscape evolution model Divide and Capture (DAC) version calculates basin projection area, upstream-averaged erosion rate, and horizontal mass flux. 

In the DAC framework, grid nodes represent rivers nodes with their topological connections dynamically stored and updated. For each grid node, the model tracks all upstream nodes that ultimately drain to it, recording their locations, elevations, local erosion rates, and local contributing areas.

## NEW in this DAC branch 
1. __basin projection area__
The surface of a drainage basin is three-dimensional. Conventional basin area in geomorphology/hydrology studies is the projection area of the 3D basin surface onto the horizontal plane. However, sometimes a projection area to a vertical plane is needed for researches of horizontal fluxes coming out of a basin. This code package implemented the method published in *__Wang and Willett 2021__* (https://doi.org/10.5194/esurf-9-1301-2021). 

2. __upstream averaged erosion rate__
The upstream-averaged erosion rate is calculated by integrating the eroded mass material and dividing it by the basin area. 

3. __horizontal mass flux__
This code package implemented the method published in *__Wang and Willett 2021__* (https://doi.org/10.5194/esurf-9-1301-2021). 

## Derivatives of this model
This model is initially published in accompany with a research paper that is recently submitted to GRL for peer-review. If you want to develop further upon this model, please contact the author: Yanyan Wang through email: wangyanyan0607@hotmail.com, or yanyan.wang@eaps.etha.ch -----July 1st, 2025

