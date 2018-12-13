

def mesh(name):
    '''Todo: add docstring '''
    # Reading in the shapefile describing the domain boundaries, and creating a gmsh file.

    domainlines = qmesh.vector.Shapes()
    domainlines.fromFile('domainlines.shp')

    domainpolygons = qmesh.vector.Shapes()
    domainpolygons.fromFile('domainpolygons.shp')

    # boundaries = qmesh.vector.Shapes()
    # boundaries.fromFile('lagoon_Swansea.shp')
    #
    # loopShapes = qmesh.vector.identifyLoops(boundaries,
    #                                         isGlobal=False, defaultPhysID=1000,
    #                                         fixOpenLoops=True)
    # polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=300,
    #                                               meshedAreaPhysID=1)

    ncresx, ncresy = 600,600
    # Create raster for mesh gradation towards full-resolution shorelines.
    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_swansea.shp')
    grad_5 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_5.setShapes(GSHHS_fine_boundaries)
    grad_5.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_5.setRasterResolution(ncresx, ncresy)
    grad_5.setGradationParameters(20.0, 5000.0, 1.0, 0.002)
    grad_5.calculateLinearGradation()
    grad_5.writeNetCDF('grad_isl.nc')

    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_isl.shp')
    grad_0 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_0.setShapes(GSHHS_fine_boundaries)
    grad_0.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_0.setRasterResolution(ncresx, ncresy)
    grad_0.setGradationParameters(100.0, 5000.0, 0.5, 0.001)
    grad_0.calculateLinearGradation()
    grad_0.writeNetCDF('grad_isl.nc')

    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_100.shp')
    grad_1 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_1.setShapes(GSHHS_fine_boundaries)
    grad_1.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_1.setRasterResolution(ncresx, ncresy)
    grad_1.setGradationParameters(150.0, 10000.0, 1)
    grad_1.calculateLinearGradation()
    grad_1.writeNetCDF('grad_100.nc')

    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_200.shp')
    grad_2 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_2.setShapes(GSHHS_fine_boundaries)
    grad_2.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_2.setRasterResolution(ncresx, ncresy)
    grad_2.setGradationParameters(300.0, 10000.0, 1)
    grad_2.calculateLinearGradation()
    grad_2.writeNetCDF('grad_200.nc')

    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('grad_500.shp')
    grad_3 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_3.setShapes(GSHHS_coarser_boundaries)
    grad_3.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_3.setRasterResolution(ncresx, ncresy)
    grad_3.setGradationParameters(1000.0, 10000.0, 1.0)
    grad_3.calculateLinearGradation()
    grad_3.writeNetCDF('grad_500.nc')



    # Calculate overall mesh-metric raster
    meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([grad_1, grad_2, grad_3, grad_0, grad_5]) #grad_0
    meshMetricRaster.writeNetCDF('meshMetric.nc')
    # Create domain object and write gmsh files.
    domain = qmesh.mesh.Domain()
    domain.setGeometry(domainlines, domainpolygons)
    # domain.setGeometry(loopShapes, polygonShapes)
    domain.setMeshMetricField(meshMetricRaster)
    domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)
    # Meshing

    domain.gmsh(geoFilename= name + '.geo', \
                fldFilename= name + '.fld', \
                mshFilename= name + '.msh', \
                )

#
def convertMesh(name):
    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh( name + '.msh', 'EPSG:32630')
    mesh.writeShapefile(name)


if __name__ == '__main__':
    import qmesh
    import os
    print os.getcwd()
    os.chdir("/data/fresh_shapes/swansea_lagoon/")
    name = "swansea_2018"
    # Initialising qgis API
    qmesh.initialise()
#    treatlines()
    mesh(name)
    convertMesh(name)