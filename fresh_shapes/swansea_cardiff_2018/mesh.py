import qmesh

def mesh(name):
    '''Todo: add docstring '''
    # Reading in the shapefile describing the domain boundaries, and creating a gmsh file.

    boundaries = qmesh.vector.Shapes()
    boundaries.fromFile('Severn_estuary_rev.shp')  # severn_cardiff.shp # Extended_BCM_lagoon_Swansea.shp

    loopShapes = qmesh.vector.identifyLoops(boundaries,
                                            isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
    polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=300,
                                                  meshedAreaPhysID=1)


    #
    inner = qmesh.vector.Shapes()
    inner.fromFile('inner.shp')
    inner_loops = qmesh.vector.identifyLoops(inner,
                                             fixOpenLoops=True, extraPointsPerVertex=10)

    inner_polygon = qmesh.vector.identifyPolygons(inner_loops, meshedAreaPhysID = 3)



    inner2 = qmesh.vector.Shapes()
    inner2.fromFile('inner_2.shp')
    inner2_loops = qmesh.vector.identifyLoops(inner2,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    inner2_polygon = qmesh.vector.identifyPolygons(inner2_loops, meshedAreaPhysID = 3)


    outer = qmesh.vector.Shapes()
    outer.fromFile('outer.shp')
    outer_loops = qmesh.vector.identifyLoops(outer,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    outer_polygon = qmesh.vector.identifyPolygons(outer_loops, meshedAreaPhysID = 4)


    outer2 = qmesh.vector.Shapes()
    outer2.fromFile('outer_2.shp')
    outer2_loops = qmesh.vector.identifyLoops(outer2,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    outer2_polygon = qmesh.vector.identifyPolygons(outer2_loops, meshedAreaPhysID = 4)


    inner_swan = qmesh.vector.Shapes()
    inner_swan.fromFile('inner_swan.shp')
    inner_swan_loops = qmesh.vector.identifyLoops(inner_swan,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    inner_swan_polygon = qmesh.vector.identifyPolygons(inner_swan_loops, meshedAreaPhysID = 5)


    outer_swan = qmesh.vector.Shapes()
    outer_swan.fromFile('outer_swan.shp')
    outer_swan_loops = qmesh.vector.identifyLoops(outer_swan,
                                             fixOpenLoops=True, extraPointsPerVertex=10)
    outer_swan_polygon = qmesh.vector.identifyPolygons(outer_swan_loops, meshedAreaPhysID = 6)

    ##################################################################################################
    # Create raster for mesh gradation towards full-resolution shorelines.

    ncresx, ncresy = 1300, 1300#700, 700#1000, 1000 #900, 900#800, 800#600, 600 # 1000, 1000 #1200, 1200 #600, 600
    # Create raster for mesh gradation towards full-resolution shorelines.


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



    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('cardiff_ref_structures.shp')
    grad_4 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_4.setShapes(GSHHS_coarser_boundaries)
    grad_4.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_4.setRasterResolution(ncresx, ncresy)
    grad_4.setGradationParameters(60., 5000.0, 1.0, 0.001)
    grad_4.calculateLinearGradation()
    grad_4.writeNetCDF('grad_card1.nc')



    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('cardiff_ref_impoundment.shp')
    grad_5 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_5.setShapes(GSHHS_coarser_boundaries)
    grad_5.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_5.setRasterResolution(ncresx, ncresy)
    grad_5.setGradationParameters(150., 5000.0, 1.0, 0.001)
    grad_5.calculateLinearGradation()
    grad_5.writeNetCDF('grad_card2.nc')


    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('grad_swansea.shp')
    grad_6 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_6.setShapes(GSHHS_fine_boundaries)
    grad_6.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_6.setRasterResolution(ncresx, ncresy)
    grad_6.setGradationParameters(40., 5000.0, 1.0, 0.001)
    grad_6.calculateLinearGradation()
    grad_6.writeNetCDF('grad_swan1.nc')


    GSHHS_coarser_boundaries = qmesh.vector.Shapes()
    GSHHS_coarser_boundaries.fromFile('lagoon_Swansea.shp')
    grad_7 = qmesh.raster.meshMetricTools.gradationToShapes()
    grad_7.setShapes(GSHHS_coarser_boundaries)
    grad_7.setRasterBounds(-8.0, -2.0, 50.0, 53.0)
    grad_7.setRasterResolution(ncresx, ncresy)
    grad_7.setGradationParameters(100., 5000.0, 1.0, 0.001)
    grad_7.calculateLinearGradation()
    grad_7.writeNetCDF('grad_swan2.nc')
    #
    # domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, inner_loops, inner_polygon)
    # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer_loops, outer_polygon)
    #
    # # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, inner_loops, inner_polygon)
    # # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer_loops, outer_polygon)
    #
    domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, inner2_loops, inner2_polygon) ###
    domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer2_loops, outer2_polygon) ###

    # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, inner2_loops, inner2_polygon)
    # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer2_loops, outer2_polygon)
    #
    # # domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, outer_swan_loops, outer_swan_polygon)
    # # domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, inner_swan_loops, inner_swan_polygon)
    #
    domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, inner_swan_loops,
                                                             inner_swan_polygon)
    domainLines, domainPolygons = qmesh.vector.insertRegions(domainLines, domainPolygons, outer_swan_loops,
                                                             outer_swan_polygon)

    # Calculate overall mesh-metric raster
    meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([grad_0, grad_1, grad_2, grad_3, grad_4, grad_5, grad_6, grad_7]) #grad_0
    meshMetricRaster.writeNetCDF('meshMetric.nc')
    # Create domain object and write gmsh files.
    domain = qmesh.mesh.Domain()
    domain.setGeometry(domainLines, domainPolygons)
    #domain.setGeometry(loopShapes, polygonShapes)
    domain.setMeshMetricField(meshMetricRaster)
    domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)
    # Meshing


    domain.gmsh(geoFilename= name + '.geo', \
                fldFilename= name + '.fld', \
                mshFilename= name + '.msh', \
                )



def convertMesh(name):
    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh( name + '.msh', 'EPSG:32630')
    mesh.writeShapefile(name)





if __name__ == '__main__':
    import qmesh
    import os

    # print (os.getcwd())
    # os.chdir("/data/cardiff_2018/inputs/")
    name = "swansea_cardiff_2018_7"
    # Initialising qgis API
    qmesh.initialise()

#    treatlines()
    mesh(name)
    convertMesh(name)


