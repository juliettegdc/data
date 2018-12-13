def polygonise(line_input, polygon_output, physID=2):
    boundaries = qmesh.vector.Shapes()
    boundaries.fromFile(line_input+'.shp')
    loopShapes = qmesh.vector.identifyLoops(boundaries,
                                            isGlobal=False, defaultPhysID=1000,
                                            fixOpenLoops=True)

    polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=400,
                                                  meshedAreaPhysID=physID)
    loopShapes.writeFile(line_input+str("_fixed"))
    polygonShapes.writeFile(polygon_output)

if __name__ == '__main__':
    import qmesh
    import os
    # Initialising qgis API
    qmesh.initialise()
    polygonise('marked_lagoon','lagoon_marked_polygon',physID=4)
    polygonise('marked_seaward','seaward_marked_polygon',physID=3)
    polygonise('lagoon','lagoon_polygon',physID=2)
    polygonise('seaward','seaward_polygon',physID=1)
