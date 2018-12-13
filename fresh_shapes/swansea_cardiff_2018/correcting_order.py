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

    #lagoon extras
    polygonise('inner', 'inner_polygon', physID=9)
    polygonise('outer', 'outer_polygon', physID=8)
    polygonise('inner_2', 'inner2_polygon', physID=7)
    polygonise('outer_2', 'outer2_polygon', physID=6)
    polygonise('marked_lagoon','lagoon_marked_polygon',physID=5) #ok
    polygonise('marked_seaward','seaward_marked_polygon',physID=4) #ok

    #lagoons basins
    polygonise('lagoon','lagoon_polygon',physID=3) #ok
    polygonise('cardiff_lagoon_lines', 'cardiff_lagoon_polygon', physID=2)

    #channel basin
    polygonise('seaward','seaward_polygon',physID=1) #need seaward polygon

