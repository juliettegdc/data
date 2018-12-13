from thetis import *
import pickle
import utm


def coriolis(mesh, lat,):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

def initialise_fields(mesh2d, inputdir, outputdir, identifier, manning = False):

    CG_2d = FunctionSpace(mesh2d, 'CG', 1)

    with timed_stage('reading bathymetry'):
        chk = DumbCheckpoint(inputdir+"/bathymetry2D", mode=FILE_READ)
        bathymetry_2d = Function(CG_2d, name="bathymetry")
        chk.load(bathymetry_2d)
        File(outputdir+"/bath.pvd").write(bathymetry_2d)
        chk.close()
    outfile = File(outputdir+"/depth.pvd")

    # viscosity
    with timed_stage('initialising viscosity'):
        chk = DumbCheckpoint(inputdir+"/viscosity", mode=FILE_READ)
        h_viscosity = Function(CG_2d, name="viscosity")
        chk.load(h_viscosity)
        File(outputdir+"/viscosity.pvd").write(h_viscosity)
        chk.close()

    # manning
    if manning == True:
        with timed_stage('initialising viscosity'):
            chk = DumbCheckpoint(inputdir+"/manning2D", mode=FILE_READ)
            n_manning = Function(CG_2d, name="manning")
            chk.load(n_manning)
            File(outputdir+"/manning.pvd").write(n_manning)
            chk.close()

    if identifier != -1:

        DG_2d = FunctionSpace(mesh2d, 'DG', 1)

        # elevation
        with timed_stage('initialising elevation'):
            chk = DumbCheckpoint(inputdir + "/elevation" + str(identifier), mode=FILE_READ)
            elev_init = Function(DG_2d, name="elevation")
            chk.load(elev_init)
            File(outputdir + "/elevation_imported.pvd").write(elev_init)
            chk.close()

        # velocity
        with timed_stage('initialising velocity'):
            chk = DumbCheckpoint(inputdir + "/velocity" + str(identifier), mode=FILE_READ)
            V = VectorFunctionSpace(mesh2d, 'DG', 1)
            uv_init = Function(V, name="velocity")
            chk.load(uv_init)
            File(outputdir + "/velocity_imported.pvd").write(uv_init)
            chk.close()

        if manning == True:
            return bathymetry_2d, h_viscosity, elev_init, uv_init, n_manning
        else:
            return bathymetry_2d, h_viscosity, elev_init, uv_init,
    else:
        if manning == True:
            return bathymetry_2d, h_viscosity, n_manning
        else:
            return bathymetry_2d, h_viscosity


def export_final_state(inputdir, identifier, uv, elev, lagoon = None):

    print_output("Exporting fields for subsequent simulation")

    chk = DumbCheckpoint(inputdir + "/velocity" + str(identifier + 1), mode=FILE_CREATE)
    chk.store(uv, name="velocity")
    File('outputs/velocityout.pvd').write(uv)
    chk.close()
    chk = DumbCheckpoint(inputdir + "/elevation" + str(identifier + 1), mode=FILE_CREATE)
    chk.store(elev, name="elevation")
    File('outputs/elevationout.pvd').write(elev)
    chk.close()

    if lagoon is not None:
        output_status = []
        for i in range(len(lagoon)):
            output_status.append(lagoon[i])
        pickle.dump(output_status, open(inputdir+"/barrage_status_"+str(identifier+1)+".p", "wb"))