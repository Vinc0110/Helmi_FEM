import sys
import gmsh


def mesh(plane: str = 'h'):
    # WR6 standard gain horn (Pasternack PEWAN1028)
    # https://www.pasternack.com/images/ProductPDF/PEWAN1028.pdf

    if plane.lower() == 'e':
        h_waveg = 0.83
        h_horn = 13
    elif plane.lower() == 'h':
        h_waveg = 1.65
        h_horn = 17.5
    else:
        raise ValueError(f'Invalid parameter {plane}.')

    unit = 1e-3
    r_freespace = 20
    l_waveg = 4.11
    l_horn = 40.85

    wavelen = 3e8 / 200e9 / unit
    meshsize = wavelen / 8

    gmsh.initialize()
    gmsh.model.add('Waveguide horn')

    tag_pt_feed_ymin = gmsh.model.geo.addPoint(-1 * l_waveg - l_horn, -0.5 * h_waveg, 0, meshSize=meshsize)
    tag_pt_feed_ymax = gmsh.model.geo.addPoint(-1 * l_waveg - l_horn, 0.5 * h_waveg, 0, meshSize=meshsize)
    tag_pt_feed_y0 = gmsh.model.geo.addPoint(-1 * l_waveg - l_horn, 0, 0, meshSize=meshsize)
    tag_pt_hornfeed_ymin = gmsh.model.geo.addPoint(-1 * l_horn, -0.5 * h_waveg, 0, meshSize=meshsize)
    tag_pt_hornfeed_ymax = gmsh.model.geo.addPoint(-1 * l_horn, 0.5 * h_waveg, 0, meshSize=meshsize)
    tag_pt_hornfeed_y0 = gmsh.model.geo.addPoint(-1 * l_horn, 0, 0, meshSize=meshsize)
    tag_pt_horn_ymin = gmsh.model.geo.addPoint(0, -0.5 * h_horn, 0, meshSize=meshsize)
    tag_pt_horn_ymax = gmsh.model.geo.addPoint(0, 0.5 * h_horn, 0, meshSize=meshsize)
    tag_pt_horn_y0 = gmsh.model.geo.addPoint(0, 0, 0, meshSize=meshsize)
    tag_pt_arc_ymin = gmsh.model.geo.addPoint(0, -1 * r_freespace, 0, meshSize=meshsize)
    tag_pt_arc_ymax = gmsh.model.geo.addPoint(0, r_freespace, 0, meshSize=meshsize)
    tag_pt_arc_y0 = gmsh.model.geo.addPoint(r_freespace, 0, 0, meshSize=meshsize)

    tag_line_feed = gmsh.model.geo.addLine(tag_pt_feed_ymin, tag_pt_feed_ymax)
    tag_line_waveg_ymax = gmsh.model.geo.addLine(tag_pt_feed_ymax, tag_pt_hornfeed_ymax)
    tag_line_horn_ymax = gmsh.model.geo.addLine(tag_pt_hornfeed_ymax, tag_pt_horn_ymax)
    tag_line_freespace_ymax = gmsh.model.geo.addLine(tag_pt_horn_ymax, tag_pt_arc_ymax)
    tag_arc_freespace_ymax = gmsh.model.geo.addCircleArc(tag_pt_arc_ymax, tag_pt_horn_y0, tag_pt_arc_y0)
    tag_arc_freespace_ymin = gmsh.model.geo.addCircleArc(tag_pt_arc_y0, tag_pt_horn_y0, tag_pt_arc_ymin)
    tag_line_freespace_ymin = gmsh.model.geo.addLine(tag_pt_arc_ymin, tag_pt_horn_ymin)
    tag_line_horn_ymin = gmsh.model.geo.addLine(tag_pt_horn_ymin, tag_pt_hornfeed_ymin)
    tag_line_waveg_ymin = gmsh.model.geo.addLine(tag_pt_hornfeed_ymin, tag_pt_feed_ymin)

    tag_loop = gmsh.model.geo.addCurveLoop([tag_line_feed, tag_line_waveg_ymax, tag_line_horn_ymax,
                                            tag_line_freespace_ymax, tag_arc_freespace_ymax, tag_arc_freespace_ymin,
                                            tag_line_freespace_ymin, tag_line_horn_ymin, tag_line_waveg_ymin])

    tag_surf = gmsh.model.geo.addPlaneSurface([tag_loop])

    gmsh.model.geo.addPhysicalGroup(2, [tag_surf], name='air')
    gmsh.model.geo.addPhysicalGroup(1, [tag_line_feed], name='bound_feed')
    gmsh.model.geo.addPhysicalGroup(1, [tag_line_freespace_ymax, tag_arc_freespace_ymax,
                                        tag_arc_freespace_ymin, tag_line_freespace_ymin], name='bound_freespace')
    gmsh.model.geo.addPhysicalGroup(1, [tag_line_waveg_ymax, tag_line_horn_ymax], name='bound_ymax')
    gmsh.model.geo.addPhysicalGroup(1, [tag_line_waveg_ymin, tag_line_horn_ymin], name='bound_ymin')

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    #gmsh.fltk.run()
    gmsh.write(f'./horn_{plane}-plane.msh')
    gmsh.finalize()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        plane = sys.argv[1]
    else:
        plane = 'e'
    mesh(plane)
