import gmsh


def mesh():
    r_core = 10 / 2
    r_cladding = 125 / 2
    r_box = 200 / 2

    mesh_fiber = 2
    mesh_box = 4

    gmsh.initialize()
    gmsh.model.add('Scatterer')

    tag_circle_center = gmsh.model.geo.addPoint(0, 0, 0, meshSize=mesh_fiber)

    tag_circle_core_x0 = gmsh.model.geo.addPoint(-1 * r_core, 0, 0, meshSize=mesh_fiber)
    tag_circle_core_x1 = gmsh.model.geo.addPoint(r_core, 0, 0, meshSize=mesh_fiber)
    tag_circle_core_y0 = gmsh.model.geo.addPoint(0, -1 * r_core, 0, meshSize=mesh_fiber)
    tag_circle_core_y1 = gmsh.model.geo.addPoint(0, r_core, 0, meshSize=mesh_fiber)

    tag_circle_cladding_x0 = gmsh.model.geo.addPoint(-1 * r_cladding, 0, 0, meshSize=mesh_fiber)
    tag_circle_cladding_x1 = gmsh.model.geo.addPoint(r_cladding, 0, 0, meshSize=mesh_fiber)
    tag_circle_cladding_y0 = gmsh.model.geo.addPoint(0, -1 * r_cladding, 0, meshSize=mesh_fiber)
    tag_circle_cladding_y1 = gmsh.model.geo.addPoint(0, r_cladding, 0, meshSize=mesh_fiber)

    tag_circle_box_x0 = gmsh.model.geo.addPoint(-1 * r_box, 0, 0, meshSize=mesh_box)
    tag_circle_box_x1 = gmsh.model.geo.addPoint(r_box, 0, 0, meshSize=mesh_box)
    tag_circle_box_y0 = gmsh.model.geo.addPoint(0, -1 * r_box, 0, meshSize=mesh_box)
    tag_circle_box_y1 = gmsh.model.geo.addPoint(0, r_box, 0, meshSize=mesh_box)

    tag_arc_core0 = gmsh.model.geo.addCircleArc(tag_circle_core_x1, tag_circle_center, tag_circle_core_y1)
    tag_arc_core1 = gmsh.model.geo.addCircleArc(tag_circle_core_y1, tag_circle_center, tag_circle_core_x0)
    tag_arc_core2 = gmsh.model.geo.addCircleArc(tag_circle_core_x0, tag_circle_center, tag_circle_core_y0)
    tag_arc_core3 = gmsh.model.geo.addCircleArc(tag_circle_core_y0, tag_circle_center, tag_circle_core_x1)

    tag_arc_cladding0 = gmsh.model.geo.addCircleArc(tag_circle_cladding_x1, tag_circle_center, tag_circle_cladding_y1)
    tag_arc_cladding1 = gmsh.model.geo.addCircleArc(tag_circle_cladding_y1, tag_circle_center, tag_circle_cladding_x0)
    tag_arc_cladding2 = gmsh.model.geo.addCircleArc(tag_circle_cladding_x0, tag_circle_center, tag_circle_cladding_y0)
    tag_arc_cladding3 = gmsh.model.geo.addCircleArc(tag_circle_cladding_y0, tag_circle_center, tag_circle_cladding_x1)

    tag_arc_box0 = gmsh.model.geo.addCircleArc(tag_circle_box_x1, tag_circle_center, tag_circle_box_y1)
    tag_arc_box1 = gmsh.model.geo.addCircleArc(tag_circle_box_y1, tag_circle_center, tag_circle_box_x0)
    tag_arc_box2 = gmsh.model.geo.addCircleArc(tag_circle_box_x0, tag_circle_center, tag_circle_box_y0)
    tag_arc_box3 = gmsh.model.geo.addCircleArc(tag_circle_box_y0, tag_circle_center, tag_circle_box_x1)

    tag_loop_box = gmsh.model.geo.addCurveLoop([tag_arc_box0, tag_arc_box1, tag_arc_box2, tag_arc_box3])
    tag_loop_core = gmsh.model.geo.addCurveLoop([tag_arc_core0, tag_arc_core1, tag_arc_core2, tag_arc_core3])
    tag_loop_cladding = gmsh.model.geo.addCurveLoop([tag_arc_cladding0, tag_arc_cladding1, tag_arc_cladding2, tag_arc_cladding3])

    tag_surf_air = gmsh.model.geo.addPlaneSurface([tag_loop_box, tag_loop_core, tag_loop_cladding])
    tag_surf_cladding = gmsh.model.geo.addPlaneSurface([tag_loop_cladding, tag_loop_core])
    tag_surf_core = gmsh.model.geo.addPlaneSurface([tag_loop_core])

    gmsh.model.geo.addPhysicalGroup(2, [tag_surf_air], name='air')
    gmsh.model.geo.addPhysicalGroup(2, [tag_surf_core], name='core')
    gmsh.model.geo.addPhysicalGroup(2, [tag_surf_cladding], name='cladding')
    gmsh.model.geo.addPhysicalGroup(1, [tag_arc_box0, tag_arc_box1, tag_arc_box2, tag_arc_box3], name='bound')
    gmsh.model.geo.addPhysicalGroup(1, [tag_arc_core0, tag_arc_core1, tag_arc_core2, tag_arc_core3], name='interface_core')
    gmsh.model.geo.addPhysicalGroup(1, [tag_arc_cladding0, tag_arc_cladding1, tag_arc_cladding2, tag_arc_cladding3], name='interface_cladding')

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    #gmsh.fltk.run()
    gmsh.write('./scatterer.msh')
    gmsh.finalize()


if __name__ == '__main__':
    mesh()
