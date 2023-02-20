import gmsh
import sys


def mesh(n_slits, w_slit, pitch):
    l_slit = 2.5
    r_box = 0.5 * pitch * (n_slits + 1)

    mesh_slit = 0.25
    mesh_box = 0.5

    if w_slit > pitch:
        raise ValueError('slit width > pitch')

    if 0.5 * n_slits * pitch > r_box:
        raise ValueError('box too small')

    gmsh.initialize()
    gmsh.model.add('Slit')

    tag_pt_center = gmsh.model.geo.addPoint(0, 0, 0)

    tag_pt_box_00 = gmsh.model.geo.addPoint(-0.5 * r_box, -1 * r_box, 0, meshSize=mesh_box)
    tag_pt_box_01 = gmsh.model.geo.addPoint(-0.5 * r_box, r_box, 0, meshSize=mesh_box)
    tag_pt_box_0c = gmsh.model.geo.addPoint(-0.5 * r_box, 0, 0, meshSize=mesh_box)
    tag_pt_box_c01 = gmsh.model.geo.addPoint(-0.5 * l_slit, r_box, 0, meshSize=mesh_box)
    tag_pt_box_c00 = gmsh.model.geo.addPoint(-0.5 * l_slit, -1 * r_box, 0, meshSize=mesh_box)
    tag_pt_box_c11 = gmsh.model.geo.addPoint(0.5 * l_slit, r_box, 0, meshSize=mesh_box)
    tag_pt_box_c10 = gmsh.model.geo.addPoint(0.5 * l_slit, -1 * r_box, 0, meshSize=mesh_box)
    tag_pt_box_1c = gmsh.model.geo.addPoint(r_box, 0, 0, meshSize=mesh_box)

    tag_arc_box_11 = gmsh.model.geo.addCircleArc(tag_pt_box_1c, tag_pt_center, tag_pt_box_c11)
    tag_arc_box_10 = gmsh.model.geo.addCircleArc(tag_pt_box_c10, tag_pt_center, tag_pt_box_1c)
    tag_line_box_c1 = gmsh.model.geo.addLine(tag_pt_box_c11, tag_pt_box_c01)
    tag_line_box_c0 = gmsh.model.geo.addLine(tag_pt_box_c00, tag_pt_box_c10)
    tag_line_box_01 = gmsh.model.geo.addLine(tag_pt_box_c01, tag_pt_box_01)
    tag_line_box_0c = gmsh.model.geo.addLine(tag_pt_box_01, tag_pt_box_00)
    tag_line_box_00 = gmsh.model.geo.addLine(tag_pt_box_00, tag_pt_box_c00)

    tag_loop_box = gmsh.model.geo.addCurveLoop([tag_arc_box_11, tag_line_box_c1, tag_line_box_01, tag_line_box_0c,
                                                tag_line_box_00, tag_line_box_c0, tag_arc_box_10])

    # slitted aperture (blind)
    height_blinds = pitch - w_slit
    height_aperture = n_slits * w_slit + (n_slits - 1) * height_blinds
    tag_pt_slit_01 = gmsh.model.geo.addPoint(-0.5 * l_slit, 0.5 * height_aperture, 0, meshSize=mesh_slit)
    tag_pt_slit_11 = gmsh.model.geo.addPoint(0.5 * l_slit, 0.5 * height_aperture, 0, meshSize=mesh_slit)
    tag_pt_slit_00 = gmsh.model.geo.addPoint(-0.5 * l_slit, -0.5 * height_aperture, 0, meshSize=mesh_slit)
    tag_pt_slit_10 = gmsh.model.geo.addPoint(0.5 * l_slit, -0.5 * height_aperture, 0, meshSize=mesh_slit)

    tag_line_blind_01 = gmsh.model.geo.addLine(tag_pt_box_c01, tag_pt_slit_01)
    tag_line_blind_h1 = gmsh.model.geo.addLine(tag_pt_slit_01, tag_pt_slit_11)
    tag_line_blind_11 = gmsh.model.geo.addLine(tag_pt_slit_11, tag_pt_box_c11)
    tag_loop_blind_1 = gmsh.model.geo.addCurveLoop([tag_line_box_c1, tag_line_blind_01, tag_line_blind_h1, tag_line_blind_11])

    tag_line_blind_00 = gmsh.model.geo.addLine(tag_pt_box_c10, tag_pt_slit_10)
    tag_line_blind_h0 = gmsh.model.geo.addLine(tag_pt_slit_10, tag_pt_slit_00)
    tag_line_blind_10 = gmsh.model.geo.addLine(tag_pt_slit_00, tag_pt_box_c00)
    tag_loop_blind_0 = gmsh.model.geo.addCurveLoop([tag_line_box_c0, tag_line_blind_00, tag_line_blind_h0, tag_line_blind_10])

    tags_loops_blind = [tag_loop_blind_0, tag_loop_blind_1]
    tags_lines_blind = [tag_line_blind_00, tag_line_blind_h0, tag_line_blind_10,
                        tag_line_blind_01, tag_line_blind_h1, tag_line_blind_11]

    for i in range(n_slits - 1):
        y_bot = -0.5 * height_aperture + i * pitch + w_slit
        y_top = y_bot + height_blinds
        tag_pt_01 = gmsh.model.geo.addPoint(-0.5 * l_slit, y_top, 0, meshSize=mesh_slit)
        tag_pt_00 = gmsh.model.geo.addPoint(-0.5 * l_slit, y_bot, 0, meshSize=mesh_slit)
        tag_pt_10 = gmsh.model.geo.addPoint(0.5 * l_slit, y_bot, 0, meshSize=mesh_slit)
        tag_pt_11 = gmsh.model.geo.addPoint(0.5 * l_slit, y_top, 0, meshSize=mesh_slit)
        tag_line_0c = gmsh.model.geo.addLine(tag_pt_01, tag_pt_00)
        tag_line_c0 = gmsh.model.geo.addLine(tag_pt_00, tag_pt_10)
        tag_line_1c = gmsh.model.geo.addLine(tag_pt_10, tag_pt_11)
        tag_line_c1 = gmsh.model.geo.addLine(tag_pt_11, tag_pt_01)
        tag_loop = gmsh.model.geo.addCurveLoop([tag_line_0c, tag_line_c0, tag_line_1c, tag_line_c1])
        tags_loops_blind.append(tag_loop)
        tags_lines_blind.append(tag_line_0c)
        tags_lines_blind.append(tag_line_c0)
        tags_lines_blind.append(tag_line_1c)
        tags_lines_blind.append(tag_line_c1)

    tag_surf_air = gmsh.model.geo.addPlaneSurface([tag_loop_box, *tags_loops_blind])

    gmsh.model.geo.addPhysicalGroup(2, [tag_surf_air], name='air')
    gmsh.model.geo.addPhysicalGroup(1, [tag_line_box_0c], name='bounds_left')
    gmsh.model.geo.addPhysicalGroup(1, [tag_arc_box_10, tag_arc_box_11], name='bounds_right')
    gmsh.model.geo.addPhysicalGroup(1, tags_lines_blind, name='bounds_aperture')

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    #gmsh.fltk.run()
    gmsh.write('./slit.msh')
    gmsh.finalize()


if __name__ == '__main__':
    n_slits = 3
    w_slit = 10
    pitch = 20
    mesh(n_slits, w_slit, pitch)
