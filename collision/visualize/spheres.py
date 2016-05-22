# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Vispy Development Team.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.
# -----------------------------------------------------------------------------
"""
Spheres Visual and shader definitions.

Raytraces a sphere at each vertex
"""

import numpy as np

from vispy.color import ColorArray
from vispy.gloo import VertexBuffer
from vispy.visuals.visual import Visual
from vispy import scene


vert = """
uniform float pointRadius;  // point size in world space
uniform float pointScale;   // scale to calculate size in pixels

attribute vec3  a_position;
attribute vec4  a_color;

varying vec4 v_color;

void main()
{
    // calculate window-space point size
    vec3 posEye = vec3($transform_scene(vec4(a_position, 1.0)));
    float dist = length(posEye);

    gl_PointSize = pointRadius * pointScale / dist;
    v_color = a_color;

    gl_Position = $transform(vec4(a_position, 1.0));
}
"""


frag = """
varying vec4 v_color;

void main()
{
    const vec3 lightDir = vec3(0.577, 0.577, 0.577);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_PointCoord.xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(N.xy, N.xy);

    if (mag > 1.0) discard;   // kill pixels outside circle

    N.z = sqrt(1.0-mag);

    // calculate lighting
    float diffuse = max(0.5, dot(lightDir, N));

    gl_FragColor = v_color * diffuse;
    gl_FragColor.a = 1;

    gl_FragDepth = gl_FragCoord.z - N.z * 0.00005f;
}
"""

class SpheresVisual(Visual):
    """ Visual displaying raytraced spheres.
    """
    def __init__(self, radius, **kwargs):
        self._vbo = VertexBuffer()
        # self._v_size_var = Variable('varying float v_size')
        self._data = None
        Visual.__init__(self, vcode=vert, fcode=frag)

        # m_window_h / tanf(m_fov * 0.5f * (float)M_PI / 180.0)
        self.view_program["pointScale"] = 600 / np.tan(60 * 0.5 * np.pi / 180)
        self.view_program["pointRadius"] = radius

        self.set_gl_state(depth_test=True, blend=False)
        self._draw_mode = 'points'
        if len(kwargs) > 0:
            self.set_data(**kwargs)
        self.freeze()

    def set_data(self, pos=None, size=10., scaling=False):
        """ Set the data used to display this visual.

        Parameters
        ----------
        pos : array
            The array of locations to display each symbol.
        size : float or array
            The symbol size in px.
        scaling : bool
            If set to True, marker scales when rezooming.
        """
        assert (isinstance(pos, np.ndarray) and
                pos.ndim == 2 and pos.shape[1] in (2, 3))

        if self._data is None:
            n = len(pos)
            data = np.zeros(n, dtype=[('a_position', np.float32, 3),
                                      ('a_color', np.float32, 4),
                                    ])
            data['a_position'][:, :pos.shape[1]] = pos

            data['a_color'] = 1
            data['a_color'][:, :3] = np.random.rand(n, 3)
            self._data = data
        else:
            self._data['a_position'][:, :pos.shape[1]] = pos

        # self.shared_program['u_antialias'] = self.antialias  # XXX make prop
        self._vbo.set_data(self._data)
        self.shared_program.bind(self._vbo)
        self.update()

    def _prepare_transforms(self, view):
        xform = view.transforms.get_transform()
        view.view_program.vert['transform'] = xform
        xform = view.transforms.get_transform(map_to='render')
        view.view_program.vert['transform_scene'] = xform

    def _prepare_draw(self, view):
        pass

    def _compute_bounds(self, axis, view):
        pos = self._data['a_position']
        if pos is None:
            return None
        if pos.shape[1] > axis:
            return (pos[:, axis].min(), pos[:, axis].max())
        else:
            return (0, 0)


Spheres = scene.visuals.create_visual_node(SpheresVisual)
