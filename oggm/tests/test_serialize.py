"""Unit tests for oggm.utils.serialize.

These exercise the codecs directly, without going through a glacier
directory. The integration side - gdir.read_npz / write_npz, the legacy
pickle path, and round trips over real data - lives in
``test_workflow.TestNpzStore``.
"""
import numpy as np
import pytest
import shapely.geometry as shpg

from oggm import cfg
from oggm.utils import serialize


@pytest.fixture(autouse=True)
def init_cfg():
    cfg.initialize_minimal()


def roundtrip(basename, obj, gdir=None):
    """Encode and decode, simulating the npz trip through plain arrays."""
    arrays, meta = serialize.encode(basename, obj, gdir=gdir)
    # everything handed to np.savez must be a plain array
    assert all(isinstance(v, np.ndarray) for v in arrays.values())
    assert not any(v.dtype == object for v in arrays.values()), \
        'object arrays would need allow_pickle=True to read back'
    arrays = {k: np.asarray(v) for k, v in arrays.items()}
    return serialize.decode(basename, arrays, meta, gdir=gdir)


class TestGeometryHelpers:

    def test_linestring(self):
        line = shpg.LineString([(0, 0), (1, 2), (3, 4)])
        back = serialize._decode_linestring(serialize._encode_linestring(line))
        assert back.equals_exact(line, 1e-9)
        assert serialize._encode_linestring(None) is None
        assert serialize._decode_linestring(None) is None

    def test_point(self):
        point = shpg.Point(3.5, -1.25)
        back = serialize._decode_point(serialize._encode_point(point))
        assert back.equals_exact(point, 1e-9)
        assert serialize._encode_point(None) is None

    def test_polygon_with_interiors(self):
        shell = [(0, 0), (10, 0), (10, 10), (0, 10)]
        hole1 = [(1, 1), (2, 1), (2, 2), (1, 2)]
        hole2 = [(5, 5), (6, 5), (6, 6), (5, 6)]
        poly = shpg.Polygon(shell, [hole1, hole2])

        arrays = {}
        serialize._encode_polygon(arrays, 'p/', poly)
        back = serialize._decode_polygon(arrays, 'p/')

        assert back.equals_exact(poly, 1e-9)
        assert len(back.interiors) == 2

    def test_multipolygon(self):
        a = shpg.Polygon([(0, 0), (2, 0), (2, 2), (0, 2)],
                         [[(0.5, 0.5), (1, 0.5), (1, 1), (0.5, 1)]])
        b = shpg.Polygon([(5, 5), (7, 5), (7, 7), (5, 7)])
        multi = shpg.MultiPolygon([a, b])

        arrays = {}
        serialize._encode_polygon(arrays, 'p/', multi)
        back = serialize._decode_polygon(arrays, 'p/')

        assert back.geom_type == 'MultiPolygon'
        assert len(back.geoms) == 2
        assert back.equals_exact(multi, 1e-9)
        assert len(back.geoms[0].interiors) == 1

    def test_polygon_rejects_other_geometries(self):
        with pytest.raises(ValueError, match='Unhandled geometry type'):
            serialize._extract_polygon_rings(shpg.LineString([(0, 0), (1, 1)]))

    def test_index_list(self):
        idx = [np.array([[1, 2], [3, 4]]), np.array([[5, 6]]),
               np.zeros((0, 2), dtype=np.int64)]
        arrays = {}
        serialize._encode_index_list(arrays, 'c/', idx)
        back = serialize._decode_index_list(arrays, 'c/')

        assert len(back) == 3
        for a, b in zip(idx, back):
            np.testing.assert_array_equal(a, b)

    def test_index_list_empty(self):
        arrays = {}
        serialize._encode_index_list(arrays, 'c/', [])
        assert serialize._decode_index_list(arrays, 'c/') == []

    def test_ragged_geometrical_widths(self):
        widths = [
            shpg.MultiLineString([[(0, 0), (1, 1)]]),
            shpg.MultiLineString([[(0, 0), (1, 1), (2, 2)],
                                  [(5, 5), (6, 6)]]),
            shpg.MultiLineString([]),
        ]
        arrays = {}
        serialize._encode_geometrical_widths(arrays, 'f/', widths)
        back = serialize._decode_geometrical_widths(arrays, 'f/')

        assert len(back) == 3
        assert len(back[0].geoms) == 1
        assert len(back[1].geoms) == 2
        assert len(back[2].geoms) == 0
        assert back[1].equals_exact(widths[1], 1e-9)


class TestGenericCodec:

    def test_dict_preserves_types(self):
        var = {'arr': np.array([1.0, 2.0]), 'scalar': 3.0, 'flag': True,
               'text': 'hello', 'missing': None, 'count': 7}
        back = roundtrip('some_unregistered_name', var)

        np.testing.assert_array_equal(back['arr'], [1.0, 2.0])
        assert back['scalar'] == 3.0 and isinstance(back['scalar'], float)
        assert back['flag'] is True
        assert back['text'] == 'hello' and isinstance(back['text'], str)
        assert back['missing'] is None
        assert back['count'] == 7 and isinstance(back['count'], int)

    def test_nested_dict(self):
        var = {'outer': 1.0, 'inner': {'a': np.arange(3), 'b': None}}
        back = roundtrip('some_unregistered_name', var)
        assert back['outer'] == 1.0
        np.testing.assert_array_equal(back['inner']['a'], [0, 1, 2])
        assert back['inner']['b'] is None

    def test_list_of_dicts(self):
        var = [{'flux': np.array([1.0]), 'is_last': False},
               {'flux': np.array([2.0, 3.0]), 'is_last': True}]
        back = roundtrip('inversion_output', var)

        assert len(back) == 2
        assert back[0]['is_last'] is False
        assert back[1]['is_last'] is True
        np.testing.assert_array_equal(back[1]['flux'], [2.0, 3.0])

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError, match='Cannot serialize'):
            serialize.encode('whatever', object())

    def test_list_of_non_dicts_raises(self):
        with pytest.raises(TypeError, match='Cannot serialize'):
            serialize.encode('whatever', [1, 2, 3])


class TestCenterlineCodec:

    @staticmethod
    def _make(n=5):
        from oggm import Centerline
        lines = []
        for k in range(2):
            coords = [(float(i), float(i + k * 10)) for i in range(n)]
            cl = Centerline(shpg.LineString(coords), dx=2.0,
                            surface_h=np.linspace(3000, 2500, n),
                            rgi_id='RGI60-11.00897', map_dx=100.)
            cl.order = k
            cl.flux = np.ones(n)
            cl.flux_out = 12.5
            cl.is_rectangular = np.zeros(n, dtype=bool)
            cl.apparent_mb = np.linspace(-1, 1, n)
            lines.append(cl)
        lines[0].set_flows_to(lines[1])
        return lines

    def test_roundtrip(self):
        orig = self._make()
        back = roundtrip('inversion_flowlines', orig)

        assert len(back) == len(orig)
        for a, b in zip(orig, back):
            assert a.line.equals_exact(b.line, 1e-9)
            np.testing.assert_allclose(a.surface_h, b.surface_h)
            np.testing.assert_allclose(a.flux, b.flux)
            np.testing.assert_allclose(a.apparent_mb, b.apparent_mb)
            np.testing.assert_array_equal(a.is_rectangular, b.is_rectangular)

    def test_scalars_are_python_scalars(self):
        back = roundtrip('inversion_flowlines', self._make())
        fl = back[0]
        assert isinstance(fl.rgi_id, str)
        assert 'RGI60' in fl.rgi_id
        assert isinstance(fl.dx, float)
        assert isinstance(fl.map_dx, float)
        assert isinstance(fl.order, int)
        assert isinstance(fl.flux_out, float)

    def test_topology_and_junction_point(self):
        orig = self._make()
        back = roundtrip('inversion_flowlines', orig)

        assert back[0].flows_to is back[1]
        assert back[1].flows_to is None
        assert back[1].inflows == [back[0]]
        # the junction point is stored, not recomputed
        assert back[0].flows_to_point.equals_exact(orig[0].flows_to_point,
                                                   1e-9)

    def test_dx_may_be_none(self):
        """compute_centerlines leaves dx unset; it must stay None."""
        from oggm import Centerline
        cl = Centerline(shpg.LineString([(0, 0), (1, 1), (2, 2)]))
        back = roundtrip('centerlines', [cl])
        assert back[0].dx is None

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError, match='Centerline'):
            serialize.encode('inversion_flowlines', [object()])


class TestModelFlowlineCodec:

    @staticmethod
    def _make(n=10):
        from oggm.core.flowline import (MixedBedFlowline,
                                        ParabolicBedFlowline,
                                        RectangularBedFlowline,
                                        TrapezoidalBedFlowline)
        surface_h = np.linspace(3000, 2500, n)
        bed_h = surface_h - 100.
        line = shpg.LineString([(float(i), 0.) for i in range(n)])
        common = dict(line=line, dx=1., map_dx=100., surface_h=surface_h,
                      bed_h=bed_h, rgi_id='RGI60-11.00897')
        # a consistent parabolic section, so MixedBedFlowline accepts it
        bed_shape = np.ones(n) * 0.003
        thick = surface_h - bed_h
        widths_m = np.sqrt(4 * thick / bed_shape)
        section = 2. / 3. * widths_m * thick

        return [
            ParabolicBedFlowline(bed_shape=bed_shape, **common),
            RectangularBedFlowline(widths=np.ones(n) * 3., **common),
            TrapezoidalBedFlowline(widths=np.ones(n) * 3.,
                                   lambdas=np.ones(n), **common),
            MixedBedFlowline(section=section, bed_shape=bed_shape,
                             is_trapezoid=np.zeros(n, dtype=bool),
                             lambdas=np.zeros(n), widths_m=widths_m,
                             **common),
        ]

    def test_every_subclass_roundtrips(self):
        orig = self._make()
        back = roundtrip('model_flowlines', orig)

        assert [type(f).__name__ for f in back] == \
               [type(f).__name__ for f in orig]
        for a, b in zip(orig, back):
            np.testing.assert_allclose(a.surface_h, b.surface_h)
            np.testing.assert_allclose(a.bed_h, b.bed_h)
            np.testing.assert_allclose(a.widths_m, b.widths_m)
            np.testing.assert_allclose(a.section, b.section)
            np.testing.assert_allclose(a.volume_m3, b.volume_m3)

    def test_mixed_cached_arrays_restored_verbatim(self):
        orig = self._make()
        back = roundtrip('model_flowlines', orig)
        mixed_a, mixed_b = orig[-1], back[-1]
        np.testing.assert_array_equal(mixed_a._w0_m, mixed_b._w0_m)
        np.testing.assert_array_equal(mixed_a._sqrt_bed, mixed_b._sqrt_bed)

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError, match='Flowline'):
            serialize.encode('model_flowlines', [object()])


class TestGeometriesCodec:

    def test_roundtrip(self):
        poly = shpg.Polygon([(0, 0), (10, 0), (10, 10), (0, 10)],
                            [[(1, 1), (2, 1), (2, 2), (1, 2)]])
        geom = {
            'polygon_hr': poly,
            'polygon_pix': shpg.Polygon([(0, 0), (5, 0), (5, 5)]),
            'polygon_area': 100.0,
            'catchment_indices': [np.array([[1, 2], [3, 4]])],
        }
        back = roundtrip('geometries', geom)

        assert back['polygon_hr'].equals_exact(poly, 1e-9)
        assert len(back['polygon_hr'].interiors) == 1
        assert back['polygon_area'] == 100.0
        assert isinstance(back['polygon_area'], float)
        np.testing.assert_array_equal(back['catchment_indices'][0],
                                      [[1, 2], [3, 4]])


class TestDownstreamLineCodec:

    def test_roundtrip_with_full_line(self):
        data = {'downstream_line': shpg.LineString([(0, 0), (1, 1)]),
                'full_line': shpg.LineString([(0, 0), (1, 1), (2, 2), (3, 3)]),
                'bedshapes': np.array([0.1, 0.2])}
        back = roundtrip('downstream_line', data)

        assert back['downstream_line'].equals_exact(data['downstream_line'],
                                                    1e-9)
        assert back['full_line'].equals_exact(data['full_line'], 1e-9)
        np.testing.assert_allclose(back['bedshapes'], [0.1, 0.2])

    def test_full_line_may_be_none(self):
        data = {'downstream_line': shpg.LineString([(0, 0), (1, 1)]),
                'full_line': None}
        back = roundtrip('downstream_line', data)
        assert back['full_line'] is None

    def test_full_line_key_always_present(self):
        data = {'downstream_line': shpg.LineString([(0, 0), (1, 1)])}
        back = roundtrip('downstream_line', data)
        assert 'full_line' in back and back['full_line'] is None

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError, match='must be a dict'):
            serialize.encode('downstream_line', [1, 2])
