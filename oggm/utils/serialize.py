"""Serialization of OGGM data objects to and from plain numpy arrays.

Historically OGGM stored flowlines, centerlines and geometries by pickling
the Python objects directly. That ties the files to the exact class layout of
OGGM and of shapely at the time they were written: shapely 2.1 refuses to
unpickle geometries written by shapely < 2.0 (see #1524).

This module removes that coupling. An object is encoded as a flat dict of
plain numpy arrays plus a small JSON-able ``meta`` dict describing the
structure, and rebuilt from those on read. Nothing in the result refers to an
OGGM or a shapely class, so the on-disk data no longer depends on either
library's internals.

Which container those arrays end up in is deliberately not this module's
concern - see :py:meth:`oggm.GlacierDirectory.read_npz` and
:py:meth:`~oggm.GlacierDirectory.write_npz`.

To support a new data product, write an ``_encode_*`` / ``_decode_*`` pair and
add one entry to ``_CODECS``. Anything not listed there falls through to a
generic mapping codec which handles dicts and lists of dicts of arrays.
"""

import numpy as np
import shapely.geometry as shpg


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------

def _put(arrays, key, value):
    """Store a value as an array, skipping None."""
    if value is None:
        return False
    arrays[key] = np.asarray(value)
    return True


def _get(arrays, key):
    """Read an array back, or None if it was never stored."""
    return arrays.get(key)


def _scalar(arrays, key):
    """Read a stored value back as a plain Python scalar (not a 0-d array)."""
    val = arrays.get(key)
    if val is None:
        return None
    val = np.asarray(val)
    return val.item() if val.ndim == 0 else val


def _encode_linestring(line):
    """A LineString as an (n, 2) coordinate array."""
    if line is None:
        return None
    return np.asarray(line.coords, dtype=np.float64)


def _decode_linestring(coords):
    if coords is None:
        return None
    return shpg.LineString(np.asarray(coords))


def _encode_point(point):
    if point is None:
        return None
    return np.asarray(point.coords, dtype=np.float64).ravel()


def _decode_point(coords):
    if coords is None:
        return None
    return shpg.Point(np.asarray(coords).ravel())


def _extract_polygon_rings(geometry):
    """One ``(coords, poly_idx, is_exterior)`` tuple per ring.

    Polygons may contain interior holes and MultiPolygons several parts, each
    with its own holes. Rings are emitted exterior-first within each part.
    """
    rings = []
    if geometry.geom_type == 'Polygon':
        rings.append((np.asarray(geometry.exterior.coords), 0, True))
        for interior in geometry.interiors:
            rings.append((np.asarray(interior.coords), 0, False))
    elif geometry.geom_type == 'MultiPolygon':
        for poly_idx, part in enumerate(geometry.geoms):
            rings.append((np.asarray(part.exterior.coords), poly_idx, True))
            for interior in part.interiors:
                rings.append((np.asarray(interior.coords), poly_idx, False))
    else:
        raise ValueError('Unhandled geometry type: '
                         + repr(geometry.geom_type))
    return rings


def _encode_polygon(arrays, prefix, geometry):
    """Flatten a (Multi)Polygon into one vertex array plus ring indices."""
    rings = _extract_polygon_rings(geometry)
    coords = (np.concatenate([r[0] for r in rings], axis=0) if rings
              else np.zeros((0, 2), dtype=np.float64))
    arrays[prefix + 'vertices'] = np.asarray(coords, dtype=np.float64)
    arrays[prefix + 'ring_lengths'] = np.array([len(r[0]) for r in rings],
                                               dtype=np.int64)
    arrays[prefix + 'ring_poly_idx'] = np.array([r[1] for r in rings],
                                                dtype=np.int64)
    arrays[prefix + 'ring_is_exterior'] = np.array([r[2] for r in rings],
                                                   dtype=bool)


def _decode_polygon(arrays, prefix):
    coords = np.asarray(arrays[prefix + 'vertices'])
    ring_lengths = np.asarray(arrays[prefix + 'ring_lengths'])
    ring_poly_idx = np.asarray(arrays[prefix + 'ring_poly_idx'])
    ring_is_exterior = np.asarray(arrays[prefix + 'ring_is_exterior'])

    splits = np.cumsum(ring_lengths)[:-1]
    rings = np.split(coords, splits) if len(ring_lengths) else []

    parts = {}
    for ring, poly_idx, is_ext in zip(rings, ring_poly_idx, ring_is_exterior):
        part = parts.setdefault(int(poly_idx), {'exterior': None, 'holes': []})
        if bool(is_ext):
            part['exterior'] = ring
        else:
            part['holes'].append(ring)

    polygons = [shpg.Polygon(parts[i]['exterior'], parts[i]['holes'])
                for i in sorted(parts)]
    if len(polygons) == 1:
        return polygons[0]
    return shpg.MultiPolygon(polygons)


def _encode_index_list(arrays, prefix, index_list):
    """Flatten a list of (n, 2) index arrays (e.g. catchment_indices)."""
    parts = [np.asarray(a, dtype=np.int64).reshape(-1, 2) for a in index_list]
    coords = (np.concatenate(parts, axis=0) if parts
              else np.zeros((0, 2), dtype=np.int64))
    arrays[prefix + 'vertices'] = coords
    arrays[prefix + 'lengths'] = np.array([len(a) for a in parts],
                                          dtype=np.int64)


def _decode_index_list(arrays, prefix):
    coords = np.asarray(arrays[prefix + 'vertices'], dtype=np.int64)
    lengths = np.asarray(arrays[prefix + 'lengths'], dtype=np.int64)
    if not len(lengths):
        return []
    splits = np.cumsum(lengths)[:-1]
    return [a.reshape(-1, 2) for a in np.split(coords, splits)]


def _encode_mapping(arrays, prefix, mapping):
    """Encode a plain dict of arrays / scalars / nested dicts.

    Returns a spec recording, per key, whether the value was an array, a
    scalar or None, so that decoding restores the original Python types.
    """
    spec = {}
    for key, val in mapping.items():
        if val is None:
            spec[key] = {'kind': 'none'}
        elif isinstance(val, dict):
            spec[key] = {'kind': 'dict',
                         'spec': _encode_mapping(arrays, f'{prefix}{key}/',
                                                 val)}
        else:
            arr = np.asarray(val)
            arrays[prefix + key] = arr
            spec[key] = {'kind': 'array' if arr.ndim else 'scalar'}
    return spec


def _decode_mapping(arrays, prefix, spec):
    out = {}
    for key, entry in spec.items():
        kind = entry['kind']
        if kind == 'none':
            out[key] = None
        elif kind == 'dict':
            out[key] = _decode_mapping(arrays, f'{prefix}{key}/',
                                       entry['spec'])
        elif kind == 'scalar':
            out[key] = _scalar(arrays, prefix + key)
        else:
            out[key] = np.asarray(arrays[prefix + key])
    return out


def _encode_topology(arrays, objects):
    """Record who flows into whom, and where the junction sits.

    ``flows_to_point`` is stored rather than recomputed on read: it depends on
    how ``set_flows_to`` was originally called (``compute_centerlines`` passes
    ``check_tail=False``), so replaying the projection would not always give
    the same point back.
    """
    idx_of = {id(obj): i for i, obj in enumerate(objects)}
    flows_to = []
    for i, obj in enumerate(objects):
        target = getattr(obj, 'flows_to', None)
        flows_to.append(idx_of.get(id(target), -1) if target is not None
                        else -1)
        _put(arrays, f'{i}/flows_to_point',
             _encode_point(getattr(obj, 'flows_to_point', None)))
    return flows_to


def _decode_topology(arrays, objects, flows_to):
    """Rebuild flows_to / inflows without re-running the projection."""
    for i, target in enumerate(flows_to):
        if not 0 <= target < len(objects):
            continue
        obj, other = objects[i], objects[target]
        point = _decode_point(_get(arrays, f'{i}/flows_to_point'))
        if point is None:
            # Nothing stored (legacy data): fall back to recomputing it.
            obj.set_flows_to(other)
            continue
        obj.flows_to = other
        obj.flows_to_point = point
        other.inflow_points.append(point)
        other.inflows.append(obj)



# ----------------------------------------------------------------------
# Centerlines (also used for inversion_flowlines)
# ----------------------------------------------------------------------

# Attributes which cannot be passed to Centerline.__init__ and are restored
# by assignment afterwards.
_CENTERLINE_ATTRS = ['order', '_widths', 'is_rectangular', 'is_trapezoid',
                     'apparent_mb', 'flux', 'flux_out',
                     'flux_needs_correction', 'orig_centerline_id']

# Which of those are scalars rather than arrays.
_CENTERLINE_SCALARS = {'order', 'flux_out', 'flux_needs_correction',
                       'orig_centerline_id'}


def _encode_geometrical_widths(arrays, prefix, widths):
    """Ragged MultiLineString widths -> flat coords + two index arrays."""
    all_coords = []      # flat coords across every member line
    line_lengths = []    # vertices per member LineString
    width_line_counts = []   # member lines per width
    for width in widths:
        if width is None:
            members = []
        elif hasattr(width, 'geoms'):
            members = list(width.geoms)
        else:
            members = [width]
        count = 0
        for member in members:
            if member is None or member.is_empty:
                continue
            coords = np.asarray(member.coords, dtype=np.float64)
            all_coords.append(coords)
            line_lengths.append(len(coords))
            count += 1
        width_line_counts.append(count)

    coords = (np.concatenate(all_coords, axis=0) if all_coords
              else np.zeros((0, 2), dtype=np.float64))
    arrays[prefix + 'gw_vertices'] = coords
    arrays[prefix + 'gw_line_lengths'] = np.asarray(line_lengths,
                                                    dtype=np.int64)
    arrays[prefix + 'gw_width_counts'] = np.asarray(width_line_counts,
                                                    dtype=np.int64)


def _decode_geometrical_widths(arrays, prefix):
    coords = _get(arrays, prefix + 'gw_vertices')
    if coords is None:
        return None
    coords = np.asarray(coords)
    line_lengths = np.asarray(arrays[prefix + 'gw_line_lengths'],
                              dtype=np.int64)
    counts = np.asarray(arrays[prefix + 'gw_width_counts'], dtype=np.int64)
    splits = np.cumsum(line_lengths)[:-1]
    lines = np.split(coords, splits) if len(line_lengths) else []
    widths = []
    idx = 0
    for count in counts:
        members = [shpg.LineString(lines[idx + j]) for j in range(count)]
        idx += count
        widths.append(shpg.MultiLineString(members))
    return widths


def _encode_centerlines(centerlines, gdir=None):
    from oggm import Centerline

    if not isinstance(centerlines, list):
        raise TypeError('Expected a list of Centerline objects, got '
                        f'{type(centerlines).__name__}.')
    if not all(isinstance(cl, Centerline) for cl in centerlines):
        raise TypeError('All items must be Centerline instances.')

    arrays = {}
    for i, cl in enumerate(centerlines):
        prefix = f'{i}/'
        _put(arrays, prefix + 'line', _encode_linestring(cl.line))
        _put(arrays, prefix + 'orig_head',
             _encode_point(getattr(cl, 'orig_head', None)))
        _put(arrays, prefix + 'dx', cl.dx)
        _put(arrays, prefix + 'map_dx', getattr(cl, 'map_dx', None))
        _put(arrays, prefix + 'surface_h', cl.surface_h)
        _put(arrays, prefix + 'rgi_id', getattr(cl, 'rgi_id', None))
        for attr in _CENTERLINE_ATTRS:
            _put(arrays, prefix + attr, getattr(cl, attr, None))
        widths = getattr(cl, 'geometrical_widths', None)
        if widths is not None:
            _encode_geometrical_widths(arrays, prefix, widths)

    meta = {'kind': 'centerlines', 'n': len(centerlines),
            'flows_to': _encode_topology(arrays, centerlines)}
    return arrays, meta


def _decode_centerlines(arrays, meta, gdir=None):
    from oggm import Centerline

    centerlines = []
    for i in range(meta['n']):
        prefix = f'{i}/'
        centerline = Centerline(
            line=_decode_linestring(_get(arrays, prefix + 'line')),
            dx=_scalar(arrays, prefix + 'dx'),
            surface_h=_get(arrays, prefix + 'surface_h'),
            orig_head=_decode_point(_get(arrays, prefix + 'orig_head')),
            rgi_id=_scalar(arrays, prefix + 'rgi_id'),
            map_dx=_scalar(arrays, prefix + 'map_dx'),
        )
        for attr in _CENTERLINE_ATTRS:
            if prefix + attr not in arrays:
                continue
            if attr in _CENTERLINE_SCALARS:
                setattr(centerline, attr, _scalar(arrays, prefix + attr))
            else:
                setattr(centerline, attr, np.asarray(arrays[prefix + attr]))
        widths = _decode_geometrical_widths(arrays, prefix)
        if widths is not None:
            centerline.geometrical_widths = widths
        centerlines.append(centerline)

    _decode_topology(arrays, centerlines, meta.get('flows_to', []))

    return centerlines


# ----------------------------------------------------------------------
# Model flowlines
# ----------------------------------------------------------------------

# Cached arrays MixedBedFlowline computes in __init__; stored and restored
# verbatim so that a round trip is bit-identical.
_MIXED_CACHED = ['_sqrt_bed', '_w0_m']


def _encode_model_flowlines(flowlines, gdir=None):
    from oggm.core.flowline import (Flowline, MixedBedFlowline,
                                    ParabolicBedFlowline,
                                    RectangularBedFlowline,
                                    TrapezoidalBedFlowline)

    if not isinstance(flowlines, list):
        raise TypeError('Expected a list of Flowline objects, got '
                        f'{type(flowlines).__name__}.')
    if not all(isinstance(fl, Flowline) for fl in flowlines):
        raise TypeError('All items must be Flowline instances.')

    arrays = {}
    classes = []
    for i, fl in enumerate(flowlines):
        prefix = f'{i}/'
        _put(arrays, prefix + 'line', _encode_linestring(fl.line))
        _put(arrays, prefix + 'dx', fl.dx)
        _put(arrays, prefix + 'map_dx', fl.map_dx)
        _put(arrays, prefix + 'surface_h', fl.surface_h)
        _put(arrays, prefix + 'bed_h', fl.bed_h)
        _put(arrays, prefix + 'rgi_id', getattr(fl, 'rgi_id', None))
        _put(arrays, prefix + 'water_level', getattr(fl, 'water_level', None))
        _put(arrays, prefix + 'order', getattr(fl, 'order', None))

        if isinstance(fl, MixedBedFlowline):
            lambdas = getattr(fl, 'lambdas', None)
            _put(arrays, prefix + 'section', fl.section)
            _put(arrays, prefix + 'bed_shape', fl.bed_shape)
            _put(arrays, prefix + 'is_trapezoid', fl.is_trapezoid)
            _put(arrays, prefix + 'widths_m', fl.widths_m)
            _put(arrays, prefix + 'lambdas',
                 getattr(fl, '_lambdas', None) if lambdas is None else lambdas)
            for attr in _MIXED_CACHED:
                _put(arrays, prefix + attr, getattr(fl, attr, None))
        elif isinstance(fl, ParabolicBedFlowline):
            _put(arrays, prefix + 'bed_shape', fl.bed_shape)
        elif isinstance(fl, RectangularBedFlowline):
            _put(arrays, prefix + 'widths', getattr(fl, '_widths', None))
        elif isinstance(fl, TrapezoidalBedFlowline):
            # Trapezoid rebuilds _w0_m from widths and lambdas, so store the
            # width property (= widths_m / map_dx) rather than _w0_m.
            _put(arrays, prefix + 'widths', fl.widths)
            _put(arrays, prefix + 'lambdas', getattr(fl, '_lambdas', None))

        classes.append(type(fl).__name__)

    meta = {'kind': 'model_flowlines', 'n': len(flowlines),
            'classes': classes,
            'flows_to': _encode_topology(arrays, flowlines)}
    if gdir is not None:
        # Flowline.__init__ resolves min_ice_thick_for_length and
        # glacier_length_method from gdir.settings, so the suffix those came
        # from has to travel with the data or the read would silently fall
        # back to cfg.PARAMS.
        meta['settings_filesuffix'] = gdir.settings_filesuffix
    return arrays, meta


def _decode_model_flowlines(arrays, meta, gdir=None):
    from oggm.core.flowline import (MixedBedFlowline, ParabolicBedFlowline,
                                    RectangularBedFlowline,
                                    TrapezoidalBedFlowline)

    settings_filesuffix = meta.get('settings_filesuffix', '')
    flowlines = []
    for i in range(meta['n']):
        prefix = f'{i}/'
        # Passing gdir through means the flowline picks up map_trafo and the
        # per-gdir settings exactly as it would when freshly built.
        base = dict(
            line=_decode_linestring(_get(arrays, prefix + 'line')),
            dx=_scalar(arrays, prefix + 'dx'),
            map_dx=_scalar(arrays, prefix + 'map_dx'),
            surface_h=_get(arrays, prefix + 'surface_h'),
            bed_h=_get(arrays, prefix + 'bed_h'),
            rgi_id=_scalar(arrays, prefix + 'rgi_id'),
            water_level=_scalar(arrays, prefix + 'water_level'),
            gdir=gdir,
            settings_filesuffix=settings_filesuffix,
        )
        cls_name = meta['classes'][i]
        if cls_name == 'ParabolicBedFlowline':
            flowline = ParabolicBedFlowline(
                bed_shape=_get(arrays, prefix + 'bed_shape'), **base)
        elif cls_name == 'RectangularBedFlowline':
            flowline = RectangularBedFlowline(
                widths=_get(arrays, prefix + 'widths'), **base)
        elif cls_name == 'TrapezoidalBedFlowline':
            flowline = TrapezoidalBedFlowline(
                widths=_get(arrays, prefix + 'widths'),
                lambdas=_get(arrays, prefix + 'lambdas'), **base)
        else:
            flowline = MixedBedFlowline(
                section=_get(arrays, prefix + 'section'),
                bed_shape=_get(arrays, prefix + 'bed_shape'),
                is_trapezoid=_get(arrays, prefix + 'is_trapezoid'),
                lambdas=_get(arrays, prefix + 'lambdas'),
                widths_m=_get(arrays, prefix + 'widths_m'), **base)
            for attr in _MIXED_CACHED:
                if prefix + attr in arrays:
                    setattr(flowline, attr,
                            np.asarray(arrays[prefix + attr]))
        flowline.order = _scalar(arrays, prefix + 'order')
        flowlines.append(flowline)

    _decode_topology(arrays, flowlines, meta.get('flows_to', []))

    return flowlines


# ----------------------------------------------------------------------
# Geometries and downstream line
# ----------------------------------------------------------------------

def _encode_geometries(geometries, gdir=None):
    arrays, spec = {}, {}
    for name, value in geometries.items():
        if 'polygon_hr' in name or 'polygon_pix' in name:
            _encode_polygon(arrays, f'{name}/', value)
            spec[name] = {'kind': 'polygon'}
        elif 'catchment_indices' in name:
            _encode_index_list(arrays, f'{name}/', value)
            spec[name] = {'kind': 'index_list'}
        elif 'downstream_line' in name:
            if _put(arrays, name, _encode_linestring(value)):
                spec[name] = {'kind': 'linestring'}
            else:
                spec[name] = {'kind': 'none'}
        else:
            spec.update(_encode_mapping(arrays, '', {name: value}))
    return arrays, {'kind': 'geometries', 'spec': spec}


def _decode_geometries(arrays, meta, gdir=None):
    out = {}
    for name, entry in meta['spec'].items():
        kind = entry['kind']
        if kind == 'polygon':
            out[name] = _decode_polygon(arrays, f'{name}/')
        elif kind == 'index_list':
            out[name] = _decode_index_list(arrays, f'{name}/')
        elif kind == 'linestring':
            out[name] = _decode_linestring(_get(arrays, name))
        else:
            out.update(_decode_mapping(arrays, '', {name: entry}))
    return out


def _encode_downstream_line(data, gdir=None):
    if not isinstance(data, dict):
        raise TypeError('downstream_line must be a dict, got '
                        f'{type(data).__name__}.')
    arrays, spec = {}, {}
    for name, value in data.items():
        if name in ('downstream_line', 'full_line'):
            if _put(arrays, name, _encode_linestring(value)):
                spec[name] = {'kind': 'linestring'}
            else:
                spec[name] = {'kind': 'none'}
        else:
            spec.update(_encode_mapping(arrays, '', {name: value}))
    return arrays, {'kind': 'downstream_line', 'spec': spec}


def _decode_downstream_line(arrays, meta, gdir=None):
    out = {}
    for name, entry in meta['spec'].items():
        if entry['kind'] == 'linestring':
            out[name] = _decode_linestring(_get(arrays, name))
        else:
            out.update(_decode_mapping(arrays, '', {name: entry}))
    # compute_downstream_line may not write full_line; callers expect the key
    out.setdefault('full_line', None)
    return out


# ----------------------------------------------------------------------
# Generic fallback: a dict, or a list of dicts (inversion_input/output, ...)
# ----------------------------------------------------------------------

def _encode_generic(obj, gdir=None):
    arrays = {}
    if isinstance(obj, dict):
        spec = _encode_mapping(arrays, '', obj)
        return arrays, {'kind': 'dict', 'spec': spec}
    if isinstance(obj, list) and all(isinstance(i, dict) for i in obj):
        specs = [_encode_mapping(arrays, f'{i}/', item)
                 for i, item in enumerate(obj)]
        return arrays, {'kind': 'list_of_dicts', 'n': len(obj),
                        'specs': specs}
    raise TypeError(
        f'Cannot serialize {type(obj).__name__}: expected a dict, a list of '
        'dicts, or a registered OGGM type. Add a codec to '
        'oggm.utils.serialize if this type needs storing.')


def _decode_generic(arrays, meta, gdir=None):
    if meta['kind'] == 'dict':
        return _decode_mapping(arrays, '', meta['spec'])
    return [_decode_mapping(arrays, f'{i}/', spec)
            for i, spec in enumerate(meta['specs'])]


# ----------------------------------------------------------------------
# Registry and public API
# ----------------------------------------------------------------------

_CODECS = {
    'centerlines': (_encode_centerlines, _decode_centerlines),
    'inversion_flowlines': (_encode_centerlines, _decode_centerlines),
    'model_flowlines': (_encode_model_flowlines, _decode_model_flowlines),
    'geometries': (_encode_geometries, _decode_geometries),
    'downstream_line': (_encode_downstream_line, _decode_downstream_line),
}


def encode(basename, obj, gdir=None):
    """Encode an OGGM object as plain arrays plus a structure description.

    Parameters
    ----------
    basename : str
        The cfg.BASENAMES key, without any filesuffix.
    obj : object
        The object to encode.
    gdir : :py:class:`oggm.GlacierDirectory`, optional
        The directory the object belongs to. Used to record context that the
        object itself does not carry, such as which settings file the
        flowlines were built against.

    Returns
    -------
    (dict, dict)
        A dict of ``{name: np.ndarray}`` and a JSON-able ``meta`` dict.
    """
    codec = _CODECS.get(basename)
    if codec is None:
        return _encode_generic(obj, gdir=gdir)
    return codec[0](obj, gdir=gdir)


def decode(basename, arrays, meta, gdir=None):
    """Rebuild an OGGM object from arrays written by :py:func:`encode`.

    Parameters
    ----------
    basename : str
        The cfg.BASENAMES key, without any filesuffix.
    arrays : dict
        The arrays as stored.
    meta : dict
        The structure description returned by :py:func:`encode`.
    gdir : :py:class:`oggm.GlacierDirectory`, optional
        Passed to flowline constructors so that they pick up the map
        transform and the per-gdir settings, exactly as they would when
        freshly built.
    """
    codec = _CODECS.get(basename)
    if codec is None:
        return _decode_generic(arrays, meta, gdir=gdir)
    return codec[1](arrays, meta, gdir=gdir)
