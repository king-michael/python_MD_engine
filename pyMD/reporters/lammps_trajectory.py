from ..misc.lammps_dump_custom import DumpFileWriter
from .base import Reporter
import numpy as np

class LammpsTrajectoryReporter(Reporter):
    default_fields = ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    map_dtype = {'%d': np.int64, '%g': np.float64, '%f': np.float64, '*': np.float64}

    def __init__(self, filename, n_dump=1000, fields=('id', 'x', 'y', 'z'), triclinic=False):

        if fields != ('id', 'x', 'y', 'z'):
            raise NotImplementedError('Other fields are not implemented yet')

        super().__init__(n_dump)

        self.filename = filename
        self.fields = fields
        self.triclinic = triclinic

        self.fp = DumpFileWriter(filename,
                                 fields=fields,
                                 triclinic=triclinic)

    def _map_fromat2dtype(self, fmt):
        return self.map_dtype.get(fmt, self.map_dtype['*'])

    def setup(self):
        # dtype = np.dtype(list(zip(self.fp.fields, map(_map_fromat2dtype, sel.fp.format_str))))

        self.n_atoms = len(self.engine.positions)
        self.output_array = np.zeros((self.n_atoms, len(self.fields)), dtype=np.float64)

        tmp = []
        for i, f in enumerate(('x', 'y', 'z')):
            try:
                c = self.fields.index(f)
                tmp.append((i, c))
            except:
                pass
        if len(tmp) > 0:
            self.write_coordinates = True
            self._axis_inp_coords, self._axis_out_coords = zip(*tmp)
        else:
            self.write_coordinates = False

        if 'id' in self.fields:
            self.output_array[:, self.fields.index('id')] = np.arange(self.n_atoms) + 1

    def report(self, step):
        # handle coordinates
        if self.write_coordinates:
            self.output_array[:,self._axis_out_coords] = self.engine.positions[:, self._axis_inp_coords]

        # Handle box
        (lx, _, _), (xy, ly, _), (xz, yz, lz) = self.engine.box
        xlo, xhi = 0.0, lx
        ylo, yhi = 0.0, ly
        zlo, zhi = 0.0, lz
        if self.triclinic:
            box = [[xlo, xhi, xy], [ylo, yhi, xz], [zlo, zhi, yz]]
        else:
            box = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]

        # write out
        self.fp.write(self.output_array, box, ts=step)

    def close(self):
        self.fp.close()

    def __del__(self):
        self.close()
