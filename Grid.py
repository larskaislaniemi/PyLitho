import copy

class Grid:

    def __init__(self, nx=0, xs=[]):
        if len(xs) != nx:
            raise Exception("Grid.__init__(): len(xs) != nx")
        self.nx = nx
        self.xs = list(xs)


class GriddedData:

    def __init__(self, nx=0, xs=[]):
        self.metadata = {}
        self.data = {}
        self.grid = Grid(nx, xs)

    def addData(self, dataname, data, datapoint=-1, overwrite=False):
        if len(data) != self.grid.nx:
            raise Exception("GriddedData.addData(): Data size mismatches grid size: len(data) != self.grid.nx")
        if not overwrite:
            if dataname in self.data.keys():
                raise Exception("GriddedData.addData(): Data exists. Use overwrite=True to overwrite")
        if datapoint < 0:
            self.data[dataname] = copy.deepcopy(data)
        else:
            self.data[dataname][datapoint] = copy.deepcopy(data)

    def addMetaData(self, dataname, data, overwrite=False):
        if not overwrite:
            if dataname in self.metadata.keys():
                raise Exception("GriddedData.addMetaData(): Meta data exists. Use overwrite=True to overwrite")
        self.metadata[dataname] = copy.deepcopy(data)

