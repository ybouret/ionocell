""" IonoCell """


class Grid:
    """ Grid for one space """

    def __init__(self, _length):
        self.length = _length

    def show(self):
        """ show Grid parameters """
        print("length=", self.length)


grid = Grid(0.1)

grid.show()
