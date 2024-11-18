from functools import cached_property

class Dimension(tuple[int, ...]):
    """ Dimensions of the ??? space """
    @cached_property
    def symmetries(self) -> tuple[int, ...]:
        """ Returns length of the symmetries in the dimensions """
        from .utils import group_by_block
        return tuple(length for _, length in group_by_block(self))

    @cached_property
    def sum(self) -> int:
        return sum(self)