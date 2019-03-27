"""
Various tools for internal use.
"""

__license__ = """
    Copyright 2015 Stephane De Mita, Mathieu Siol

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

import tempfile, os

########################################################################

class TempFile:

    """
    Each instance of this class creates a temporary file that is
    guaranteed to be deleted at destruction time or when :meth:`~.clean`
    is called. The file is closed immediately after creation. The is no
    error if the file is deleted before destruction of the corresponding
    instance.
    """

    def __init__(self):
        fd, self._fname = tempfile.mkstemp()
        os.close(fd)

    def __del__(self):
        self.clean()

    def clean(self):

        """
        Remove the underlying temporary file if it exists.
        """

        if os.path.isfile(self._fname):
            os.remove(self._fname) 

    @property
    def fname(self):

        """
        Name of the temporary file.
        """

        return self._fname

########################################################################

class Pool:

    """
    Minimal object pool implementation. This class is aimed to manage
    object recycling. Objects are yielded by :meth:`~.Pool.get` and returned
    after use by :meth:`~.put`. ``len(pool)`` returns the number of
    objects currently store in the pool (although you normally don't
    need to worry about it). Similarly, one can clear the stored objects
    with :meth:`~.Pool.clear`.

    :param generator: method allowing to create new objects.
    :param arguments: arguments to be passed to the generator method.
        Use ``None`` or an empty sequence to pass no arguments.
    :param recycle: method to recycle an object returning to the pool.
        This method should be take this type of objects as unique
        argument. If ``None``, no recycling method is applied.
    :param recycle_args: arguments to be passed to the recycling method.
        Use ``None`` or an empty sequence to pass no arguments (which is
        the default).

    Here is a demonstration with :class:`~.Align` objects:

    .. code-block:: python

       import egglib

       P = egglib.misc.Pool(egglib.Align, (), egglib.Align.reset)
       print 'pool size:', len(P)
       aln1 = P.get()
       aln2 = P.get()
       print 'pool size:', len(P)
       print 'fresh objects:', aln1.ns, aln1.ls, aln2.ns, aln2.ls
       aln1.add_sample('', 'AAAAAAAAAA')
       aln1.add_sample('', 'AAAAAAAAAA')
       aln1.add_sample('', 'AAAAAAAAAA')
       aln2.add_sample('', 'AAAAAAAAAAAAAAAAAAAA')
       aln2.add_sample('', 'AAAAAAAAAAAAAAAAAAAA')
       print 'modified objects:', aln1.ns, aln1.ls, aln2.ns, aln2.ls
       P.put(aln1)
       P.put(aln2)
       print 'objects back to pool:', aln1.ns, aln1.ls, aln2.ns, aln2.ls
       print 'one is not supposed to keep a reference of recycled objects!'
       print 'pool size:', len(P)
       aln3 = P.get()
       print 'recycled object:', aln3.ns, aln3.ls

    This produces the following output::

        pool size: 0
        pool size: 0
        fresh objects: 0 0 0 0
        modified objects: 3 10 2 20
        objects back to pool: 0 0 0 0
        one is not supposed to keep a reference of recycled objects!
        pool size: 2
        recycled object: 0 0
    """

    def __init__(self, generator, arguments, recycle, recycle_args = None):
        self._pool = []
        self._g = generator
        if arguments == None: self._args = ()
        else: self._args = arguments
        if recycle == None: self._rec = (lambda obj: None)
        else: self._rec = recycle
        if recycle_args == None: self._rec_args = ()
        else: self._rec_args = recycle_args


    def get(self):

        """
        Get an object from the pool.
        """

        if len(self._pool) == 0: self._pool.append(self._g(* self._args))
        return self._pool.pop()

    def put(self, obj):

        """
        Return an object to the pool.
        """

        self._rec(obj, * self._rec_args)
        self._pool.append(obj)

    def __len__(self):
        return len(self._pool)

    def clear(self):
        del self._pool[:]

########################################################################

def _get_struct(struct, nsam, nout):

    # (1) get a Structure or None, (2) check type, (3) check that max
    # indexes are <nsam and/or nout, (4) check that has (a) pop
    # structure only, (b) indiv structure only, or (c) both with pop
    # obligatory mapped to indiv, (5) return a (struct, has_genotypes,
    # has_structure) tuple where struct is a _eggwrapper Structure
    # instance or None.

    if struct == None:
        has_genotypes = False
        has_structure = False
    else:
        try:
            has_structure = struct._obj.pop_flag()
            has_genotypes = struct._obj.indiv_flag()
        except AttributeError:
            raise TypeError, 'invalid type: {0} (a `Structure` is expected)'.format(type(struct))

        if (has_structure == False or struct._obj.num_pop() < 2) and has_genotypes == False:
            raise ValueError, 'at least one (either individual or population) level of structure must be defined'
        if has_structure and has_genotypes and struct._obj.pop_to_indiv() == False:
            raise ValueError, 'if both individual and population levels are defined, individuals are required to map to populations'

        max_ingroup = struct.max_ingroup
        if max_ingroup != None and max_ingroup >= nsam: raise ValueError, 'index defined in structure instance out of range'
        max_outgroup = struct.max_outgroup
        if max_outgroup != None and max_outgroup >= nout: raise ValueError, 'outgroup index defined in structure instance out of range'

        struct = struct._obj

    return struct, has_genotypes, has_structure
