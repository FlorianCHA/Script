"""
This module contains the classes :class:`~.Tree` and its helper
:class:`~.Node`.
"""

__license__ = """
    Copyright 2008-2011,2015 Stephane De Mita, Mathieu Siol

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

import _eggwrapper
import misc
import sys

_random_pool = misc.Pool(_eggwrapper.Random, (), None)

########################################################################

class Node(object):

    """
    This class provides an interface to a :class:`~.Tree` instance's
    nodes and allows access and modification of data attached to a given
    node as well as the tree descending from that node. A node must be
    understood as the point *below* a branch. Branches (connections
    between nodes) have a direction: they go *from* a node *to* another
    node. Nodes have therefore *children*  and *parents* (a given node
    can have one, or possibly no, parent). Connecting a node to itself,
    making a two-way branch (to branches connecting the same two nodes
    in opposite directions) or duplicate branches (between the same two
    nodes and in the same direction) are illegal.

    :param label: node label (in case of a terminal node, its leaf
        label), if needed. Labels, if provided, are expected to be
        strings for terminal nodes, and numeric values for internal
        nodes, but technically, all user-supplied values are accepted
        (however, some :class:`~.Tree` methods require proper types).
    """

    ####################################################################

    def __init__(self, label=None):
        self._parent = None
        self._parent_branch = None
        self._children = []
        self._child_branches = {}
        self._label = label

    ####################################################################

    def newick(self, labels=True, brlens=True):

        """
        Formats the node and the subtree descending from is as a
        newick string. If *labels* is ``False``, omit internal branch
        labels. If *brlens* is ``False``, omit branch lengths.
        """

        return ''.join(self._newick(labels, brlens) + [';'])

    ####################################################################

    def _newick(self, labels, brlens):
        string_bits = []
        if len(self._children):
            string_bits.append('(')
            for child in self._children:
                string_bits.extend(child._newick(labels, brlens))
                brlen = self._child_branches[child]
                if brlen != None and brlens: string_bits.append(':' + str(brlen))
                string_bits.append(',')
            string_bits[-1] = ')'
        if self._label != None and (labels or not len(self._children)): string_bits.append(str(self._label))
        return string_bits

    ####################################################################

    def leaves_down(self):

        """
        Recursively gets all leaf labels descending from that node. If
        this is a terminal node, returns its label in a one-item list.
        """

        leaves = []
        if len(self._children) > 0:
            for child in self._children: leaves += child.leaves_down()
            return leaves
        else:
            return [self._label]

    ####################################################################
    
    def leaves_up(self):

        """
        Recursively gets all leaf labels contained on the other side of
        the tree (that is, all leaves except those descending from this
        node). If this is the root node, returns an empty list).
        """

        leaves = []
        if self._parent == None: return []
        else:
            for brother in self._parent._children:
                if brother != self: leaves += brother.leaves_down()
            leaves += self._parent.leaves_up()
        return leaves

    ####################################################################

    def _unlink(self):

        # Clears all references to other :class:`~.Node` instances
        # (parent and children) from this node. The nodes in question
        # are not affected (they are not deleted, and references they
        # hold to the current node are not erased).

        self._parent = None
        self._parent_branch = None
        self._children = []
        self._child_branches = {}

    ####################################################################

    def _add_child(self, label, brlen):

        # Generates a new :class:`~.Node` instance descending from the
        # current instance (add a child to this instance).

        # :param label: label of the new node.
        # :param brlen: length of the branch connecting this node to the
        #     new node.
        # :return: The new node.

        child = Node(label)
        child._parent = self
        child._parent_branch = brlen
        self._children.append(child)
        self._child_branches[child] = brlen
        return child

    ####################################################################

    def _connect(self, node, brlen):

        # Connect this node to an other, existing, node. The orientation
        # of the link is *from* the current instance *to* the passed
        # instance. If the passed node has children, they will be
        # retained, allowing to connect a subtree to this tree. It is
        # not allowed to connect a node that has a (non-``None`) parent.
        # To do so, one must first disconnect the node from its parent
        # and then connect it to a new one.

        # :param node: an existing node.
        # :param brlen: length of the branch between this node and the
        #     passed one.

        # tests are not performed for a hidden method
        # if node._parent != None: raise ValueError, 'cannot connect node: already has a parent (disconnect first)'
        self._children.append(node)
        self._child_branches[node] = brlen
        node._parent = self
        node._parent_branch = brlen
        
    ####################################################################

    def has_descendant(self, node):

        """
        Return a boolean indicating whether the :class:`~.Node` passed
        as argument can be found amongst the descendants (children,
        children of children, and so on down to the leaves) of this
        node.
        """

        for child in self._children:
            if child == node: return True
            if child.has_descendant(node): return True
        return False

    ####################################################################

    def is_parent(self, node):
        
        """
        Return ``True`` if the :class:`~.Node` instance passed as
        argument is the parent of the current node. As a side effect,
        passing ``None`` to an instance that has no parent will return
        ``True``.
        """
    
        return self._parent == node
    
    ####################################################################

    def is_child(self, node):

        """
        Return ``True`` if the :class:`~.Node` instance passed as
        argument is one of the children of the current node.
        """

        for child in self._children:
            if node == child: return True
        return False

    ####################################################################

    @property
    def parent(self):

        """
        Reference to the parent of the current node (``None`` if the
        node has no parent).
        """

        return self._parent

    ####################################################################

    def children(self):

        """
        Return an iterator over this node's children.

        .. versionchanged:: 3.0.0
            Returns an iterator instead of a list.
        """

        for node in self._children: yield node

    ####################################################################

    def child(self, idx):

        """
        Return a given child, as a :class:`~.Node` instance.
        """

        try: return self._children[idx]
        except IndexError: raise IndexError, 'node index out of range'

    ####################################################################

    def siblings(self):

        """
        List of other children of this node's parent. It is required
        that this node has a parent.
        """

        if self._parent == None: raise ValueError, 'cannot access siblings: this node has no parent'
        return [node for node in self._parent._children if node != self]

    ####################################################################

    def _remove_parent(self):

        # Removes the branch between this node and its parent. Note that
        # this method removes also this node from its parent's children
        # (both sides of the branch are removed).

        self._parent._children.remove(self)
        del self._parent._child_branches[self]
        self._parent = None
        self._parent_branch = None

    ####################################################################

    def _remove_child(self, node):

        # Removes the branch between this node and one of its children.
        # Note that this method also removes the parent of child in
        # question (both sides of the branch are removed).

        self._children.remove(node)
        del self._child_branches[node]
        node._parent = None
        node._parent_branch = None

    ####################################################################

    def set_branch_to(self, child, brlen):

        """
        Set the length of the branch to one of this node's children.

        :param child: node whose branch should be resized. It can be
            represented by a direct reference (as a :class:~.Node`), or
            by its index in the child list. In that case, ensure that
            the index is currently valid.

        :param brlen: new branch length (it is allowed to pass ``None``).
        """

        if isinstance(child, Node):
            if child not in self._children: raise ValueError, 'invalid node: not part of the list of children'
            self._child_branches[child] = brlen
            child._parent_branch = brlen
        else:
            if child >= len(self._children): raise ValueError, 'invalid node index'
            self._child_branches[self._children[child]] = brlen
            self._children[child]._parent_branch = brlen

    ####################################################################

    def branch_to(self, child):
        
        """
        Get the length of the branch to one of this node's children.
        Non-specified branch lengths are represented by ``None``.

        :param child: nodes whose branch should be returned. It can be
            represented by a direct reference (as a :class:~.Node`), or
            by its index in the children list. In that case, ensure that
            the index is currently valid.
        """

        if isinstance(child, Node):
            if child not in self._children: raise ValueError, 'invalid node: not part of the list of children'
            return self._child_branches[child]
        else:
            if child >= len(self._children): raise ValueError, 'invalid node index'
            return self._child_branches[self._children[child]]

    ####################################################################

    @property
    def parent_branch(self):
        
        """
        Length of the branch to this node's parent. Non-specified branch
        lengths are represented by ``None``. An exception is thrown if
        this node has no parent. This attribute can be modified.
        """

        if self._parent == None: raise ValueError, 'cannot get parent\'s branch (this node has no parent)'
        return self._parent_branch

    @parent_branch.setter
    def parent_branch(self, brlen):
        if self._parent == None: raise ValueError, 'cannot set parent\'s branch (this node has no parent)'
        self._parent.set_branch_to(self, brlen) # will set both brlen at both sides

    ####################################################################

    @property
    def label(self):
        """ Node's label (modifiable). """
        return self._label

    @label.setter
    def label(self, label):
        self._label = label

    ####################################################################

    @property
    def num_children(self):
        """ Number of children connected to this node. """
        return len(self._children)

########################################################################

class Tree(object):

    """
    This class Handle trees. A tree is a linked collection of nodes
    which all have one parent (except the ultimate base node) and any
    number of children.

    Nodes are implemented as :class:`~.Node` instances. A node without
    children is a *leaf* (others are *internal*). A node with exactly
    one child is generally meaningless, but is allowed. All nodes
    (internal nodes as well as leaves) have a *label* which in the case
    of leaves can be used as leaf name. It is not possible to apply a
    name *and* a label to leaf node, in agreement with the newick
    format. All connections between nodes (*branches*) are oriented and
    can have a length (although the lengths can be omitted) but note
    that labels are applied to nodes, not branches. All :class:`~.Tree`
    instances have at least one base node which is the only one allowed
    not to have a parent. Network-like structures are not allowed
    (because nodes must have exactly one parent).

    Import and export to/from strings and files are in the bracket-based
    `newick <http://evolution.genetics.washington.edu/phylip/newicktree.html>`_
    format (the parser treats terminal node labels as strings, and
    internal node labels as integers or, by default, floats).
    :class:`~.Tree` instances can be exported using the built-in
    :class:`str` function, and the method :meth:`~.Tree.newick`. Nodes
    also have a :meth:`~.Node.newick` method.

    :class:`~.Tree` instances are iterable. Three iterators are
    provided: one is depth-first (:meth:`~.depth_iter`), another is
    breath-first (:meth:`~.breadth_iter`) and one iterates on terminal
    nodes (leaves) only (:meth:`~.iter_leaves`).

    The instance can be initialized as an empty tree (with only a root
    node), or from a newick-formatted string. By default, the string is
    read from the file name passed as the *fname* argument to the
    constructor, but it can be passed directly through the constructor
    argument *string*. It is not allowed to set both *fname* and
    *string* at the same time. The newick parser expects a well-formed
    newick string (including the trailing semicolon).

    .. versionchanged:: 2.0.1
        (This change concerns the constructor.) Imports directly from a
        file. If a string is passed, it is interpreted as a file name by
        default.

    .. versionchanged:: 3.0.0
       Several interface changes. Trees are not allowed to be network at
       all anymore.
    """

    ####################################################################
    
    def __init__(self, fname=None, string=None):

        self._base = Node()
        self._leaves = []
        self._nodes = [self._base]

        if fname != None and string != None:
            raise ValueError, 'cannot use both arguments of Tree\'s constructor simulateneously'

        if fname != None:
            f = open(fname)
            string = f.read()
            f.close()

        if string != None:
            string = string.translate(None, ' \n\r\t')
            if len(string) == 0: raise ValueError, 'empty newick string'
            if string[0] != '(' or string[-2:]!= ');': raise ValueError, 'invalid newick string'
            self._parse(string[:-1], self._base)

    ####################################################################

    @staticmethod
    def _from_coalesce(coal, idx, fun):

        # Create a :class:`~.Tree` instance from the results of a
        # coalescent simulation.
        #
        # :param coal: a :cpp:class:`~.Coalesce` object which has run a
        #     valid simulation.
        # :param idx: the index of a tree (must be within range of
        #     :cpp:meth:`~.Coalesce.number_of_trees`).
        # :param fun: executable expression converting the label to a
        #     string.
        #
        # :return: A new :class:`~.Tree` instance.

        source = coal.tree(idx)
        root = source.root()
        tree = Tree.__new__(Tree)
        tree._base = Node()
        tree._leaves = []
        tree._nodes = [tree._base]
        tree._from_coalesce_helper(tree._base, source, root, fun)
        return tree

    #####                                                          #####

    def _from_coalesce_helper(self, node_dst, tree_src, node_src, fun):

        # Helper for _from_coalesce
        # :param node_dst: current node in self.
        # :param tree_src: the reference C++ tree.
        # :param node_src: current C++ node in reference tree
        # :param fun: like for :meth:`~._from_coalesce`.

        if node_src.is_terminal():
            self._leaves.append(node_dst)
            node_dst._label = fun(node_src.label())
        else:
            son1 = tree_src.node(node_src.son1())
            child1 = node_dst._add_child(None, son1.get_L())
            self._from_coalesce_helper(child1, tree_src, son1, fun)
            self._nodes.append(child1)

            son2 = tree_src.node(node_src.son2())
            child2 = node_dst._add_child(None, son2.get_L())
            self._from_coalesce_helper(child2, tree_src, son2, fun)
            self._nodes.append(child2)

    ####################################################################

    def copy(self, node=None):

        """
        Create a new instance of :class:`~.Tree` that is a deep copy of
        a subtree of the the current tree.

        :param node: a :class:`~.Node` instance (one of the nodes of
            the current tree) at the base of the subtree that should be
            copied. By default, or if ``None`` or the base of the tree
            is passed, the whole tree is copied. It is not allowed to
            pass a leaf.
        """

        if node == None or node == self._base: return self._copy_full()
        if node not in self._nodes: raise ValueError, 'cannot copy tree: node is not part of tree'
        if node in self._leaves: raise ValueError, 'cannot copy tree: node is terminal'

        clone = Tree.__new__(Tree)
        clone._base = Node()
        clone._nodes = [clone._base]
        clone._leaves = []
        self._copy_helper(clone, clone._base, node)

        return clone

    #####                                                          #####

    def _copy_helper(self, clone, dest_node, node):

        dest_node._label = node._label
        if node in self._leaves:
            clone._leaves.append(dest_node)
        for child in node._children:
            new_node = dest_node._add_child(child._label, child._parent_branch)
            self._copy_helper(clone, new_node, child)
            clone._nodes.append(new_node)

    #####                                                          #####

    def _copy_full(self):

        # create a new bare tree (all nodes are also disconnected)
        clone = Tree.__new__(Tree)
        clone._nodes = [Node() for i in xrange(len(self._nodes))]
        clone._leaves = []

        # locate the base by its index
        base_index = self._nodes.index(self._base)
        clone._base = clone._nodes[base_index]

        # connect nodes bases on their indexes
        for node, copy in zip(self._nodes, clone._nodes):
            copy._label = node._label                  # copy label
            copy._parent_branch = node._parent_branch    # copy parent brlen
            if node._parent == None:
                copy._parent = None
            else:
                idx = self._nodes.index(node._parent)
                copy._parent = clone._nodes[idx]       # copy parent
            for child in node._children:
                idx = self._nodes.index(child)
                child_copy = clone._nodes[idx]
                copy._children.append(child_copy)     # copy a child
                copy._child_branches[child_copy] = node._child_branches[child] # copy child brlen
            if node in self._leaves: clone._leaves.append(copy)

        return clone

    ####################################################################

    def extract(self, node, label=None):

        """
        Remove a subtree of the current tree and return it as a new
        instance of :class:`~.Tree`. All nodes of the subtree descending
        from the requested node will now belong to the new tree. The
        label of the node at the base of the selected clade is deleted.
        In the original tree, the extracted clade is replaced by a
        terminal node which has, by default, the label of the node at
        the base of the extracted clade, or the label passed as
        argument.

        :param node: a :class:`~.Node` instance (one of the nodes of
            the current tree) at the base of the subtree that should be
            extracted. It is not allowed to pass a leaf or the base of
            the tree.

        :param label: label to affect to the terminal node that is
            introduced to replace the extracted clade in the orginal
            tree. By default (or if ``None``), use the label of the node
            passed as first argument (note that this label should be in
            principle a number).
        """

        if node not in self._nodes: raise ValueError, 'cannot extract subtree: node is not part of tree'
        if node in self._leaves: raise ValueError, 'cannot extract subtree: node is terminal'
        if node == self._base: raise ValueError, 'cannot extract subtree: node is base'

        parent = node._parent
        brlen = node._parent_branch
        node._remove_parent()
        if label == None: label = node._label
        parent._add_child(label, brlen)

        clone = Tree.__new__(Tree)
        clone._base = node
        clone._nodes = []
        clone._leaves = []
        self._extract_helper(clone, node)
        return clone

    #####                                                          #####

    def _extract_helper(self, clone, node):
        clone._nodes.append(node)
        self._nodes.remove(node)
        if node in self._leaves:
            clone._leaves.append(node)
            self._leaves.remove(node)
        else:
            for child in node._children:
                self._extract_helper(clone, child)

    ####################################################################

    def __del__(self):
        for node in self._nodes: node._unlink()
        
    ####################################################################

    def _parse(self, string, cur):

        # This parser expects a string starting with an open round
        # bracket and ending with a closed one.

        if not len(string): raise IOError, 'invalid newick string (empty string)'

        string= string[1:-1]
        i=0

        while(True):

            # case of a subtree
            subtree = []
            if string[i] == '(':

                # we gather everything until the matching closing bracket
                subtree.append(string[i])
                acc = 1
                while acc:
                    i+=1
                    if i == len(string): raise ValueError, 'invalid newick string'
                    subtree.append(string[i])
                    if string[i] == '(': acc+=1
                    if string[i] == ')': acc-=1
                i+=1

            # get all until the next comma or end of string
            buff = ''
            while i < len(string) and string[i] != ',':
                buff += string[i]
                i+=1

            # convert to label and brlen
            if not len(buff):
                label = None
                brlen = None
            else:
                tbuff = buff.split(':')
                if len(tbuff) == 1:
                    label = buff
                    brlen = None
                elif len(tbuff) == 2:
                    if len(tbuff[0]) == 0: label = None
                    else: label = tbuff[0]
                    try: brlen = float(tbuff[1])
                    except ValueError: raise ValueError, 'invalid newick string (invalid branch length: {0}'.format(tbuff[1])
                else: raise ValueError, 'invalid newick string'

            # add the node (whichever internal or terminal)
            node = cur._add_child(label=label, brlen=brlen)
            self._nodes.append(node)

            # in case there is a subtree, it is recursively parsed
            if len(subtree):
                self._parse(''.join(subtree), node)
                # try to convert label to int of float (silent otherwise)
                if node._label != None:
                    try: node._label = int(node._label)
                    except ValueError:
                        try: node._label = float(node._label)
                        except ValueError: pass

            # else, cache the leaf
            else:
                self._leaves.append(node)

            # if we didn't reach the end of the string, we continue parsing
            if i == len(string): break
            else:
                i+=1
                continue

    ####################################################################

    def __str__(self):
        return self._base.newick()

    ####################################################################

    @property
    def num_nodes(self):

        """
        Total number of nodes in the tree (the number is never smaller
        than 1, even for empty trees).
        """

        return len(self._nodes)

    ####################################################################

    @property
    def num_leaves(self):

        """
        Number of terminal nodes in the tree. If the tree is empty (only
        the default base node), the number of leaves is 0.
        """

        return len(self._leaves)

    ####################################################################

    def newick(self, labels=True, brlens=True):

        """
        Return the newick-formatted string representing the instance.

        :param labels: if ``False``, omit the internal branch labels.
        :param brlens: if ``False``, omit the branch lengths.
        """

        return self._base.newick(labels, brlens)

    ####################################################################

    @property
    def base(self):

        """
        Basal node of the tree (if the tree is unrooted, it is a
        trifurcation whose location should be considered as arbitrary;
        if the tree is rooted, it is the root). This attribute is a
        :class:`~.Node` instance which can be modified, but it cannot be
        replaced.
        """

        return self._base

    ####################################################################

    def add_node(self, parent, label=None, brlen=None):

        """
        Add a node to the tree.

        :param parent: one of the nodes of this instance, as a
            :class:`~.Node` reference.

        :param label: node label which will be internal node label or
            leaf name according to the final structure of the tree. The
            new node has initially no children and is therefore a leaf
            until it is itself connected to a child (if ever).

        :param brlen: length of the branche connecting *parent* to the
            new node.

        :return: The new node as a :class:`~.Node` reference.
        """

        if parent not in self._nodes: raise ValueError, 'cannot add node to a node that does not belong to this tree'
        child = parent._add_child(label, brlen)
        self._leaves.append(child)
        self._nodes.append(child)
        if parent in self._leaves:
            self._leaves.remove(parent)

    ####################################################################

    def iter_leaves(self):

        """
        Return an iterator over the leaf nodes.

        .. versionadded:: 3.0.0
            Replaces the method :meth:`get_terminal_nodes` (the user
            needs to get labels of those to replace the method
            :meth:`all_leaves`).
        """

        for node in self._leaves: yield node

    ####################################################################

    def get_leaf(self, label):

        """
        Return the terminal node (as a :class:`~.Node` instance) that
        has the requested leaf label. If several nodes have this label,
        returns the first one. If no nodes have this label, returns
        ``None``.
        """

        for leaf in self._leaves:
            if leaf._label == label: return leaf
        else:
            return None

    ####################################################################

    def depth_iter(self, start=None):

        """
        Return a depth-first iterator. Iterate over the :class:`~.Node`
        instances of the trees, starting from the base but, then,
        following a depth-first order.

        :param start: start point of the iteration, as a :class:`.Node`
            instance of this tree. By default, start from the base of
            the tree.
        """

        if start is None: start = self._base
        elif start not in self._nodes: raise ValueError, 'node does not belong to this tree'
        return _tree_depth_first_iterator(self, start)

    ####################################################################

    def breadth_iter(self, start=None):

        """
        Return a breadth-first iterator. Iterate over the
        :class:`~.Node` instances of the trees, starting from the base
        but, then, following a breadth-first order.

        :param start: start point of the iteration, as a :class:`.Node`
            instance of this tree. By default, start from the base of
            the tree.
        """

        if start is None: start = self._base
        elif start not in self._nodes: raise ValueError, 'node does not belong to this tree'
        return _tree_breath_first_iterator(self, start)

    ####################################################################

    def total_length(self):

        """
        Compute the sum of all branch lengths across all nodes of the
        tree. All branch lengths must be defined (non-``None``),
        otherwise a :class:`ValueError` will be raised.
        """

        L = 0
        for node in self._nodes:
            for child in node.children():
                l = node.branch_to(child)
                if l == None:
                    raise ValueError, 'cannot compute tree\'s length: at least one branch has no defined length'
                L += l
        return L

    ####################################################################

    def find_clade(self, names, ancestral=False, both_sides=False):

        """
        Check whether a group is one of the clades defined by the tree.

        The leaf names must be provided as a list, or any other iterable
        (ideally, a :class:`set`). Leaf names are normally :class:`str`
        instances. All leaves must be present in the tree.

        With the option ``ancestral=False``, search for the clade that
        contains the provided list of names as descendant. There must
        not be any other name amongst descendant. If the tree is
        unrooted and the representation uses a base within the clade, it
        will not be detected. It possible (with the option
        ``both_sides=True`` to allow searching for the complement of the
        clade, thereby detecting the right clade even if it is at the
        root. By default, if this situation occurs, the clade will not
        be detected. Return ``None`` if the clade is not found.

        With the option ``ancestral=True``, search for the most recent
        common ancestor of all leaves specified in *names*. Use of this
        method necessarily supposes that the tree is rooted (however,
        there is no requirements regarding its shape such as bifurcation
        at the base) and it is not allowed to specify
        ``both_sides=True``. With this option, it is not possible to
        have a ``None`` return value (since it is required that all
        leaves are present in the tree, in the worse case the base of
        the tree is returned).

        :param name: a :class:`set` (or compatible) specifying the
            requested leaves (as node labels, normally :class:`str`
            instances).

        :param ancestral: whether to look for the most recent common
            ancestral clade containing requested leaves (by default,
            looking for the clade containing the exact same list of
            leaves).

        :param both_sides: only allowed with ``ancestral=False``. Look
            for both the requested list of leaves and its complement,
            allowing to detect a clade even if it is spanning the base
            of the tree.

        :return: The :class:`~.Node` instance, if it exists, which has
            the exact same list of descendants than *taxa*. If no such
            clade is found, returns ``None``.

        .. warning::
            This method assumes that all leaf names of the tree are
            unique, as well as the list of names provided as argument.
            If this condition is not fulfilled, the right clade might
            not be found even if it exists.

        .. versionchanged:: 3.0.0
            Replaces previous methods :meth:`findGroup`,
            :meth:`findMonophyleticGroup`, :meth:`smallest_group` and
            :meth:`smallest_monophyleticGroup` with a modification of
            the underlying algorithm.
        """

        query = set(names)

        # capture errors or trivial cases
        if ancestral == True and both_sides == True: raise ValueError, 'cannot combine `ancestral` and `both_sides` options'
        if len(names) == 0: raise ValueError, 'cannot find clade: empty list of names'
        leaves_set = set([leaf._label for leaf in self._leaves])

        if query < leaves_set: pass
        elif query == leaves_set: return self._base # clade is the tree itself
        else: raise ValueError, 'cannot find clade: names are not all present in the tree'
        if both_sides: compl = leaves_set - query

        # find any leaf belonging to query
        for node in self._leaves:
            if node._label in query:
                query.discard(node._label)
                break

        # at least one of the names should be found
        else: raise RuntimeError, 'unexpected case in `Tree.find_clade()`: please report is as a bug'

        # climbs up the tree collecting other leaves from the query
        while node != None:

            # the descendant of the current node must have been all processed
            # they are all found to be part of the query (and removed from it)
            # if query completed, it is the right node
            if len(query) == 0: return node

            # there should be a parent because if node is the root, then we should have collected all queried leaves
            if node.parent == None: raise RuntimeError, 'unexpected case in `Tree.find_clade()`: please report is as a bug'

            # collect descendants of the parent (except those already processed through the current node)
            siblings = node.siblings()

            # move one level up
            node = node.parent

            for sib in siblings:
                for label in sib.leaves_down():

                    # discard leaves that are in query
                    if label in query:
                        query.discard(label)

                    # if non-query leaf found
                    else:

                        # if we are looking for the most recent common ancestor, it is fine
                        if ancestral == True: pass

                        # if the current node is the base, the clade might be spanning the base
                        # try the complement (if user allowed it)
                        elif node.parent == None and both_sides == True:
                            return self.find_clade(names=compl, both_sides=False)

                        # otherwise there is no such clade in this tree
                        else: return None

        # in principle, we never reach this point if the query
        raise RuntimeError, 'unexpected case in `Tree.find_clade()`: please report is as a bug'

    ####################################################################

    def collapse(self, node, ignore_len=False, ignore_label=False):

        # note: this API method is also used by remove_node() with explicit option values

        """
        Collapse a branch of the tree.

        :param node: :class:`~.Node` representing the branch to remove.
        :param ignore_len: don't try to transfer branch lengths to
            children.
        :param ignore_label: don't transfer label of the destroyed node
            to its parent.

        *node* represents the branch that must be removed from the tree
        (this node is destroyed in the process). It must be one of the
        nodes contained in the tree (as a :class:`~.Node` instance), but
        not the base of the tree. It cannot be an terminal node (leaf). 

        If *ignore_label* is not set to ``True``, the label of the
        destroyed node is transferred to the parent based on the
        following procedure: (1) if the destroyed node's label is
        ``None``, nothing is done; (2) if the destroyed node's parent's
        label is ``None``, the destroyed node's label is copied to its
        parent as is; (3) otherwise both labels are converted to strings
        (if they are not yet) and concatenated as in the string ``a;b``
        (where ``a`` is the parent's label and ``b`` is the destroyed
        node's label), even if the two labels are identical.

        If *ignore_len* is not set to ``True``, and if the length of the
        removed branch (branch from the specified node to its parent) is
        specified, it will be spread equally among the branches to its
        children (see example below). This requires that the branch
        length to all children are specified. If the removed branch has
        no specified length, nothing is done.

        Collapsing node ``[4]`` on the following tree::
        
             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]             /-------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \-------------->[6]
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------->[11]
                                          
        will generate the following tree, with the correction of edge
        lengths as  depicted::

             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]
             |             |
             |             |-------------------->[5]        L5 = L5+L4/2
            [0]            |
             |             \-------------------->[6]        L6 = L6+L4/2
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------->[11]

        Although the total edge length of the tree is not modified, the
        relationships will be altered: the distance between the
        descendants of the collapsed node (nodes 5 and 6 in the example
        above) will be artificially increased.
        """

        # sanity checks
        if node not in self._nodes != self: raise ValueError, 'cannot collapse a node that does not belong to this tree'
        if node == self._base: raise ValueError, 'cannot collapse the base of the tree'
        if node.num_children == 0: raise ValueError, 'cannot collapse a terminal node'
        parent = node.parent
        children = list(node.children())

        # process branch lengths
        p_len = node.parent_branch
        s_len = [node.branch_to(child) for child in children]
        if not ignore_len and p_len != None:
            if None in s_len: raise ValueError, 'cannot collapse node: branches to children are required'
            s_len = [length + float(p_len)/node.num_children for length in s_len]

        # disconnects the node to be collapsed and remove it
        node._remove_parent()
        for child in children: node._remove_child(child)
        self._nodes.remove(node)

        # connects the sons their grandparent
        for child, length in zip(children, s_len):
            parent._connect(child, length)

        # saves the label information
        if not ignore_label:
            plab = parent._label
            nlab = node._label
            if nlab != None:
                if plab == None: parent.label = nlab
                else: parent.label = '%s;%s' %(plab, nlab)

    ####################################################################

    def root(self, outgroup, branch_split=0.5, reoriente=False):

        """
        Roots or reoriente the tree. By default, a new node is created
        to represent the root and is placed on the branch leading to
        the provided outgroup node (the second argument determines where
        the new node is placed on this branch). Otherwise, the tree is
        reoriented such as its base is placed at the location of the
        provided outgroup. In the former case, its ends with a
        bifurcation at the root; in the latter case, a trifurcation.

        :param outgroup: :class:`~.Node` instance contained in this
            tree. It can be a leaf or any internal node, but not the
            current base of the tree (unless *reoriente* is ``True``: in
            that case, it *might* be the base of the tree [it will not
            change anything] and it *cannot* be a leaf).

        :param branch_split: where to cut the branch leading to the
            outgroup.

        :param reoriente: don't create any root node (*branch_site* is
            therefore not considered) and only place the node provided
            as the *outgroup* argument at the base of the tree, thereby
            merely changing the representation of the tree.

        The information below describes the case where
        ``reoriented=False`` (proper rooting).

        If this branch to the provided outgroup doesn't have a branch
        length, the *branch_split* argument is ignored. Otherwise,
        *branch_split* must be a real number between 0 and 1 and gives
        the proportion of the branch that must be allocated to the basal
        branch leading to the outgroup, the complement being allocated
        to the branch leading to the rest of the tree. If *branch_split*
        is either 0 or 1, one of the branch will have a length of 0, but
        it will exist anyway.

        It is illegal to call this method on trees that are already
        rooted (have a difurcation at the root).

        If the original tree has this structure::

             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]             /-------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \-------------->[6]
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \---[ROOT]-->[9]
                                          |
                                          \------------->[11]

        And rooting is requested at node ``[9]``, the root will be
        placed on the branch marked by ``[ROOT]``.  The outcome will be
        as depicted below, with the introduction of a new node (marked
        ``[ROOT]``) and the reorientation of the tree to place it at the
        base::

                                /-------------------------[1]
                                |
                         /-----[0]     /------------------[3]
                         |      |      |
                         |      \-----[2]      /----------[5]
                         |             |       |
             /---[E2]---[7]            \------[4]
             |           |                     |
             |           |                     \----------[6]
             |           |
           [ROOT]        \--------------------------------[8]
             |
             |                     /---------------------[10]
             |                     |
             \--------[E1]--------[9]
                                   |
                                   \---------------------[11]

        In this example, the relationship between nodes ``[7]`` and
        ``[0]`` (the previous base of the tree) is reverted. The label
        of node ``[7]`` is automatically transferred to node ``[0]``.
        This is consistent with the idea that internal node labels
        describe a property of the branch. The original label of the
        base, if it exists, is discarded. Since the branch between
        ``[7]`` and ``[9]`` is cut in two, the original label of node
        ``[9]`` is copied to node ``[7]``, leaving them both with the
        same label. However, if the outgroup is a terminal node, the
        label is not copied and the other basal branch is left without
        label.

        Let ``L`` be the length of the branch from ``[7]`` and ``[9]``
        in the original tree, and ``r`` the value of the parameter
        ``branch_split``. The length of the branch ``[E1]`` will be set
        to ``rL``, and the branch ``[E2]`` to ``(1-r)L``. Overall, the
        length of the tree will not be modified.

        In the case that ``reoriente=True``, the final tree is rather::

            /------------------------------------------->[10]
            |
            |------------------------------------------->[11]
            |n
           [9]       /----------------------------------->[8]
            |        |
            \------>[7]       /-------------------------->[1]
                     |        |
                     \------>[0]      /------------------>[3]
                              |       |
                              \----->[2]        /-------->[5]
                                      |         |
                                      \------->[4]
                                                |
                                                \-------->[6]

        The topology of the tree is the same as the initial one, except
        that the base is now ``[9]``. The lengths of all branches are
        conserved. However, node labels between the old and the new base
        are reverted: the node label of the new base (``[9]`` in the
        example) is affected to the next node (``[7]`` in the example)
        and so on until the old base (``[0]`` in the example), whose
        labe, if it exists, is discarded.

        .. versionchanged:: 3.0.0
            Merged with :meth:`reoriente`.
        """

        # sanity checking
        if outgroup not in self._nodes: raise ValueError, 'cannot root tree: outgroup node is not part of the tree'
        if outgroup == self._base and reoriente == False: raise ValueError, 'cannot root tree: outgroup node is the current base of the tree'
        if outgroup.num_children == 0 and reoriente == True: raise ValueError, 'cannot reoriente tree: outgroup node is a terminal node'
        if self._base.num_children == 2: raise ValueError, 'cannot root tree: tree is already rooted'
        if self._base.num_children < 2: raise ValueError, 'cannot root tree: tree has an invalid structure'
        L = outgroup.parent_branch
        if not reoriente and L != None and (branch_split < 0.0 or branch_split > 1.0): raise 'cannot root tree: invalid `branch_split` value'

        # collect all nodes on the path from the new to the old root
        path = []
        labels = [] # collect also labels
        brlens = [] # and branch lengths

        if reoriente: cur = outgroup
        else: cur = outgroup.parent

        while cur != None:
            path.append(cur)
            labels.append(cur._label) # last one will be discarded
            brlens.append(cur._parent_branch) # last one is None
            cur = cur.parent
            # path[0] = node after new root (or new base if reoriente)
            # path[...] = intermediate nodes
            # path[-1] = base of the tree

        if not reoriente:

            # remove the branch where the root is placed
            outgroup._remove_parent()

            # create root and add it to nodes
            root = Node()
            self._nodes.append(root)
            self._base = root

        else:
            self._base = outgroup

        # disconnect all nodes on the way from new base to old base (old base not included)
        for node in path[:-1]: node._remove_parent()

        # reconnect in the other order
        for i in xrange(len(path)-1):
            path[i]._connect(path[i+1], brlens[i])  # note that brlens[-1] is discarded / is None anyway
            path[i+1]._label = labels[i]            # labels[-1] is discarded and path[0]._label unset yet

        if not reoriente:

            # compute length of basal branches
            if L != None:
                L1 = L * branch_split
                L2 = L * (1.0 - branch_split)
            else:
                L1 = None
                L2 = None

            # connect the new root to the two subtrees
            root._connect(outgroup, L1)
            root._connect(path[0], L2)

        if not reoriente:
            # copy the root branch label to the new branch
            if outgroup.num_children > 0: path[0]._label = outgroup._label
            else: path[0]._label = None
        else:
            outgroup._label = None

    ####################################################################

    def map_descendants(self):

        """
        Generate and return a :class:`dict` where keys are all nodes (as
        :class:`~.Node` instances) of the tree excepting its base and
        all terminal node, and values are the :class:`tuple` of leaf
        labels that descend from this node.

        .. versionadded:: 3.0.0
        """

        mapping = {}
        queue = {}

        # add the parent of all leaves to the queue
        for leaf in self._leaves:
            if leaf._parent not in queue: queue[leaf._parent] = []
            queue[leaf._parent].append([leaf._label])

        # process queue'd node and their parent until exhaustion
        while len(queue):
            to_del = set()
            to_add = {}

            for node, items in queue.iteritems():

                # this checks if the node has been fully processed
                if len(items) == len(node._children):

                    # if the node is the base of the tree, just trash the labels
                    if node._parent != None:

                        # otherwise, send the data to the mapping...
                        mapping[node] = reduce(list.__add__, items)

                        # ... and propagate labels to higher level
                        if node._parent in queue: queue[node._parent].append(mapping[node])
                        elif node._parent in to_add: to_add[node._parent].append(mapping[node])
                        else: to_add[node._parent] = [mapping[node]]

                    # remove the node from the queue
                    to_del.add(node)

            for node in to_del: del queue[node]
            queue.update(to_add)

        for node in mapping: mapping[node] = tuple(mapping[node])
        return mapping

    ####################################################################

    def frequency_nodes(self, trees, relative=False):
        
        """
        Labels all nodes of the current instances by integers counting
        the number of trees where the same node exists among the trees
        in the iterable *trees*. It is required that all leaf labels are
        unique.

        :param trees: an iterable containing :class:`~.Tree` instances
            with exactly the same set of leaf labels.
        :param relative: node frequencies are expressed as fractions.
            The use of this option requires that at least one tree is
            provided.
        :return: Nothing (result is available as node labels).

        .. warning::
            With the exception of the base of the tree (which is ignored
            by this function) and leaf labels, all previously set labels
            are erased.

        .. versionchanged:: 3.0.0
            Labels are not converted to strings anymore.
        """

        # prepare a set for safety checking
        leaves = [node._label for node in self._leaves]
        n = len(leaves)
        leaves = set(leaves)
        if len(leaves) != n: raise ValueError, 'cannot compute node frequencies: not all leaf labels are unique'

        # get own list of nodes and reverse relationships
        nodes = dict([(clade, node) for (node, clade) in self.map_descendants().iteritems()])

        # initialize counter dict
        counts = dict.fromkeys(nodes.itervalues(), 0)

        # process all provided trees
        n = 0
        for tree in trees:
            n += 1

            # check consistency
            if set([node._label for node in tree._leaves]) != leaves: raise ValueError, 'cannot compute node frequencies: the list of leaves differ between trees!'

            # count if any set clade of provided tree matches with a clade of the reference tree
            for clade in tree.map_descendants().itervalues():
                if clade in nodes: counts[nodes[clade]] += 1

        # compute relative frequencies
        if relative:
            if n == 0: raise ValueError, 'cannot compute relative node frequencies: no provided trees!'
            for node in counts: counts[node] = float(counts[node]) / n

        # annotate nodes
        for node, count in counts.iteritems(): node._label = count

    ####################################################################

    def clean_internal_labels(self):

        """
        Remove all internal node labels, including the base of the tree.
        In practice, they are set to ``None``.
        """

        for node in self._nodes:
            if len(node._children) > 0: node._label = None

    ####################################################################

    def clean_branch_lengths(self):

        """
        Remove all branch lengths. In practice, they are set to
        ``None``.
        """

        for node in self._nodes:
            node._parent_branch = None # also the base 
            for child in node._children: node._child_branches[child] = None

    ####################################################################

    def remove_node(self, node, drop_parent=True):

        """
        Remove a node from the tree, as well as all its descendants.
        Since this operation may create a node with a single child, this
        method may remove the parent or the brother of the removed node
        depending on the structure of the tree (see *drop_parent*),
        unless specified otherwise.

        :param node: the node to remove, as a :class:`~.Node` instance
            belowing to the current tree. Terminal nodes can be removed,
            but not the base of the tree.

        :param drop_parent: if ``True``, remove the parent of the
            removed node if it is left with only one child. If the
            parent is the base of the tree, remove the other descendant
            if it is not terminal (see below).

        Assume we remove node ``[3]`` from the tree with this
        structure::

                               /---------------------------->[2]
                               |
                 /----------->[1]           /--------------->[4]
                 |             |            |
                 |             \---------->[3]
                 |                          |
                [0]                         \--------------->[5]
                 |
                 |              /--------------------------->[7]
                 |              | 
                 \------------>[6]            /------------->[9]
                                |             |
                                \----------->[8]
                                              |
                                              \------------>[10]

        Then, we would end up with the following tree::

                 /----------->[1]--------------------------->[2]
                 |
                 |
                [0]
                 |              /--------------------------->[7]
                 |              | 
                 \------------>[6]            /------------->[9]
                                |             |
                                \----------->[8]
                                              |
                                              \------------>[10]

        The default behaviour is then to remove node ``[1]`` (and delete
        its label if it exists) and to set the length of the branch
        from ``[0]`` to ``[2]`` to the sum of the ``[0]`` to ``[1]`` and
        ``[1]`` to ``[2]``. But, with ``drop_parent=False``, the tree is
        left as is.

        There is a special case with the base of the tree. Assume that
        we remove node ``[1]`` from the original tree above. We then
        would have a non-standard structure with a single child to the
        base of the tree::

                                /--------------------------->[7]
                                | 
                [0]----------->[6]            /------------->[9]
                                |             |
                                \----------->[8]
                                              |
                                              \------------>[10]

        In that case, the base is not removed, but node ``[6]`` is
        removed using the :meth:`~.collapse` method (using options
        ``ignore_len=False`` but ``ignore_label=True`` since the base is
        not suppose to bear a label). We end up with the following
        structure::

                 /--------------------------------->[7]
                 | 
                [0]                  /------------->[9]
                 |                   |
                 \----------------->[8]
                                     |
                                     \------------>[10]

        The length of the branch from ``[0]`` to ``[6]`` is spread
        equally between the branch from ``[0]`` to ``[7]`` and the
        branch from ``[0]`` to ``[8]`` (and so on if there are actually
        more than one descendants). if ``drop_parent=False`` or if
        ``[6]`` is a terminal node, it is not removed.
        """

        # sanity check
        if node not in self._nodes: raise ValueError, 'cannot delete a node that is not part of the tree'
        if node == self._base: raise ValueError, 'cannot delete the base of the tree'

        # record parent
        parent = node._parent

        # disconnect the subtree
        node._remove_parent()

        # remove nodes of the subtree
        for node in Tree._gather_nodes(node):
            node._unlink()
            self._nodes.remove(node)
            if node in self._leaves: self._leaves.remove(node)

        # collapse the ascendant/brother if needed
        if drop_parent:
            if parent.num_children == 1:
                if parent == self._base:
                    brother = parent._children[0]
                    if brother.num_children != 0:
                        self.collapse(brother, ignore_len=False, ignore_label=True)
                else:
                    self.collapse(parent, ignore_len=False, ignore_label=True)

    ####################################################################

    @staticmethod
    def _gather_nodes(node):

        # collect all nodes from a subtree: return list of this node and
        # all its descendants

        nodes = [node]
        for child in node._children:
            nodes.extend(Tree._gather_nodes(child))
        return nodes

    ###################################################################

    def lateralize(self, reverse=False):

        """
        Modify the order of children of all node of the trees in such a
        way such that the children are sorted from the smallest to the
        largest number of descending leaves.

        :param reverse: sort from in the more-descendants to
            less-descendants order instead.
        """

        # get node mapping (only store number of descendants)
        nodes = self.map_descendants()
        for node in nodes: nodes[node] = len(nodes[node])

        # add leaves
        for leaf in self._leaves: nodes[leaf] = 1

        # sort all nodes of the tree
        for node in self._nodes:
            numdict = dict(zip(node._children, [nodes[child] for child in node._children]))
            node._children.sort(key=numdict.get, reverse=reverse)

    ####################################################################

    def midroot(self):

        """
        Automatic midpoint rooting of the tree. The tree must be
        initially unrooted (trifurcation at the root). This method
        identifies the most distant pair of terminal nodes (in case of
        a draw, one is picked randomly) and the root of the tree (as a
        new node) placed at the middle point of this path.
        """

        # sanity check
        if len(self._nodes) < 3: raise ValueError, 'cannot perform automatic rooting: not enough branches in tree'
        if self._base.num_children < 3: raise ValueError, 'cannot perform automatic rooting: tree must have a trifurcation at the root'

        # collect the path to root for all leaves
        paths = []
        for node in self._leaves:
            path = []
            while node.parent != None:
                L = node.parent_branch
                if L == None: raise ValueError, 'cannot perform automatic rooting: all branch lengths must be specified'
                path.append((node, L))
                node = node.parent
            paths.append(path)

        # find shortest path(s)
        best = []
        best_dist = - sys.float_info.max

        for i in xrange(len(paths)):
            for j in xrange(i+1, len(paths)):

                part1 = list(paths[i]) # deep copies to allow editing
                part2 = list(paths[j])
                while part1[-1] == part2[-1]: # the first of each at least is different
                    del part1[-1]
                    del part2[-1]
                path = [(node, L, False) for (node, L) in part1] + [(node, L, True) for (node, L) in part2[::-1]]
                d = sum([L for (node, L, reverse) in path])

                if (d - best_dist) > -0.000000001:
                    if (d - best_dist) > +0.000000001:
                        best_dist = d
                        best = []
                    best.append(path)

        # if more than several best, pick one randomly
        if len(best) > 1:
            rnd = _random_pool.get()
            best = best[rnd.irand(len(best))]
            _random_pool.put(rnd)
        else: best = best[0]

        # find midpoint
        lim = best_dist / 2.0
        acc = 0.0
        for node, L, reverse in best:
            acc += L
            if acc - lim > -0.000000001:
                outgroup = node
                if reverse: offset = min(1.0, (acc - lim) / L)
                else: offset = max(0.0, (acc - lim) / L)
                break
        else:
            raise RuntimeError, 'midpoint rooting: cannot find root location (please report bug)'

        # root itself
        self.root(outgroup, offset)

    ####################################################################

    def unroot(self, reverse=False):

        """
        Remove the root. The tree must be initially unrooted
        (bifurcation at the root). This method removes the root node and
        places the base of the tree at one of the two basal nodes (the
        nodes that are ancestral to the two basal groups). This method
        does not change to total length of the tree. And error is raised
        if only one of the two basal branches has a length. If the
        initial basal node has a label, it is lost. If the node that
        becomes the base has a label, it is left there (it will appear
        a the base of the tree).

        :param reverse: if ``True``, place the base of the tree at the
            second basal node (by default, the first basal node is
            used).
        """

        if self._base.num_children != 2: raise ValueError, 'cannot unroot tree: tree is not rooted'
        node1 = self._base.child(0)
        node2 = self._base.child(1)
        L1 = node1.parent_branch
        L2 = node2.parent_branch
        if L1 == None and L2 == None: L = None
        elif L1 != None and L2 != None: L = L1 + L2
        else: raise ValueError, 'cannot unroot tree: only one of the two basal branches has a length'
        self._base._remove_child(node1)
        self._base._remove_child(node2)
        self._nodes.remove(self._base)
        if reverse:
            node2._connect(node1, L)
            self._base = node2
        else:
            node1._connect(node2, L)
            self._base = node1

########################################################################

class _tree_breath_first_iterator(object):

    def __init__(self, tree, start):
        self._this_level = [start]
        self._next_level = list(tree._base._children) # must be a deep copy

    def __iter__(self):
        return self

    def next(self):
        if len(self._this_level) == 0:
            if len(self._next_level) == 0:
                raise StopIteration
            else:
                self._this_level = self._next_level
                self._next_level = reduce(list.__add__, [node._children for node in self._this_level])
        node = self._this_level.pop(0)
        return node

########################################################################

class _tree_depth_first_iterator(object):

    def __init__(self, tree, start):
        self._cur = start
        self._pile = []

    def __iter__(self):
        return self

    def next(self):
        if self._cur == None: raise StopIteration
        node = self._cur
        children = list(node.children())
        if len(children) == 0:
            if len(self._pile) == 0: self._cur = None
            else: self._cur = self._pile.pop()
        else:
            self._cur = children[0]
            self._pile.extend(children[1:])
        return node

########################################################################
