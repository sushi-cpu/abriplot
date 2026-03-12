"""
abriplot.newick
===============
Pure-Python Newick tree parser.  No Biopython required.
"""

import numpy as np


class TreeNode:
    __slots__ = ("name", "length", "children", "parent", "x", "y", "leaf_index")

    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.x = 0.0
        self.y = 0.0
        self.leaf_index = -1

    def is_leaf(self):
        return len(self.children) == 0

    def leaves(self):
        if self.is_leaf():
            yield self
        for c in self.children:
            yield from c.leaves()


def _tokenise(s):
    i, n = 0, len(s)
    while i < n:
        ch = s[i]
        if ch in "(),;":
            yield ch
            i += 1
        elif ch == ":":
            i += 1
            j = i
            while j < n and s[j] not in "(),;":
                j += 1
            yield ("BL", s[i:j].strip())
            i = j
        elif ch in " \t\n\r":
            i += 1
        else:
            j = i
            while j < n and s[j] not in "():,;":
                j += 1
            yield ("LB", s[i:j].strip())
            i = j


def parse_newick(s):
    """Parse a Newick string and return the root TreeNode."""
    s = s.strip().rstrip(";")
    tokens = list(_tokenise(s))
    root = TreeNode()
    cur = root
    for tok in tokens:
        if tok == "(":
            child = TreeNode()
            child.parent = cur
            cur.children.append(child)
            cur = child
        elif tok == ")":
            cur = cur.parent
        elif tok == ",":
            sib = TreeNode()
            sib.parent = cur.parent
            cur.parent.children.append(sib)
            cur = sib
        elif isinstance(tok, tuple) and tok[0] == "LB":
            cur.name = tok[1]
        elif isinstance(tok, tuple) and tok[0] == "BL":
            try:
                cur.length = float(tok[1])
            except ValueError:
                cur.length = 0.0
    return root


def assign_x(node, cum=0.0):
    """Assign cumulative x (branch-length) coordinates."""
    node.x = cum + node.length
    for c in node.children:
        assign_x(c, node.x)


def assign_y(node):
    """Assign y coordinates: leaves get integer indices, internals get midpoints."""
    leaves = list(node.leaves())
    for idx, lf in enumerate(leaves):
        lf.y = idx
        lf.leaf_index = idx

    def _iy(n):
        if n.is_leaf():
            return
        for c in n.children:
            _iy(c)
        n.y = np.mean([c.y for c in n.children])

    _iy(node)
    return leaves


def tip_order(tree_string, strip_trailing_underscore=True):
    """
    Parse a Newick string and return:
        names  : ordered list of tip label strings
        root   : root TreeNode
        leaves : list of leaf TreeNodes in display order
    """
    root = parse_newick(tree_string)
    assign_x(root)
    leaves = assign_y(root)
    names = [lf.name for lf in leaves]
    if strip_trailing_underscore:
        names = [n.rstrip("_") for n in names]
    return names, root, leaves
