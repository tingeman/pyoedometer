from PyQt5 import QtGui, QtCore, QtWidgets


class TreeItem(object):
    """ A python object used to return row/column data, and keep note of
    it's parents and/or children in a Tree structure
    """

    def __init__(self, data, parent=None):
        self.item_data = data
        self.parentItem = parent
        self.childItems = []

    def appendChild(self, item):
        self.childItems.append(item)

    def child(self, row):
        return self.childItems[row]

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return len(self.item_data)

    def data(self, column):
        try:
            return self.item_data[column]
        except IndexError:
            return None

    def parent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


class TreeModel(QtCore.QAbstractItemModel):
    """ A Tree model to display a few names, ordered by sex
    """

    def __init__(self, headers=[], parent=None):
        """ Initialize the tree model and add data """
        super(TreeModel, self).__init__(parent)

        self.headers = headers
        self.rootItem = TreeItem(headers, None)
        #self.rootItem = QtGui.QStandardItem(headers,)
        self.parents = {0: self.rootItem}
        self.setupModelData()

    def columnCount(self, parent=None):
        """ Return the number of colums or fields in the tree model """
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        """ Adapted method to retrieve the data of a specific index in the model
        NOT SURE EXACTLY WHAT THIS DOES!
        ROLE determines whether the item is viewable or editable etc.
        """
        if not index.isValid():
            return None

        if role != QtCore.Qt.DisplayRole:
            return None

        item = index.internalPointer()
        return item.data(index.column())

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.NoItemFlags

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self, column, orientation, role):
        """ STANDARD METHOD to retrieve the header of a column """
        if (orientation == QtCore.Qt.Horizontal and
                    role == QtCore.Qt.DisplayRole):
            try:
                return self.headers[column]
            except IndexError:
                pass

        return None

    def index(self, row, column, parent):
        """ STANDARD METHOD to retrieve the index of an item """
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        """ STANDARD METHOD to retrieve the parent of an item with the given index """
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        if not childItem:
            return QtCore.QModelIndex()

        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QtCore.QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent=QtCore.QModelIndex()):
        """ STANDARD METHOD to retrieve number of rows in the model """
        if parent.column() > 0:
            return 0
        if not parent.isValid():
            p_Item = self.rootItem
        else:
            p_Item = parent.internalPointer()
        return p_Item.childCount()

    def addItem(self, data, parentItem=None):
        if parentItem is None:
            parentItem = self.rootItem

        try:
            item_data = [data[key] for key in self.headers]
        except:
            if len(data) == len(self.headers):
                item_data = data
            else:
                raise ValueError('The data passed does not have the right number of elements!')

        child = TreeItem(item_data, parentItem)
        parentItem.appendChild(child)

        return child

    def setupModelData(self, data=None, parent=None):
        """ CUSTOM METHOD to setup the tree model of the supplied data """
        pass

    def searchModel(self, key=None, value=None):
        """ CUSTOM METHOD to search the tree and return the index of the specified person """

        def searchNode(node):
            """
            a function called recursively, looking at all nodes beneath node
            """
            for child in node.childItems:
                if key in child.item_data:
                    index = self.createIndex(child.row(), 0, child)
                    return index

                if child.childCount() > 0:
                    result = searchNode(child)
                    return result
                else:
                    return None

        retarg = searchNode(self.rootItem)
        return retarg


class CheckableTreeModel(TreeModel):
    """ A Tree model to display a few names, ordered by sex
    """

    def addItem(self, data, parentItem=None):
        child = super(CheckableTreeModel, self).addItem(data, parentItem)
        if parentItem is None:
            parentItem = self.rootItem

        parentItem.setFlags(parentItem.flags() | QtCore.Qt.ItemIsTristate | QtCore.Qt.ItemIsUserCheckable)
        child.setFlags(child.flags() | QtCore.Qt.ItemIsUserCheckable)
        child.setCheckState(0, QtCore.Qt.Unchecked)

        return child