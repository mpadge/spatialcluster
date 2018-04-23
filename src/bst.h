#pragma once

#include <iostream>
#include <cstdlib>
using namespace std;

// Templating with recursive pointers is much harder than simply changing the
// typedef, and makes the code much less readable
typedef double data_type;

struct tree_node
{
    tree_node * lo;
    tree_node * hi;
    tree_node * parent;
    data_type data;
};

class BinarySearchTree
{
    private:
        tree_node * root;
        data_type tmin (tree_node * node);
        void clear_node (tree_node * node);
        tree_node * removeNode (tree_node * node, data_type value);

    public:
        BinarySearchTree ()
        {
            root = nullptr;
        }
        ~BinarySearchTree ()
        {
            treeClear ();
        }
        void insert (data_type);
        void remove (data_type value);
        data_type treeMin ();

        tree_node * getRoot ();
        tree_node * getNode (tree_node * node, data_type value);
        tree_node * treeMinTree ();
        tree_node * tminTree (tree_node * node);

        tree_node * nextHi (tree_node * node);

        void treeClear ();
};

void BinarySearchTree::insert (data_type d)
{
    tree_node * t = new tree_node;
    tree_node * parent;
    t->data = d;
    t->lo = nullptr;
    t->hi = nullptr;
    parent = nullptr;

    if (root == nullptr)
        root = t;
    else
    {
        tree_node * node;
        node = root;
        while (node != nullptr)
        {
            parent = node;
            if (t->data > node->data)
                node = node->hi;
            else
                node = node->lo;
        }

        t->parent = parent;
        if (t->data < parent->data)
            parent->lo = t;
        else
            parent->hi = t;
    }
}

void BinarySearchTree::remove (data_type value)
{
    root = removeNode (root, value);
}

// recursive private member function:
tree_node * BinarySearchTree::removeNode (tree_node * node, data_type value)
{
    if (node == nullptr)
        return node;

    if (value < node->data) {
        node->lo = removeNode (node->lo, value);
    } else if (value > node->data) {
        node->hi = removeNode (node->hi, value);
    } else {
        if (node->lo == nullptr && node->hi == nullptr) { // no children
            delete node;
            node = nullptr;
        }
        else if (node->lo == nullptr) { // 1 child: hi
            tree_node * temp = node;
            node->hi->parent = node->parent;
            node = node->hi;
            delete temp;
        }
        else if (node->hi == nullptr) { // 1 childe: lo
            tree_node * temp = node;
            node->lo->parent = node->parent;
            node = node->lo;
            delete temp;
        }
        else // 2 children
        {
            tree_node * temp = tminTree (node->hi);
            node->data = temp->data;
            node->hi = removeNode (node->hi, temp->data);
        }
    }
    return node; // then the root node which needs to be updated
}

data_type BinarySearchTree::treeMin ()
{
    return tmin (root);
}

data_type BinarySearchTree::tmin (tree_node * node)
{
    while (node->lo != nullptr)
        node = node->lo;

    return node->data;
}

tree_node * BinarySearchTree::treeMinTree ()
{
    return tminTree (root);
}

tree_node * BinarySearchTree::tminTree (tree_node * node)
{
    while (node->lo != nullptr)
        node = node->lo;

    return node;
}

void BinarySearchTree::treeClear ()
{
    clear_node (root);
}

void BinarySearchTree::clear_node (tree_node * node)
{
    if (node != nullptr)
    {
        clear_node (node->lo);
        clear_node (node->hi);
        delete node;
    }
}

tree_node * BinarySearchTree::getRoot ()
{
    tree_node * node = root;
    return node;
}

tree_node * BinarySearchTree::getNode (tree_node * node, data_type value)
{
    if (node == nullptr)
    {
        //std::cout << "value = " << value <<
        //    " does not exist in the tree" << std::endl;
        return node;
    }

	if (value < node->data) 
		return getNode (node->lo, value);
	else if (value > node->data)
		return getNode (node->hi, value);
	else
		return node;
}

tree_node * BinarySearchTree::nextHi (tree_node * node)
{
    if (node->hi != nullptr)
        return tminTree (node->hi);

    tree_node * y = node->parent;
    while (y != nullptr && node == y->hi)
    {
        node = y;
        y = y->parent;
    }
    //if (y == nullptr)
    //    std::cout << "already at max value" << std::endl;
    return y;
}
